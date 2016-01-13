function [tau, E, Dtau, objMaskT] = lern2frmtau(src, ref, para, grpParam)
% Yuanyi Xue @ NYU-Poly
% This function solve the following problem using ADMM:
% argmin_{E,\tau} |E|
% s.t. D\circ \tau + E = C
% D is the current frame, C is the previous frame, and \tau is the global
% warping parameter.

% Inputs:
% 1) src: source frame, grayscale (C), range from [-0.5 0.5]
% 2) ref: reference frame (D), must be the same size/range as src
% 3) para: a struct contains the following fields:
%  - tauModel: define the parameterized model using for global warping, can be {'affine', 'projective'}
%  - tau0: initial \tau, if all zero, the internal SIFT feature based method is used.
%  - numPts: number of SIFT feature points to derive initial tau0, not used if tau0 is not all zero
%  - tolInner/maxInnerIter: inner loop convergence check and maximum iterations
%  - tolOuter/maxOuterIter: outer loop convergence check and maximum iterations
%  - lambda: ADMM dual parameter 
%  - objMask: a binary mask of same frame size, where =1 indicates an object pixel
% 4) grpParam: for group sparsity only, optional
%  - graph: struct containing the description of graph, check SPAM doc
%  - param: struct containing the parameters for solving the L1\inf problem

% Outputs:
% 1) tau: learned transform
% 2) E: sparse residual frame
% 3) Dtau: transformed frame (from ref)

% Pre-requisites: RASL toolbox, SPAM, and VLfeat (if using internal initialization)

% Regular parameters
tauModel = para.tauModel;
tau0 = para.tau0;
numPts = para.numPts;
tolInner = para.tolInner;
maxInnerIter = para.maxInnerIter;
tolOuter = para.tolOuter;
maxOuterIter = para.maxOuterIter;
tolDeltaXi = para.tolDeltaXi;
lambda = para.lambda;
TVmode = para.TVmode; % Minimize TV+L1 mode
if isfield(para,'objMask')
    objMask = para.objMask;
else
    objMask = ones(size(src)); % no mask in fact
end
if isfield(para,'verbose')
    verbose = para.verbose;
else
    verbose = 0;
end
% Group sparsity related parameters
if nargin < 4,
    grpParam = [];
end
    

imgSize = size(src);

if sum(tau0(:)) == 0,
    % Use SIFT feature points to find the initial tau
    % Transform matrix from ref -> current
    % Reference frame SIFT features
    [refSIFT,refDist] = vl_sift(single(ref));
    % Find SIFT features
    [currSIFT,currDist] = vl_sift(single(src));
    [matches,scores] = vl_ubcmatch(refDist,currDist);
    % Only use the top 10 scored matches
    [junk,scoreIdx] = sort(scores,'ascend');
    numMatch = min([length(scores) numPts]);
    matchIdx = scoreIdx(1:numMatch);
    % Get the feature center correspondences
    refMatchingPt = refSIFT(1:2,matches(1,matchIdx));
    currMatchingPt = currSIFT(1:2,matches(2,matchIdx));
    % Derive affine mapping through LS
    % Mapping reference frame to current frame
    A = [refMatchingPt' ones(length(matchIdx),1)];
    x = currMatchingPt(1,:)';
    D = currMatchingPt(2,:)';
    % affine parameters {a1 a2 a3} and {b1 b2 b3}
    a = pinv(A)*x;
    b = pinv(A)*D;
    transformations = [transpose(a)
        transpose(b)
        0 0 1];
else
    transformations = tau0;
end

% image derivatives for the reference frame
% refYx = imfilter( ref, (-fspecial('sobel')') / 8 );
% refYy = imfilter( ref,  -fspecial('sobel')   / 8 );
[refYx, refYy] = imgradientxy( ref, 'CentralDifference');
%refYy = refYy/norm(refYy);
% xdata = [1 imgSize(2)];
% ydata = [1 imgSize(1)];
xdata = [1 imgSize(2)] - ceil(imgSize(2)/2);
ydata = [1 imgSize(1)] - ceil(imgSize(1)/2);

%% Main loop
T_in = transformations;
Tfm = fliptform(maketform('projective',T_in'));
% Tfm = maketform('projective',T_in');
% Initial mask
% objMaskT = vec(imtransform(objMask, Tfm, 'nearest', 'XData', xdata, 'YData', ydata, 'UData', xdata, 'VData', ydata, 'Size', imgSize, 'fill', 1));
objMaskT = vec(objMask);
maskIdx = find(objMaskT==1);

numIterInner = 0;
numIterOuter = 0;
iterNum = 0;
isconverged = false;
prevObj = inf;   

% Jacobian indices
% u = vec(repmat(1:imgSize(2),imgSize(1),1));
% v = vec(repmat((1:imgSize(1))',1,imgSize(2)));
u = vec(repmat([xdata(1):xdata(2)],imgSize(1),1));
v = vec(repmat([ydata(1):ydata(2)]',1,imgSize(2)));

while ~isconverged,
    iterNum = iterNum + 1;
    numIterOuter = numIterOuter + 1;
    
    if verbose
        disp(['Iter ' num2str(iterNum)]);
    end
        
    % TO-DO: exclude non-existing pixels
    I   = vec(imtransform(ref, Tfm,'bicubic','XData', xdata, 'YData', ydata, 'UData', xdata, 'VData', ydata,'Size',imgSize));
    I(I<0) = 0; I(I>255) = 255; % clipping the I at every iteration
    Iu  = vec(imtransform(refYx,Tfm,'bicubic','XData', xdata, 'YData', ydata, 'UData', xdata, 'VData', ydata,'Size',imgSize));
    Iv  = vec(imtransform(refYy,Tfm,'bicubic','XData', xdata, 'YData', ydata, 'UData', xdata, 'VData', ydata,'Size',imgSize));
    % transform the object mask as well
    D   = I;
    % Now mask the D, Iu, and Iv using the transformed mask
    Dmask = D.*objMaskT;
    Iumask = Iu.*objMaskT;
    Ivmask = Iv.*objMaskT;
    
    % transformation matrix to parameters
    xi = projective_matrix_to_parameters(tauModel,T_in); 
    
    % Compute Jacobian for affine parameters
    if strcmp(tauModel,'AFFINE'),
        J = [ Iumask.*u,   Iumask.*v,   Iumask,   Ivmask.*u,   Ivmask.*v,   Ivmask ];
    elseif strcmp(tauModel,'HOMOGRAPHY'),
        T = ones(3,3);
        T(1,:) = xi(1:3);
        T(2,:) = xi(4:6);
        T(3,1:2) = xi(7:8);
        X = T(1,1)*u + T(1,2)*v + T(1,3);
        Y = T(2,1)*u + T(2,2)*v + T(2,3);
        N = T(3,1)*u + T(3,2)*v + 1 + eps;
        
        J = [ Iumask .* u ./ N, Iumask .* v ./ N, Iumask ./ N, ...
            Ivmask .* u ./ N, Ivmask .* v ./ N, Ivmask ./ N, ...
            ( -Iumask .* X .* u ./ ( N.^2 ) - Ivmask .* Y .* u ./ (N.^2) ), ...
            ( -Iumask .* X .* v ./ ( N.^2 ) - Ivmask .* Y .* v ./ (N.^2) ) ];
    else
        error('tauModel is not supported!');
    end
    
    % Temporary remedy, just drop the masked region
    srcTmp = vec(src).*objMaskT;
    srcTmp = srcTmp(maskIdx);
    DmaskTmp = Dmask(maskIdx);
    Dprime = srcTmp - DmaskTmp;
    J = J(maskIdx,:);
    
    % Inner loop based on RASL
    % Solve the following:
    % min ||S||_1 s.t. J*\delta_\tau + S = b - D (Dprime)
    % -----------------------------------------------------------------
    % -----------------------------------------------------------------
%     Dprime = vec(src).*objMaskT - Dmask;
    if isempty(grpParam) && ~TVmode,
        [S, delta_xi, numIterInnerEach] = twoFrm_inner_ialm(Dprime, J, lambda, tolInner, maxInnerIter);
    elseif TVmode,
        [S, delta_xi, numIterInnerEach] = twoFrm_inner_ialm_TV(Dprime, J, lambda, tolInner, maxInnerIter, imgSize);
    else
        [S, delta_xi, numIterInnerEach] = twoFrm_inner_ialm(Dprime, J, lambda, tolInner, maxInnerIter, grpParam);
    end
    % -----------------------------------------------------------------
    % -----------------------------------------------------------------
    numIterInner = numIterInner + numIterInnerEach ;
    
    % Parameter update
    xi = xi + delta_xi;
    T_in = parameters_to_projective_matrix(tauModel,xi);    
    % Transformed image and derivatives with respect to affine parameters
    Tfm = fliptform(maketform('projective',T_in'));
%     Tfm = maketform('projective',T_in');
    % Update the mask
%     objMaskT = vec(imtransform(objMask, Tfm, 'nearest', 'XData', xdata, 'YData', ydata, 'UData', xdata, 'VData', ydata, 'Size', imgSize, 'fill', 0));
    
    % Use the updated transform to mask the S
    if TVmode,
        curObj = tvcost(reshape(S.*objMaskT,imgSize(1),imgSize(2)),lambda) ;
    else
%         curObj = norm(S.*objMaskT,1);
        curObj = norm(S,1);
    end
    
    if verbose
        disp(['previous objective function: ' num2str(prevObj) ]);
        disp([' current objective function: ' num2str(curObj) ]);
    end
      
    % Convergence check
    if ( (prevObj - curObj < tolOuter) || iterNum >= maxOuterIter || mean(abs(delta_xi)) <= tolDeltaXi)
        isconverged = true;
        if (prevObj - curObj) < 0
            % revert back to the last iteration xi
            xi = xi - delta_xi;
            T_in = parameters_to_projective_matrix(tauModel,xi);
        end
        if ( iterNum >= maxOuterIter )
            disp('Maximum iterations reached') ;
        end
    else
        prevObj = curObj;
    end
end

if verbose
    disp(['total number of iterations: ' num2str(numIterInner) ]);
    disp(['number of outer loop: ' num2str(numIterOuter) ]);
end

%% Save data
Tfm = fliptform(maketform('projective',T_in'));
% Tfm = maketform('projective',T_in');

Dtau   = vec(imtransform(ref, Tfm,'bicubic','XData', xdata, 'YData', ydata, 'UData', xdata, 'VData', ydata,'Size',imgSize));
Dtau(Dtau>255) = 255; Dtau(Dtau<0) = 0; % clipping Dtau
Dtau = Dtau.*objMaskT;
tau = T_in;
E = vec(src).*objMaskT - Dtau; % sparse frame

% E_org = reshape(E,imgSize(1),imgSize(2));
E = reshape(E,imgSize(1),imgSize(2));
Dtau = reshape(Dtau,imgSize(1),imgSize(2));
end

% Nested TV objective
function c = tvcost(x, lambda)
% TVCOST computes the TV-functional value
difh = [x(:, 2:end), x(:, 1)] - x;
difv = [x(2:end, :); x(1, :)] - x;
c = norm(x(:),1) + lambda*sum(sqrt((difh(:).^2) + (difv(:).^2)));
end