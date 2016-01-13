function [S_dual, dt_dual, iter, Y] = twoFrm_inner_ialm_TV(Dprime, J, lambda, tol, maxIter, imgSize)
% Inner loop of solving a TV+L1 minimization algorithm
% argmin_{S, dt} ||S||_1 + \lambda * ||S||_TV  s.t. J*dt + S = Dprime

if nargin < 3
    error('Too few arguments') ;
end

if nargin < 4
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 5
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

DISPLAY_EVERY = 30 ;

% initialize
Y = Dprime;
% [U S V] = svd(Dprime);
% Y = U * V';
norm_two = norm(Y, 2);
norm_inf = norm( Y(:), inf);
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;
Z = Y; % second multiplier
dt_dual = zeros(size(J,2),1);
M_dual = zeros(size(Dprime,1),1);

mu1 = 1/norm(Dprime) ;
mu2 = 0.01;
rho = 1.2;

d_norm = norm(Dprime);

iter = 0;
converged = false;

while ~converged       
    iter = iter + 1;
    
    % Solve subproblem 1 w.r.t S    
    temp_T = (1/(mu1+mu2))*(mu1*(Dprime - J*dt_dual - (1/mu1)*Y) + mu2*(M_dual + (1/mu2)*Z));
    % Regular mode, just soft-thresholding
    S_dual = sign(temp_T) .* pos(abs(temp_T) - 2/(mu1+mu2));
    
    % Solve subproblem 2 w.r.t. M
    temp_Timg = reshape(S_dual-(1/mu2)*Z,imgSize(1),imgSize(2));
    temp_T = wavetr(temp_Timg);
    M_dualCoef = temp_T;
    % soft-thresholding on the wavelet transform coefficients (difference coefficients only)
    subtemp_T1 = temp_T(2*imgSize(1)+1:3*imgSize(1),:);
    subtemp_T2 = temp_T(3*imgSize(1)+1:end,:);
    [M_dualCoef(2*imgSize(1)+1:3*imgSize(1),:), M_dualCoef(3*imgSize(1)+1:end,:)] = threshold(subtemp_T1, subtemp_T2, 4*sqrt(2)*lambda/mu2);
    M_dual = vec(reshape(iwavetr(M_dualCoef),imgSize(1),imgSize(2)));
    
    % Solve subproblem 3 w.r.t dt_dual
    temp_T = Dprime - S_dual - (1/mu1)*Y;
    dt_dual =  pinv(J)*temp_T;
    
    % Update Lagrangian (Y/X)
    S = -Dprime + J*dt_dual + S_dual;
    S2 = M_dual - S_dual;
    Y = Y + mu1*S;
    Z = Z + mu2*S2;
        
    mu1 = mu1*rho;
    mu2 = mu2*rho;
    stoppingCriterion = norm(S) / d_norm;
    stoppingCriterion2 = norm(S2) / d_norm;
    
    if mod( iter, DISPLAY_EVERY) == 0
        disp(['#Iteration ' num2str(iter)  ...
            ' ||S||_0 ' num2str(length(find(abs(S_dual)>0)))...
            '  Stopping Criterion ' num2str(stoppingCriterion) ...
            ' Stopping Criterion 2 ' num2str(stoppingCriterion2)]);
        disp(['1/mu=' num2str(1/mu1)]);
    end    
    
    if stoppingCriterion <= tol && stoppingCriterion <= tol,
        disp('RASL inner loop is converged at:');
        disp(['Iteration ' num2str(iter) ...
            ' ||S||_0 ' num2str(length(find(abs(S_dual)>0)))  '  Stopping Criterion ' ...
            num2str(stoppingCriterion) ' Stopping Criterion 2 ' ...
            num2str(stoppingCriterion2)]) ;
        converged = true ;
    end
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end
