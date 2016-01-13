%% Testing recursiveME
close all; clear all;
% load('Stefan_3_Frames.mat');
load('Office1_QP22.mat');
src = double(Ymtx_Cur);
ref = double(Ymtx_Ref);
inc = double(Ymtx_Inc); 
[height, width] = size(src);
totPix = height*width;

% L1-solver para
para.tauModel = 'AFFINE';
para.tau0 = eye(3);
para.numPts = 10;
para.tolInner = 1e-4;
para.maxInnerIter = 100;
para.tolOuter = 1e-5;
para.maxOuterIter = 200;
para.tolDeltaXi = 1e-5;
para.lambda = 1;
para.objMask = ones(height,width);
para.TVmode = 0; % enable TV-L1 solver
para.verbose = 1; % show debug info

% segmentation parameters, used only when maskMode==1
segpara.Debugging_Enabled = 0;
segpara.Conn_area = round(0.01*height*width); % now is the absolute pixel count
segpara.QP = 22;
segpara.Displacement = 5;
segpara.channel = 1;
segpara.ByPassFilling = 1;
segpara.QuadTreeMode = 0;

% recursive parameters
recurPara.MElayer = 0;
recurPara.maskMode = 1;
recurPara.thresh_outlier = 10; % threshold on error level for determining outlier pixels
recurPara.inlier_cnt_percent = 0.85; % percentage of inlier pixels
recurPara.max_recur = 5; 
recurPara.Y_channel = 1;
recurPara.IsDebug = 1;

%% Deriving ObjMasks
% first layer ObjMask
[predFrameSrc, errFrameSrc, ~, objMaskL1, ~] = SolveL1MErecursive(src, ref, para, segpara, recurPara);

% second layer ObjMask
recurPara.MElayer = 1;
recurPara.maskMode = 1;
recurPara.objMask = objMaskL1;
% assign smaller percentage at object level
recurPara.inlier_cnt_percent = 0.85; recurPara.max_recur = 15; % recurPara.thresh_outlier = 8;
segpara.QuadTreeMode = 1;
[predFrameSrcL2, errFrameSrcL2, ~, objMaskL2, qTreeCent] = SolveL1MErecursive(src, predFrameSrc, para, segpara, recurPara);

% Now derives the predicted frame for the incoming frames
%% Deriving first layer tau
recurPara.maskMode = 0; 
recurPara.MElayer = 0;
[predFrame, ~, tauL0, ~, ~] = SolveL1MErecursive(inc, src, para, segpara, recurPara);

% simple inpainting for the out-of-boundary region by src 
trans = floor(abs(tauL0(1:2,3))).*sign(tauL0(1:2,3)); % only the translation part in tau
if trans(1)~=0 && trans(2)==0,
    % translation along horizontal direction
    if trans(1) < 0,
        % patch left part
        predFrame(:,1:abs(trans(1))) = src(:,1:abs(trans(1)));
    else
        % patch right part
        predFrame(:,end-trans(1):end) = src(:,end-trans(1):end);
    end
elseif trans(1)==0 && trans(2)~=0,
    % trasnlation along vertical direction
    if trans(2) < 0,
        % patch top part
        predFrame(1:abs(trans(2)),:) = src(1:abs(trans(2)),:);
    else
        % patch bottom part
        predFrame(end-trans(2):end,:) = src(end-trans(2):end,:);
    end
end
errFrame = inc - predFrame;

%% Deriving second layer tau
% Now, check qTreeCent, if qTreeCent is non-empty, define quad-tree split
% mode, apply objMask in four quadrants and apply objMaskL2 within each
% quadrant separately.
for objIdx = 1:length(objMaskL1),
    qTree = qTreeCent{objIdx};
    objMaskL1sub = objMaskL1{objIdx};
    objMaskL2sub = objMaskL2{objIdx};
    if isempty(qTree),
        % non split mode, directly apply objMask and objMaskL2
        recurPara.MElayer = 1;
        recurPara.objMask = objMaskL1sub;
        [predFrameL1, errFrameL1, tauL1, ~, ~] = SolveL1MErecursive(inc, predFrame, para, segpara, recurPara);
        
        recurPara.MElayer = 2;
        recurPara.objMask = objMaskL2sub;
        [predFrameL2, errFrameL2, tauL2, ~, ~] = SolveL1MErecursive(inc, predFrameL1, para, segpara, recurPara);
    else
        % split mode
        Xcentroid = qTree(1); Ycentroid = qTree(2);
        Xstart = [1 1 Xcentroid+1 Xcentroid+1];
        Xend = [Xcentroid Xcentroid height height];
        Ystart = [1 Ycentroid+1 1 Ycentroid+1];
        Yend = [Ycentroid width Ycentroid width];
        
        % per quadrant mask
        curr_objMaskQuad = cell(4,1);
        [curr_objMaskQuad{:}] = deal(zeros(height,width));
        curr_objMaskQuad{1}(1:Xcentroid,1:Ycentroid) = objMaskL1sub(1:Xcentroid,1:Ycentroid);
        curr_objMaskQuad{2}(1:Xcentroid,Ycentroid+1:end) = objMaskL1sub(1:Xcentroid,Ycentroid+1:end);
        curr_objMaskQuad{3}(Xcentroid+1:end,1:Ycentroid) = objMaskL1sub(Xcentroid+1:end,1:Ycentroid);
        curr_objMaskQuad{4}(Xcentroid+1:end,Ycentroid+1:end) = objMaskL1sub(Xcentroid+1:end,Ycentroid+1:end);
        
        % save tauL1 per quadrant as well
        tauL1 = cell(4,1);
        predFrameL1 = predFrame;
        errFrameL1 = zeros(height,width);
        
        % Now apply each one successively
        recurPara.MElayer = 1;
        recurPara.objMask = curr_objMaskQuad;
        [predFrameL1, errFrameL1, tauL1, ~, ~] = SolveL1MErecursive(inc, predFrame, para, segpara, recurPara);
        
        % Once L1 is done, go for L2, and apply for each quadrant as well
        tauL2 = cell(4,1);
        predFrameL2 = predFrameL1;
        
        for qIdx = 1:4,
            recurPara.MElayer = 2;
            if ~isempty(objMaskL2sub{qIdx}),
                recurPara.objMask = objMaskL2sub{qIdx};
                
                [predFrameL2q, errFrameL2q, tauL2{qIdx}, ~, ~] = SolveL1MErecursive(inc, predFrameL1, para, segpara, recurPara);
                predFrameL2(Xstart(qIdx):Xend(qIdx),Ystart(qIdx):Yend(qIdx)) = predFrameL2q(Xstart(qIdx):Xend(qIdx),Ystart(qIdx):Yend(qIdx));
            end
        end
        
        errFrameL2 = inc - predFrameL2;
        
    end
end