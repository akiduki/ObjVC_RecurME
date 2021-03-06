function [predFrame, errFrame, tau, objMask, qTreeCent] = SolveL1MErecursive(srcFrame, refFrame, refFrame_org, L1solverPara, segPara, recurPara)
% Author : Yuanyi Xue @ NYU-Poly
% The upper warpper function for solving the L1 ME recursively.
% The function takes in the referece frame and source frame, with
% parameters for L1 and recursive criteron, on three structs, respectively,
% 
% Inputs:
% 1) srcFrame/refFrame - the source frame, reference frame in matrices
% 1b) refFrame_org - original reference frame, only used when object ME is inverse estimated.
% 2) L1solverPara - Parameter struct for L1 solver, check L1 solver for
%                   parameter definination
% 3) segPara      - Parameter struct for object segmentation, used when
%                   maskMode==1 only
% 4) recurPara    - Parameter controls the recursive ME logic
%    .MElayer            - decite which layer it is, =0 means gloabl layer
%    .maskMode           - =1 to enable deriving object masks
%    .thresh_outlier     - threshold for determining outlier pixels
%    .inlier_cnt_percent - threshold on percentage of inlier pixels
%    .max_recur          - maximum number of recursives allowed
%    .Y_channel          - ==1 means Y_channel, for UV use 0
%    .objMask            - object binary mask, same size as input frame
%    .transform          - transform is only used for U/V channels
% Outputs:
% 1) predFrame    - the motion compensated frame
% 2) errFrame     - the error frame of the predicted frame
% 3) tau          - motion parameter
% 4) objMask      - object mask derived from segmentation procedure, a
%                   cell array contains the objMask per element
% 5) qTreeCentroid  - centroid coordinates if quad-tree flag is on.
% *Note the predFrame, errFrame are with respect to the input size.
% *transform is derived with respect to the center of input images.
% 4-5 are only used when maskMode==1 

% recursive L1 ME parameters
MElayer = recurPara.MElayer;
maskMode = recurPara.maskMode; % maskMode=1 for deriving object masks only
thresh_outlier = recurPara.thresh_outlier; % threshold on error level for determining outlier pixels
inlier_cnt_percent = recurPara.inlier_cnt_percent; % percentage of inlier pixels
max_recur = recurPara.max_recur; % maximal recursive time
Y_channel = recurPara.Y_channel;
IsDebug = recurPara.IsDebug; % debug flag, used to visualize the intermediate result

EnableInvTrans = para.EnableInvTrans;

qTreeCent = []; % quad-tree flag

% Object segmentation parameters
if isfield(segPara,'dilate_width'),
    min_dilate_width = segPara.dilate_width;
else
    min_dilate_width = 2;
end
QuadTreeMode = segPara.QuadTreeMode; % QuadTreeMode enabling flag

% Image size
[height, width] = size(srcFrame);
totPix = height*width;
% Make sure frames are in double
srcFrame = double(srcFrame);
refFrame = double(refFrame);

if Y_channel == 1,
    if MElayer == 0,
        % Global ME level, the mask is set to all 1 by default.
        L1solverPara.objMask = recurPara.GlobalMask; % GlobalMask is the compliment of the first layer object mask
        
        % Calling the lern2frmtau for the first time to derive outliers
        [transform, err, pred] = lern2frmtau(srcFrame, refFrame, L1solverPara);
        
        % for deriving bounding box, shall do ME recursively
        if maskMode,
            % Global ME level, the mask is set to all 1 by default.
            L1solverPara.objMask = ones(height, width);
            
            % Define the transform origin
            imgSize = [height, width];
            
            [transform, err, pred] = computeRecursiveME(srcFrame, refFrame, err, L1solverPara, thresh_outlier, inlier_cnt_percent, max_recur, imgSize, QuadTreeMode);
            % bounding box mode, needs derive bounding box from err
            warning('DERIVING BOUNDING BOX MODE ENABLED!');
            warning('CHECK WHETHER THE INPUTS ARE RECONSTRUCTED FRAME!');
            Mask = postprocessing(srcFrame, err, segPara);
            dilate_width = round(max([min_dilate_width max(abs(transform(1:2,3)'))]));
            se = strel('disk', dilate_width);
            for i=1:length(Mask);
                objMask{i} = imdilate(Mask{i}, se);
            end
        else
            objMask = [];
        end
        
        % outputs the result
        predFrame = pred;
        errFrame = err;
        tau = transform;
    else
        predFrame = refFrame;
        % sub-frame object layers
        % crop the sub-frames from boundingbox
        for objID=1:length(recurPara.objMask),
            curr_objMask = recurPara.objMask{objID};
            imgSize = size(srcFrame);
            xdata = [1 imgSize(2)] - ceil(imgSize(2)/2);
            ydata = [1 imgSize(1)] - ceil(imgSize(1)/2);
            
            L1solverPara.objMask = curr_objMask;
            
            % for deriving bounding box, use recursive ME
            if maskMode,
                if IsDebug,
                    errOld = err; % store the old error for comparison only
                end
                
                % Calling the lern2frmtau, old version for maskMode
                [transform, err, pred] = lern2frmtau(srcFrame, refFrame, L1solverPara);
                
                % Put each piece back to the reconstructed whole frame
                predFrame(curr_objMask(:)==1) = pred(curr_objMask(:)==1);
                
                predFrame = refFrame;
                [transform, err, pred, inlier_objMask, quadtreeFlag] = ...
                    computeRecursiveME(srcFrame, refFrame, err, L1solverPara, thresh_outlier, inlier_cnt_percent, max_recur, imgSize, QuadTreeMode);
                
                % bounding box mode, needs derive bounding box from err
                warning('DERIVING BOUNDING BOX MODE ENABLED!');
                warning('CHECK WHETHER THE INPUTS ARE RECONSTRUCTED FRAME!');
                
                if ~isempty(quadtreeFlag),
                    % use quad-tree decomposition
                    Xcentroid = quadtreeFlag(1);
                    Ycentroid = quadtreeFlag(2);
%                     srcFrameObjMask{1} = srcFrame(1:Xcentroid,1:Ycentroid);
%                     srcFrameObjMask{2} = srcFrame(1:Xcentroid,Ycentroid+1:end);
%                     srcFrameObjMask{3} = srcFrame(Xcentroid+1:end,1:Ycentroid);
%                     srcFrameObjMask{4} = srcFrame(Xcentroid+1:end,Ycentroid+1:end);
                    % decompose the original mask as well
                    curr_objMaskQuad = cell(4,1);
                    errMask = cell(4,1);
                    [curr_objMaskQuad{:}] = deal(zeros(imgSize(1),imgSize(2)));
                    curr_objMaskQuad{1}(1:Xcentroid,1:Ycentroid) = curr_objMask(1:Xcentroid,1:Ycentroid);
                    curr_objMaskQuad{2}(1:Xcentroid,Ycentroid+1:end) = curr_objMask(1:Xcentroid,Ycentroid+1:end);
                    curr_objMaskQuad{3}(Xcentroid+1:end,1:Ycentroid) = curr_objMask(Xcentroid+1:end,1:Ycentroid);
                    curr_objMaskQuad{4}(Xcentroid+1:end,Ycentroid+1:end) = curr_objMask(Xcentroid+1:end,Ycentroid+1:end);
                    
                    for i=1:4, %<---- hard coded 4 here
                        errMask{i} = err{i}.*curr_objMaskQuad{i};
                        currFrameObjMask = srcFrame.*curr_objMaskQuad{i};
                        
                        MaskQuad = postprocessing(currFrameObjMask, errMask{i}, segPara);
%                         for ii=1:length(boundingboxQuad),
%                             boundingboxQuad{ii}(1) = boundingboxQuad{ii}(1) + curr_boundingbox(1);
%                             boundingboxQuad{ii}(3) = boundingboxQuad{ii}(3) + curr_boundingbox(3);
%                         end
                        se = strel('disk', min_dilate_width);
                        for ii=1:length(MaskQuad);
                            MaskQuad{ii} = imdilate(MaskQuad{ii}, se);
                        end
                        
                        MaskSub{i} = MaskQuad;
                        
                        % update predFrame by current quadrants
                        predFrame(curr_objMaskQuad{i}(:)==1) = pred(curr_objMaskQuad{i}(:)==1);
                    end
                    
                    if IsDebug,
                        % Assemble the masked error term from Quad-tree
                        errMaskAgg = errMask{1}+errMask{2}+errMask{3}+errMask{4};
                        figure;
                        subplot(1,2,1);imshow(errOld,[-255,255]);
                        title(['Error before Quad-tree, L1=' num2str(norm(errOld(:),1)) ' L2=' num2str(norm(errOld(:),2))]);
                        subplot(1,2,2);imshow(errMaskAgg,[-255 255]);
                        title(['Error after Quad-tree, L1=' num2str(norm(errMaskAgg(:),1)) ' L2=' num2str(norm(errMaskAgg(:),2))]);
                    end
                    
                else
                    errMask = zeros(imgSize(1),imgSize(2));
                    % further layers shall do on masked err only
                    errMask = err.*curr_objMask;
                    srcFrameObjMask = srcFrame.*curr_objMask;
                    MaskSub = postprocessing(srcFrameObjMask, errMask, segPara);
                    
                    se = strel('disk', min_dilate_width);
                    for i=1:length(MaskSub);
                        MaskSub{i} = imdilate(MaskSub{i}, se);
                    end
                    
                    % update predFrame by new inlier mask
                    predFrame(inlier_objMask(:)==1) = pred(inlier_objMask(:)==1);
                    
                end
                
                % save it according to the original organizations
                tau{objID} = transform;
                objMask{objID} = MaskSub;
                qTreeCent{objID} = quadtreeFlag;
            else
                % Calling the lern2frmtau for coding mode
                if EnableInvTrans,
                    % Inverse mode, compute the tau inversely
                    [rr_tran, ~, pred] = lern2frmtau(refFrame, srcFrame, L1solverPara);
                    
                    % this is the forward transform
                    fw_trm = maketform('projective',rr_tran');
                    curr_objMaskT = imtransform(curr_objMask, fw_trm, 'nearest', 'XData', xdata, 'YData', ydata, 'UData', xdata, 'VData', ydata, 'Size', imgSize, 'fill', 0);
                    
                    % warp the reference frame to the source by using fw_trm as well
                    predFrame_objMaskT = imtransform(refFrame, fw_trm, 'bicubic', 'XData', xdata, 'YData', ydata, 'UData', xdata, 'VData', ydata, 'Size', imgSize, 'fill', 0);
                    % now update the predFrame only within the masked region
                    predFrame(curr_objMaskT(:)==1) = predFrame_objMaskT(curr_objMaskT(:)==1);
                    
                    % store the tau
                    tau{objID} = rr_tran;
                else
                    % Regular mode
                    [tran, ~, pred] = lern2frmtau(srcFrame, refFrame, L1solverPara);
                    
                    % Update the predFrame only within the masked region
                    predFrame(curr_objMask(:)==1) = pred(curr_objMask(:)==1);
                    tau{objID} = tran;
                end
                objMask = [];
            end
        end
        % output error frame
        errFrame = srcFrame - predFrame;
    end
else
    % for U/V frames, apply transform only
    imgSize = [height, width];
    xdata = [1 imgSize(2)] - ceil(imgSize(2)/2);
    ydata = [1 imgSize(1)] - ceil(imgSize(1)/2);
    
    warning('U/V channel mode, transform and bounding box shall be scaled already!');
    if MElayer==0,
        transform = recurPara.transform;
        % apply transform to all image pixels
        Tfm = fliptform(maketform('projective',transform'));
        predFrame = imtransform(refFrame, Tfm,'bicubic','XData', xdata, 'YData', ydata, ...
            'UData', xdata, 'VData', ydata, 'Size', imgSize);
        errFrame = srcFrame - predFrame;
        tau = [];
        objMask = [];
    else
        predFrame = refFrame;
        % sub-frame object layers
        % crop the sub-frames from boundingbox
        for objID=1:length(recurPara.objMask),
            curr_objMask = recurPara.objMask{objID};
            curr_transform = recurPara.transform{objID};  
            
            Tfm = fliptform(maketform('projective',curr_transform'));
            pred = imtransform(refFrame, Tfm,'bicubic','XData', xdata, 'YData', ydata, ...
                'UData', xdata, 'VData', ydata, 'Size', imgSize);
            
            predFrame(curr_objMask==1) = pred(curr_objMask==1);
        end
        errFrame = srcFrame - predFrame;
        tau = [];
        objMask = [];
    end
end
end

% Nested function for recursive ME
function [transform, err, predFrame, inlier_objMask, quadtreePart] = computeRecursiveME(srcFrame, refFrame, err, L1solverPara, thresh_outlier, inlier_cnt_percent, max_recur, imgSize, QuadTreeMode)
% TO-DO: update inlier_objMask within the area object mask defined.
predFrame = refFrame;
originalObjMask = L1solverPara.objMask; % sub-frame object region
obj_idx = find(originalObjMask(:)==1);
obj_pix_cnt = length(obj_idx);
xdata = [1 imgSize(2)] - ceil(imgSize(2)/2);
ydata = [1 imgSize(1)] - ceil(imgSize(1)/2);

% Deriving the outlier pixels and update the mask
inlier_objMask = zeros(imgSize(1),imgSize(2));
inlier_idx = find(abs(err(:))<thresh_outlier);
inlier_obj_idx = intersect(obj_idx,inlier_idx,'legacy');
inlier_objMask(inlier_obj_idx) = 1;
IsConverged = false;
recur_cnt = 0;
prev_inlier_ratio = 0;
quadtreePart = [];
while ~IsConverged,
    % for each recursive iteration, do the following:
    % 1) use the inlier_objMask to find the motion through lern2frmtau
    % 2) apply learned ME for all pixels
    % 3) evaluate the new inlier pixels, update inlier_objMask
    % 4) if not converged and recur<max_recur, go 1)
    L1solverPara.objMask = inlier_objMask;
    [transform, ~, ~, ~] = lern2frmtau(srcFrame, refFrame, L1solverPara);
    L1solverPara.tau0 = transform;
    
    % apply transform to all image pixels
    Tfm = fliptform(maketform('projective',transform'));
    pred = imtransform(refFrame, Tfm,'bicubic','XData', xdata, 'YData', ydata, ...
        'UData', xdata, 'VData', ydata, 'Size', imgSize);
    % new error image
    predFrame(obj_idx) = pred(obj_idx);
    err = srcFrame - predFrame;
    
    % Convergence check
    inlier_idx = find(abs(err(:))<thresh_outlier);
    inlier_obj_idx = intersect(obj_idx,inlier_idx,'legacy');
    inlier_cnt = length(inlier_obj_idx);
    inlier_ratio = inlier_cnt/obj_pix_cnt;
    if inlier_ratio >= inlier_cnt_percent,
        IsConverged = true;
    else
        % not converged
        inlier_objMask = zeros(imgSize(1),imgSize(2));
        inlier_objMask(inlier_obj_idx) = 1;
        if (inlier_ratio-prev_inlier_ratio)*obj_pix_cnt <= 10 && QuadTreeMode && obj_pix_cnt >= 1600,
            % check whether there is a significant amount of pixel changes
            % if not, partition into quad-tree
            
            % get the centroid coordinates
            XidxMat = repmat([1:imgSize(1)]',1,imgSize(2));
            YidxMat = repmat([1:imgSize(2)],imgSize(1),1);
            Xcentroid = round(sum(sum(originalObjMask.*XidxMat))/obj_pix_cnt);
            Ycentroid = round(sum(sum(originalObjMask.*YidxMat))/obj_pix_cnt);
            quadtreePart = [Xcentroid Ycentroid];
            err = cell(4,1);
            transform = cell(4,1);
            inlier_objMask_quad = cell(4,1);
            originalObjMask_quad = cell(4,1);
            [inlier_objMask_quad{:}] = deal(zeros(imgSize(1),imgSize(2)));
            [originalObjMask_quad{:}] = deal(zeros(imgSize(1),imgSize(2)));
            
            predFrameQuad = refFrame;
            
            % Now spilt into four quadrants
            disp('Quad-tree split occured!')
            inlier_objMask_quad{1}(1:Xcentroid,1:Ycentroid) = inlier_objMask(1:Xcentroid,1:Ycentroid);
            inlier_objMask_quad{2}(1:Xcentroid,Ycentroid+1:end) = inlier_objMask(1:Xcentroid,Ycentroid+1:end);
            inlier_objMask_quad{3}(Xcentroid+1:end,1:Ycentroid) = inlier_objMask(Xcentroid+1:end,1:Ycentroid);
            inlier_objMask_quad{4}(Xcentroid+1:end,Ycentroid+1:end) = inlier_objMask(Xcentroid+1:end,Ycentroid+1:end);
            originalObjMask_quad{1}(1:Xcentroid,1:Ycentroid) = originalObjMask(1:Xcentroid,1:Ycentroid);
            originalObjMask_quad{2}(1:Xcentroid,Ycentroid+1:end) = originalObjMask(1:Xcentroid,Ycentroid+1:end);
            originalObjMask_quad{3}(Xcentroid+1:end,1:Ycentroid) = originalObjMask(Xcentroid+1:end,1:Ycentroid);
            originalObjMask_quad{4}(Xcentroid+1:end,Ycentroid+1:end) = originalObjMask(Xcentroid+1:end,Ycentroid+1:end);
            
            % do it for each quadrant
            for i=1:4,
                % per-quadrant constants
                curr_inlier_objMask = inlier_objMask_quad{i};
                imgSize = size(srcFrame);
                obj_idx = find(originalObjMask_quad{i}==1); % mask w.r.t original pixels
                obj_pix_cnt = length(obj_idx);
                xdata = [1 imgSize(2)] - ceil(imgSize(2)/2);
                ydata = [1 imgSize(1)] - ceil(imgSize(1)/2);
                L1solverParaQuad = L1solverPara;
                prev_inliner_ratio = 0;
                IsConverged = false;
                recur_cnt_quad = 0;
                
                
                while ~IsConverged,
                    L1solverParaQuad.objMask = curr_inlier_objMask;
                    [tQuad, ~, ~, ~] = lern2frmtau(srcFrame, refFrame, L1solverParaQuad);
                    L1solverParaQuad.tau0 = tQuad;
                    
                    % apply transform to all image pixels
                    Tfm = fliptform(maketform('projective',tQuad'));
                    predQuad = imtransform(refFrame, Tfm,'bicubic','XData', xdata, 'YData', ydata, ...
                        'UData', xdata, 'VData', ydata, 'Size', imgSize);
                    % new error image
                    predFrameQuad(obj_idx) = predQuad(obj_idx);
                    errQuad = srcFrame - predFrameQuad;
                    
                    % Convergence check
                    inlier_idx = find(abs(errQuad(:))<thresh_outlier);
                    inlier_obj_idx = intersect(obj_idx,inlier_idx,'legacy');
                    inlier_cnt = length(inlier_obj_idx);
                    inlier_ratio = inlier_cnt/obj_pix_cnt;
                    
                    if inlier_ratio >= inlier_cnt_percent || (prev_inlier_ratio-inlier_ratio) < 0.005, % <======= less than 0.5% change
                        IsConverged = true;
                    else
                        % not converged, update the inlier objMask
                        curr_inlier_objMask = zeros(imgSize(1),imgSize(2));
                        curr_inlier_objMask(inlier_obj_idx) = 1;
                        prev_inlier_ratio = inlier_ratio;
                    end
                    % Increment recursive counter
                    recur_cnt_quad = recur_cnt_quad + 1;
                    % if maximum recursive count is reached
                    if recur_cnt_quad > max_recur,
                        disp('Maximum recursive number reached!');
                        break;
                    end
                end
                
                % Save per-quadrant data back
                err{i} = errQuad.*originalObjMask_quad{i};
                transform{i} = tQuad;
                inlier_objMask_quad{i} = curr_inlier_objMask;
            end
            inlier_objMask = cell(4,1);
            inlier_objMask = inlier_objMask_quad;
            
            predFrame = predFrameQuad;
            break; % after processed the four quadrants, quit the loop
        end
        
        prev_inlier_ratio = inlier_ratio;
    end
    
    % Increment recursive counter
    recur_cnt = recur_cnt + 1;
    % if maximum recursive count is reached
    if recur_cnt > max_recur,
        disp('Maximum recursive number reached!');
        break;
    end
end
end