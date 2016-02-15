%% RecursiveME for Object Video Coding Project
% Yuanyi Xue @ NYU-Poly

% TO-DO:
% 1) revamp the global layer ME so that the outlier is excluded
% 2) add an option to do non-hierachical parametric ME
% 3) apply the inverse ME when driving the outliers as well (object masks)

% close all; 
clear all;

% seqName = {'Stefan', 'Office1', 'City', 'Football'};
seqName = {'Stefan'};
% QP = [22 27 32 37];
QP = 22;

for seqID = 1:length(seqName),
    for QPidx = 1:length(QP),
        currData = ['./Data/' seqName{seqID} '_QP' int2str(QP(QPidx)) '.mat'];
        load(currData);
        src = double(Ymtx_Cur);ref = double(Ymtx_Ref);inc = double(Ymtx_Inc);
        srcU = double(Umtx_Cur);incU = double(Umtx_Inc);
        srcV = double(Vmtx_Cur);incV = double(Vmtx_Inc);
        [height, width] = size(src);
        totPix = height*width;
        
        % Global control parameters
        Enable3rdObjLayer = 0;
        EnableGlobalMask  = 1; % global mask, used to exclude the first object layer when predicting global motion
        
        % L1-solver para
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
        segpara.QP = QP(QPidx);
        segpara.Displacement = 5;
        segpara.channel = 1;
        segpara.ByPassFilling = 1;
        
        % recursive parameters
        recurPara.thresh_outlier = 10; % threshold on error level for determining outlier pixels
        recurPara.inlier_cnt_percent = 0.85; % percentage of inlier pixels
        recurPara.max_recur = 5;
        recurPara.Y_channel = 1;
        recurPara.IsDebug = 0;
        recurPara.GlobalMask = ones(height, width); % default global ME mask
        
        %% Deriving ObjMasks
        % first layer ObjMask
        recurPara.MElayer = 0;
        recurPara.maskMode = 1;
        segpara.QuadTreeMode = 0; % No need to do quadtree for the first layer
        segpara.ByPassAddMorpFilt = 0;
        % segpara.Conn_area = round(0.005*height*width); % or 0.01/0.02, here is the absolute number of pixels now
        segpara.Conn_area = 500;
        segpara.ConfThrLo =  10; % 8 for city 6;
        % para.tauModel = 'AFFINE';
        para.tauModel = 'HOMOGRAPHY';
        [predFrameSrc, errFrameSrc, ~, objMaskL1, ~] = SolveL1MErecursive(src, ref, ref, para, segpara, recurPara);
        
        % second layer ObjMask
        recurPara.MElayer = 1;
        recurPara.maskMode = 1;
        recurPara.objMask = objMaskL1;
        % assign smaller percentage at object level
        recurPara.inlier_cnt_percent = 0.80; recurPara.max_recur = 15; % recurPara.thresh_outlier = 8;
        segpara.QuadTreeMode = 1;
        segpara.ByPassAddMorpFilt = 0; % at finer object level, bypass additional morphological filtering
        segpara.Conn_area = 150; % for Quad-tree, needs siginificantly smaller Conn_area
        % segpara.ConfThrLo = 4;
        para.tauModel = 'AFFINE';
        [predFrameSrcL2, errFrameSrcL2, ~, objMaskL2, qTreeCent] = SolveL1MErecursive(src, predFrameSrc, ref, para, segpara, recurPara);
        
        % third layer ObjMask
        if Enable3rdObjLayer,
            recurPara.MElayer = 2;
            recurPara.maskMode = 1;
            % assign smaller percentage at object level
            recurPara.inlier_cnt_percent = 0.85; recurPara.max_recur = 15;
            segpara.Conn_area = 100;
            para.tauModel = 'HOMOGRAPHY';
            % initialize the objMaskL3 and qTreeCentL2
            objMaskL3 = cell(length(objMaskL2),1);
            qTreeCentL2 = objMaskL3;
            
            for objIdx=1:length(objMaskL2),
                objMaskL2sub = objMaskL2{objIdx};
                for j=1:length(objMaskL2sub),
                    if ~isempty(objMaskL2sub{j}),
                        recurPara.objMask = objMaskL2sub{j};
                        [predFrameSrcL3, errFrameSrcL3, ~, objMaskL3sub, qTreeCentL2sub] = SolveL1MErecursive(src, predFrameSrcL2, ref, para, segpara, recurPara);
                        objMaskL3{objIdx}{j} = objMaskL3sub;
                        qTreeCentL2{objIdx}{j} = qTreeCentL2sub;
                    else
                        objMaskL3{objIdx}{j} = [];
                        qTreeCentL2{objIdx}{j} = [];
                    end
                end
            end
        end
        
        % Now derives the predicted frame for the incoming frames
        %% Deriving first layer tau
        recurPara.maskMode = 0;
        recurPara.MElayer = 0;
        para.tauModel = 'HOMOGRAPHY';
        if EnableGlobalMask,
            tmpMask = zeros(height, width);
            for objIdx = 1:length(objMaskL1),
                tmpMask = tmpMask + objMaskL1{objIdx};
            end
            recurPara.GlobalMask = (tmpMask==0);
        end
        
        [pred, ~, tauL0, ~, ~] = SolveL1MErecursive(inc, src, src, para, segpara, recurPara);
        
        % get the frame back in case the global mask is enabled
        if EnableGlobalMask,
            predFrame = src;
            predFrame(recurPara.GlobalMask==1) = pred(recurPara.GlobalMask==1);
        else
            predFrame = pred;
        end
        
        % simple inpainting for the out-of-boundary region by src
        trans = ceil(abs(tauL0(1:2,3))).*sign(tauL0(1:2,3)); % only the translation part in tau
        if trans(1)~=0,
            % translation along horizontal direction
            if trans(1) < 0,
                % patch left part
                predFrame(:,1:abs(trans(1))) = src(:,1:abs(trans(1)));
            else
                % patch right part
                predFrame(:,end-trans(1):end) = src(:,end-trans(1):end);
            end
        end
        
        if trans(2)~=0,
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
        tauL1 = cell(length(objMaskL1),1); tauL2 = tauL1;
        para.tauModel = 'AFFINE';
        for objIdx = 1:length(objMaskL1),
            qTree = qTreeCent{objIdx};
            objMaskL1sub = objMaskL1{objIdx};
            
            % Copy whatever has already done
            if objIdx==1,
                predL1 = predFrame;
            else
                predL1 = predFrameL1;
            end
            
            % Check whether it is a quad-tree split
            if isempty(qTree),
                % non split mode, directly apply objMask
                recurPara.MElayer = 1;
                recurPara.objMask = {objMaskL1sub};
                [predFrameL1, ~, tauL1{objIdx}, ~, ~] = SolveL1MErecursive(inc, predL1, src, para, segpara, recurPara);
                
            else
                % quad-tree split mode
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
                
                % Now apply each one successively
                recurPara.MElayer = 1;
                recurPara.objMask = curr_objMaskQuad;
                [predFrameL1, ~, currL1, ~, ~] = SolveL1MErecursive(inc, predL1, src, para, segpara, recurPara);
                tauL1{objIdx} = currL1;
                
            end
        end
        errFrameL1 = inc - predFrameL1;
        
        for objIdx = 1:length(objMaskL1),
            % Go on for the second object layer
            qTree = qTreeCent{objIdx};
            objMaskL2sub = objMaskL2{objIdx};
            
            % For first object in L2, initiate the frames
            if objIdx==1,
                predL2 = predFrameL1;
                predFrameL2 = predFrameL1;
            else
                predL2 = predFrameL2;
            end
            
            if isempty(qTree),
                % non split mode, directly apply objMaskL2
                if ~isempty(objMaskL2sub),
                    recurPara.MElayer = 2;
                    recurPara.objMask = objMaskL2sub;
                    [predFrameL2, ~, tauL2{objIdx}, ~, ~] = SolveL1MErecursive(inc, predL2, src, para, segpara, recurPara);
                end
                
            else
                % split mode, go for each quadrant as well
                currL2 = cell(1,4);
                % quad-tree split mode
                Xcentroid = qTree(1); Ycentroid = qTree(2);
                Xstart = [1 1 Xcentroid+1 Xcentroid+1];
                Xend = [Xcentroid Xcentroid height height];
                Ystart = [1 Ycentroid+1 1 Ycentroid+1];
                Yend = [Ycentroid width Ycentroid width];
                
                for qIdx = 1:4,
                    recurPara.MElayer = 2;
                    if ~isempty(objMaskL2sub{qIdx}),
                        recurPara.objMask = objMaskL2sub{qIdx};
                        
                        [predFrameL2q, ~, currL2{qIdx}, ~, ~] = SolveL1MErecursive(inc, predL2, src, para, segpara, recurPara);
                        predFrameL2(Xstart(qIdx):Xend(qIdx),Ystart(qIdx):Yend(qIdx)) = predFrameL2q(Xstart(qIdx):Xend(qIdx),Ystart(qIdx):Yend(qIdx));
                    end
                end
                tauL2{objIdx} = currL2;
                
            end
        end
        errFrameL2 = inc - predFrameL2;
        
        %% Processing U/V channels
        recurPara.Y_channel = 0; % 0=UV channels
        
        % global layer
        recurPara.MElayer = 0;
        recurPara.transform = tauL0;
        % scale down only the translation part
        recurPara.transform(1,3) = recurPara.transform(1,3)/2; recurPara.transform(2,3) = recurPara.transform(2,3)/2;
        
        predFrameU = SolveL1MErecursive(incU, srcU, srcU, para, segpara, recurPara);
        predFrameV = SolveL1MErecursive(incV, srcV, srcV, para, segpara, recurPara);
        errFrameU = incU - predFrameU;
        errFrameV = incV - predFrameV;
        
        % first object layer
        recurPara.MElayer = 1;
        for objIdx = 1:length(objMaskL1),
            qTree = qTreeCent{objIdx};
            objMaskL1sub = objMaskL1{objIdx};
            tauL1sub = tauL1{objIdx};
            
            % Copy whatever has already done
            if objIdx==1,
                predL1U = predFrameU;
                predL1V = predFrameV;
            else
                predL1U = predFrameL1U;
                predL1V = predFrameL1V;
            end
            
            % Check whether it is a quad-tree split
            if isempty(qTree),
                % non split mode, directly apply objMask
                recurPara.MElayer = 1;
                recurPara.objMask = {imresize(objMaskL1sub, 0.5, 'nearest')};
                recurPara.transform = tauL1sub;
                % scale down only the translation part
                recurPara.transform{1}(1,3) = recurPara.transform{1}(1,3)/2; recurPara.transform{1}(2,3) = recurPara.transform{1}(2,3)/2;
                
                predFrameL1U = SolveL1MErecursive(incU, predL1U, srcU, para, segpara, recurPara);
                predFrameL1V = SolveL1MErecursive(incV, predL1V, srcV, para, segpara, recurPara);
                
            else
                % quad-tree split mode
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
                % resize for U/V
                curr_objMaskQuad = cellfun(@(x) imresize(x,0.5,'nearest'), curr_objMaskQuad, 'UniformOutput', false);
                
                % Now apply each one successively
                recurPara.MElayer = 1;
                recurPara.objMask = curr_objMaskQuad;
                recurPara.transform = tauL1sub;
                % needs to scale down the translational components in each element
                recurPara.transform{1}(1,3) = recurPara.transform{1}(1,3)/2; recurPara.transform{1}(2,3) = recurPara.transform{1}(2,3)/2;
                recurPara.transform{2}(1,3) = recurPara.transform{2}(1,3)/2; recurPara.transform{2}(2,3) = recurPara.transform{2}(2,3)/2;
                recurPara.transform{3}(1,3) = recurPara.transform{3}(1,3)/2; recurPara.transform{3}(2,3) = recurPara.transform{3}(2,3)/2;
                recurPara.transform{4}(1,3) = recurPara.transform{4}(1,3)/2; recurPara.transform{4}(2,3) = recurPara.transform{4}(2,3)/2;
                
                predFrameL1U = SolveL1MErecursive(incU, predL1U, srcU, para, segpara, recurPara);
                predFrameL1V = SolveL1MErecursive(incV, predL1V, srcV, para, segpara, recurPara);
                
            end
        end
        errFrameL1U = incU - predFrameL1U;
        errFrameL1V = incV - predFrameL1V;
        
        % second object layer
        for objIdx = 1:length(objMaskL1),
            % Go on for the second object layer
            qTree = qTreeCent{objIdx};
            objMaskL2sub = objMaskL2{objIdx};
            tauL2sub = tauL2{objIdx};
            
            recurPara.MElayer = 2;
            
            % For first object in L2, initiate the frames
            if objIdx==1,
                predL2U = predFrameL1U;
                predL2V = predFrameL1V;
                predFrameL2U = predFrameL1U;
                predFrameL2V = predFrameL1V;
            else
                predL2U = predFrameL2U;
                predL2V = predFrameL2V;
            end
            
            if isempty(qTree),
                % non split mode, directly apply objMaskL2
                if ~isempty(objMaskL2sub),
                    recurPara.objMask = cellfun(@(x) imresize(x,0.5,'nearest'), objMaskL2sub, 'UniformOutput', false);
                    for i=1:length(objMaskL2sub),
                        recurPara.transform{i} = tauL2sub{i};
                        recurPara.transform{i}(1,3) = recurPara.transform{i}(1,3)/2;
                        recurPara.transform{i}(2,3) = recurPara.transform{i}(2,3)/2;
                    end
                    
                    predFrameL2U = SolveL1MErecursive(incU, predL2U, srcU, para, segpara, recurPara);
                    predFrameL2V = SolveL1MErecursive(incV, predL2V, srcV, para, segpara, recurPara);
                end
                
            else
                % quad-tree split mode
                Xcentroid = qTree(1); Ycentroid = qTree(2);
                Xstart = round([1 1 Xcentroid+1 Xcentroid+1]/2);
                Xend = round([Xcentroid Xcentroid height height]/2);
                Ystart = round([1 Ycentroid+1 1 Ycentroid+1]/2);
                Yend = round([Ycentroid width Ycentroid width]/2);
                
                
                % split mode, go for each quadrant as well
                for qIdx = 1:4,
                    if ~isempty(objMaskL2sub{qIdx}),
                        recurPara.objMask = cellfun(@(x) imresize(x,0.5,'nearest'), objMaskL2sub{qIdx}, 'UniformOutput', false);
                        for i=1:length(objMaskL2sub{qIdx}),
                            recurPara.transform{i} = tauL2sub{qIdx}{i};
                            recurPara.transform{i}(1,3) = recurPara.transform{i}(1,3)/2;
                            recurPara.transform{i}(2,3) = recurPara.transform{i}(2,3)/2;
                        end
                        
                        predFrameL2Uq = SolveL1MErecursive(incU, predL2U, srcU, para, segpara, recurPara);
                        predFrameL2Vq = SolveL1MErecursive(incV, predL2V, srcV, para, segpara, recurPara);
                        predFrameL2U(Xstart(qIdx):Xend(qIdx),Ystart(qIdx):Yend(qIdx)) = predFrameL2Uq(Xstart(qIdx):Xend(qIdx),Ystart(qIdx):Yend(qIdx));
                        predFrameL2V(Xstart(qIdx):Xend(qIdx),Ystart(qIdx):Yend(qIdx)) = predFrameL2Vq(Xstart(qIdx):Xend(qIdx),Ystart(qIdx):Yend(qIdx));
                    end
                end
                
            end
        end
        errFrameL2U = incU - predFrameL2U;
        errFrameL2V = incV - predFrameL2V;
        
        % re-cast the error frame back to uint8 for HEVC coding
        errFrame_uint8 = uint8( dzquantRecon(errFrame,2)/2 + 128 );
        errFrameL1_uint8 = uint8( dzquantRecon(errFrameL1,2)/2 + 128 );
        errFrameL2_uint8 = uint8( dzquantRecon(errFrameL2,2)/2 + 128 );
        errFrameU_uint8 = uint8( dzquantRecon(errFrameU,2)/2 + 128 );
        errFrameL1U_uint8 = uint8( dzquantRecon(errFrameL1U,2)/2 + 128 );
        errFrameL2U_uint8 = uint8( dzquantRecon(errFrameL2U,2)/2 + 128 );
        errFrameV_uint8 = uint8( dzquantRecon(errFrameV,2)/2 + 128 );
        errFrameL1V_uint8 = uint8( dzquantRecon(errFrameL1V,2)/2 + 128 );
        errFrameL2V_uint8 = uint8( dzquantRecon(errFrameL2V,2)/2 + 128 );
        
        % Write the frame back to YUV
        % 1) global layer
        video_dest = [seqName{seqID} '_QP' int2str(QP(QPidx)) '_GlobalME.yuv'];
        yuv_fid = fopen(video_dest, 'a');
        fwrite(yuv_fid, errFrame_uint8', 'uint8');
        fwrite(yuv_fid, errFrameU_uint8', 'uint8');
        fwrite(yuv_fid, errFrameV_uint8', 'uint8');
        fclose(yuv_fid);
        % 2) first object layer
        video_dest = [seqName{seqID} '_QP' int2str(QP(QPidx)) '_ObjL1.yuv'];
        yuv_fid = fopen(video_dest, 'a');
        fwrite(yuv_fid, errFrameL1_uint8', 'uint8');
        fwrite(yuv_fid, errFrameL1U_uint8', 'uint8');
        fwrite(yuv_fid, errFrameL1V_uint8', 'uint8');
        fclose(yuv_fid);
        % 3) second object layer
        video_dest = [seqName{seqID} '_QP' int2str(QP(QPidx)) '_ObjL2.yuv'];
        yuv_fid = fopen(video_dest, 'a');
        fwrite(yuv_fid, errFrameL2_uint8', 'uint8');
        fwrite(yuv_fid, errFrameL2U_uint8', 'uint8');
        fwrite(yuv_fid, errFrameL2V_uint8', 'uint8');
        fclose(yuv_fid);
        
        % That's it, save the data
        save([seqName{seqID} '_QP' int2str(QP(QPidx)) '_ObjVC.mat'], 'tauL0', 'tauL1', 'tauL2', 'qTreeCent',...
            'predFrame','predFrameL1','predFrameL2','errFrame','errFrameL1','errFrameL2',...
            'predFrameU','predFrameL1U','predFrameL2U','errFrameU','errFrameL1U','errFrameL2U',...
            'predFrameV','predFrameL1V','predFrameL2V','errFrameV','errFrameL1V','errFrameL2V');
        % and clear up the parameters
        clear('para','segpara','recurPara');
    end
end
