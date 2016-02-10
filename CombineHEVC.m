clear all;

% seqName = {'Stefan', 'Office1', 'City', 'Football'};
seqName = {'Stefan', 'Office1', 'City'};
QP = [22 27 32 37];
QP = 22;

% for seqID = 1:length(seqName),
%     for QPidx = 1:length(QP),
%         currPath = ['C:\Users\akidu\Dropbox\ObjVCData\' seqName{seqID} '\'];
%         currData = [currPath seqName{seqID} '_QP' int2str(QP(QPidx)) '_ObjVC.mat'];
%         load(currData);
% %         re-cast the error frame back to uint8 for HEVC coding
%         errFrame_uint8 = uint8( dzquantRecon(errFrame,2)/2 + 128 );
%         errFrameL1_uint8 = uint8( dzquantRecon(errFrameL1,2)/2 + 128 );
%         errFrameL2_uint8 = uint8( dzquantRecon(errFrameL2,2)/2 + 128 );
%         errFrameU_uint8 = uint8( dzquantRecon(errFrameU,2)/2 + 128 );
%         errFrameL1U_uint8 = uint8( dzquantRecon(errFrameL1U,2)/2 + 128 );
%         errFrameL2U_uint8 = uint8( dzquantRecon(errFrameL2U,2)/2 + 128 );
%         errFrameV_uint8 = uint8( dzquantRecon(errFrameV,2)/2 + 128 );
%         errFrameL1V_uint8 = uint8( dzquantRecon(errFrameL1V,2)/2 + 128 );
%         errFrameL2V_uint8 = uint8( dzquantRecon(errFrameL2V,2)/2 + 128 );
%         
%         % just shift version
%         errFrame_uint16 = uint16(floor( errFrame + 255 ));
%         errFrameL1_uint16 = uint16(floor( errFrameL1 + 255 ));
%         errFrameL2_uint16 = uint16(floor( errFrameL2 + 255 ));
%         errFrameU_uint16 = uint16(floor( errFrameU + 255 ));
%         errFrameL1U_uint16 = uint16(floor( errFrameL1U + 255 ));
%         errFrameL2U_uint16 = uint16(floor( errFrameL2U + 255 ));
%         errFrameV_uint16 = uint16(floor( errFrameV + 255 ));
%         errFrameL1V_uint16 = uint16(floor( errFrameL1V + 255 ));
%         errFrameL2V_uint16 = uint16(floor( errFrameL2V + 255 ));
%         errFrame_uint16(errFrame_uint16>511) = 511;
%         errFrameL1_uint16(errFrameL1_uint16>511) = 511;
%         errFrameL2_uint16(errFrameL2_uint16>511) = 511;
%         errFrameU_uint16(errFrameU_uint16>511) = 511;
%         errFrameL1U_uint16(errFrameL1U_uint16>511) = 511;
%         errFrameL2U_uint16(errFrameL2U_uint16>511) = 511;
%         errFrameV_uint16(errFrameV_uint16>511) = 511;
%         errFrameL1V_uint16(errFrameL1V_uint16>511) = 511;
%         errFrameL2V_uint16(errFrameL2V_uint16>511) = 511;
%         
% %         Write the frame back to YUV
% %         1) global layer
%         video_dest = [currPath '9bits\' seqName{seqID} '_QP' int2str(QP(QPidx)) '_GlobalME.yuv'];
%         yuv_fid = fopen(video_dest, 'a');
%         fwrite(yuv_fid, errFrame_uint16', 'uint16');
%         fwrite(yuv_fid, errFrameU_uint16', 'uint16');
%         fwrite(yuv_fid, errFrameV_uint16', 'uint16');
%         fclose(yuv_fid);
% %         2) first object layer
%         video_dest = [currPath '9bits\' seqName{seqID} '_QP' int2str(QP(QPidx)) '_ObjL1.yuv'];
%         yuv_fid = fopen(video_dest, 'a');
%         fwrite(yuv_fid, errFrameL1_uint16', 'uint16');
%         fwrite(yuv_fid, errFrameL1U_uint16', 'uint16');
%         fwrite(yuv_fid, errFrameL1V_uint16', 'uint16');
%         fclose(yuv_fid);
% %         3) second object layer
%         video_dest = [currPath '9bits\' seqName{seqID} '_QP' int2str(QP(QPidx)) '_ObjL2.yuv'];
%         yuv_fid = fopen(video_dest, 'a');
%         fwrite(yuv_fid, errFrameL2_uint16', 'uint16');
%         fwrite(yuv_fid, errFrameL2U_uint16', 'uint16');
%         fwrite(yuv_fid, errFrameL2V_uint16', 'uint16');
%         fclose(yuv_fid);
%     end
% end
%% PSNR calculation
allQP = QP;
DecPath = 'C:\Users\akidu\Dropbox\ObjVCData\Decoded\9 bit\QP';
OriPath = 'C:\Users\akidu\Documents\ObjectVideoCoding\Data\';
HEVCPath = 'C:\Users\akidu\Downloads\ObjVC Raw Data\';
height = [288 480 288];
width = [352 640 352];
Hfil = [2/3; 1/6; 1/6];
for seqIdx = 1:length(seqName),
    h = height(seqIdx); w = width(seqIdx);
    for qpIdx = 1:length(allQP),
        QP = allQP(qpIdx);
        % data
        currPath = ['C:\Users\akidu\Dropbox\ObjVCData\' seqName{seqIdx} '\'];
        currSeq = [currPath seqName{seqIdx} '_QP' int2str(QP) '_ObjVC.mat'];
        currOri = [OriPath seqName{seqIdx} '_QP' int2str(QP) '.mat'];
        currHEVC = [HEVCPath seqName{seqIdx} '\' seqName{seqIdx} '_dec_QP' int2str(QP) '.yuv'];
        load(currSeq); load(currOri,'Ymtx_Inc','Umtx_Inc','Vmtx_Inc');
        currL1Dec = [DecPath int2str(QP) '\dec_' seqName{seqIdx} '_QP' int2str(QP) '_ObjL1.yuv'];
        currL2Dec = [DecPath int2str(QP) '\dec_' seqName{seqIdx} '_QP' int2str(QP) '_ObjL2.yuv'];
        
        [HEVCrecon, HEVCreconU, HEVCreconV] = yuv_read(currHEVC,8,w,h,2,'yv12');
        
        [errL1_recon, errL1_reconU, errL1_reconV] = yuv_read(currL1Dec,16,w,h,0,'yv12');
        [errL2_recon, errL2_reconU, errL2_reconV] = yuv_read(currL2Dec,16,w,h,0,'yv12');
        % convert back
        errL1_recon_final = (double(errL1_recon) - 255);
        errL1_reconU_final = (double(errL1_reconU) - 255);
        errL1_reconV_final = (double(errL1_reconV) - 255);
        errL2_recon_final = (double(errL2_recon) - 255);
        errL2_reconU_final = (double(errL2_reconU) - 255);
        errL2_reconV_final = (double(errL2_reconV) - 255);
        
        disp(norm(errL2_recon_final(:),1));
        disp(norm(errFrameL2(:),1));
        disp(norm(errL2_recon_final(:)));
        disp(norm(errFrameL2(:)));
        
        % reconstructed frames
        L1recon_final = errL1_recon_final + predFrameL1;
        L1reconU_final = errL1_reconU_final + predFrameL1U;
        L1reconV_final = errL1_reconV_final + predFrameL1V;
        L2recon_final = errL2_recon_final + predFrameL2;
        L2reconU_final = errL2_reconU_final + predFrameL2U;
        L2reconV_final = errL2_reconV_final + predFrameL2V;
        
        % MSE and PSNR
        MSE_L1y = mean((L1recon_final(:) - double(Ymtx_Inc(:))).^2);
        MSE_L1u = mean((L1reconU_final(:) - double(Umtx_Inc(:))).^2);
        MSE_L1v = mean((L1reconV_final(:) - double(Vmtx_Inc(:))).^2);
        MSE_L1(qpIdx,seqIdx) = [MSE_L1y MSE_L1u MSE_L1v]*Hfil;
        
        MSE_L2y = mean((L2recon_final(:) - double(Ymtx_Inc(:))).^2);
        MSE_L2u = mean((L2reconU_final(:) - double(Umtx_Inc(:))).^2);
        MSE_L2v = mean((L2reconV_final(:) - double(Vmtx_Inc(:))).^2);
        MSE_L2(qpIdx,seqIdx) = [MSE_L2y MSE_L2u MSE_L2v]*Hfil;
        
        PSNR_L1(qpIdx,seqIdx) = 10*log10(255^2./MSE_L1(qpIdx,seqIdx));
        PSNR_L2(qpIdx,seqIdx) = 10*log10(255^2./MSE_L2(qpIdx,seqIdx));
        
        % Y SSIM
        SSIM_L1(qpIdx,seqIdx) = ssim(L1recon_final,double(Ymtx_Inc));
        SSIM_L2(qpIdx,seqIdx) = ssim(L2recon_final,double(Ymtx_Inc));
        SSIM_HEVC(qpIdx,seqIdx) = ssim(double(HEVCrecon),double(Ymtx_Inc));
    end
end