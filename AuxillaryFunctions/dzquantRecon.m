function reconQBlk = dzquantRecon(currCoef,qfact)
% Dead-zone quantizer, derived from quantRecon
% Author: Yuanyi Xue @ NYU-Poly
% Date: 01-13-15
% This function put an uniform quantization to the meanBlk and coefficients
% returns a reconstructed block after quantization
% 1) currCoef: currrent coefficients for the current block
% 2) meanBlk: mean value of the current block
% WARNING: NOT USED IN CURRENT VER.
% 3) qStep: quantization stepsize
% 4) mode: 1=only quantization, 0=quantize and reconstruct
% mode=3 used for the debug
% Outputs:
% a) reconQBlk: the reconstructed block after quantization if mode=0; or
% the quantized value (index) if mode=1
% this version uses the imquantize function

if sum(isnan(currCoef(:)))==0,
    % quantized index
    qIdx = floor(abs(currCoef)./qfact + 0.5).*sign(currCoef);
    % reconstruction levels
    reconQBlk = qfact.*qIdx;
else
    error('NaN happens, check the input!!');
end
% qIdx(find(qIdx>length(dLev))) = length(dLev);
% if mode == 0,
%     reconQBlk = (dLev(qIdx)+dLev(qIdx-1))/2;
% elseif mode == 1,
%     reconQBlk = qIdx;
% end