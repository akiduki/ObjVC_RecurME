function alpha = wavetr(x)
% WAVETR overcomplete representation of the image.
%
% alpha = wavetr(x)
%
% Input:
% - x: signal
%
% Output:
% - alpha: coefficients
%
% U. S. Kamilov, E. Bostan, BIG, EPFL, 2012.

C = 1/realsqrt(2);

xh = [x(:, 2:end), x(:, 1)];
xv = [x(2:end, :); x(1, :)];

avgh = C*(xh + x);
avgv = C*(xv + x);

difh = C*(xh - x);
difv = C*(xv - x);

alpha = [avgh; avgv; difh; difv];