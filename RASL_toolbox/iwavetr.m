function x = iwavetr(alpha)
% IWAVETR applies dual operator H* to invert the overcomplete signal.
% 
% x = iwavetr(alpha)
%
% Input:
% - alpha: coefficients
%
% Output:
% - x: reconstruction
%
% U. S. Kamilov, E. Bostan, BIG, EPFL, 2011.

% size
[nr, nc] = size(alpha);
nr = 0.25*nr;


% Get individual coefs
w1 = zeros(nr, nc);
w2 = zeros(nr, nc);
w3 = zeros(nr, nc);
w4 = zeros(nr, nc);

m = 0.5*nc;

% Put coeffs
w1(:, 1:m) = alpha(1:nr, 1:2:end);
w1(:, m+1:end) = alpha(2*nr+1:3*nr, 1:2:end);
w1 = iht(w1, 'h', 0);

w2(:, 1:m) = alpha(1:nr, 2:2:end);
w2(:, m+1:end) = alpha(2*nr+1:3*nr, 2:2:end);
w2 = iht(w2, 'h', 1);

m = 0.5*nr;

w3(1:m, :) = alpha(nr+1:2:2*nr, :); 
w3(m+1:end, :) = alpha(3*nr+1:2:end, :); 
w3 = iht(w3, 'v', 0);

w4(1:m, :) = alpha(nr+2:2:2*nr, :);
w4(m+1:end, :) = alpha(3*nr+2:2:end, :);
w4 = iht(w4, 'v', 1);

% Average
x = 0.25*(w1 + w2 + w3 + w4);