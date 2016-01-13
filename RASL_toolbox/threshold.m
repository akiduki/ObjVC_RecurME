function [xhat1, xhat2] = threshold(y1, y2, prox)
% THRESHOLD implements 2D threshold.
%
% This function solves [m x n] problems of type:
% xhat = ARG_MIN{ 0.5*norm(y - x)^2  + lambda*norm(x) } where minimization
% is over x in R^2.
%
% Input:
% - y1, y2: [m x n] matrices of measurements
% - lambda: scalar representing the threshold
%
% Output:
% - xhat1, xhat2: [m x n] matrices of solutions.
%
% Ulugbek Kamilov, Emrah Bostan, BIG @ EPFL, 2011.

% Dimensions of the input
[nr, nc] = size(y1);

% Norm of the input
ny = sqrt(y1.^2 + y2.^2);

% Prevent division by zero
ny(ny <= 0) = 1;

% Compute threshold
% t = prox(ny(:)) ./ ny(:);
t = max(ny(:)-prox,0) ./ ny(:);

xhat1 = reshape(t .* y1(:), nr, nc);
xhat2 = reshape(t .* y2(:), nr, nc);