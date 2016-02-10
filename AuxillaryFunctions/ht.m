function w = ht(x, dir, shift)
% HT 1D single scale haar-transform for a 2D signal x.
%
% Input:
% - x: image of size [n x n] where n must be even.
% - dir: string representing the direction of transform ['h' or 'v'].
% - shift: if non-zero shifts the input before taking the transform.
%
% Output:
% - w: transformed image [n x n] where half of pixels are approximation
% coefficients and the rest are detail coefficients.
%
% Ulugbek Kamilov, Emrah Bostan, BIG @ EPFL, 2011.


% Dimensions of the input signal
[nr, nc] = size(x);

% Initialize wavelet domain output
w = zeros(nr, nc);

% Haar-normalization
C = 1/realsqrt(2);

% Process
switch(dir)
    case{'h'}
        % Shift
        if(shift)
            x = [x(:, 2:end), x(:, 1)];            
        end
        
        % Mid-point
        m = nc/2;
        
        w(:, 1:m) = C*(x(:, 2:2:end) + x(:, 1:2:end));
        w(:, m+1:end) = C*(x(:, 2:2:end) - x(:, 1:2:end));
        
    case{'v'}
        % Shift
        if(shift)
            x = [x(2:end, :); x(1, :)];
        end
        
        % Mid-point
        m = nr/2;
        
        w(1:m, :) = C*(x(2:2:end, :) + x(1:2:end, :));
        w(m+1:end, :) = C*(x(2:2:end, :) - x(1:2:end, :));
        
    otherwise
        error('Only horizontal and verical shifts supported');
end

