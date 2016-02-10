function y = iht(w, dir, shift)
% IHT 1D single scale inverse haar-transform for 2D signals.
%
% Input:
% - w: transformed image [n x n] where half of the pixels are approximation
% coefficients and the rest are detail coefficients.
% - dir: string representing the direction of transform ['h' or 'v'].
% - shift: if non-zero shifts the output before returning.
%
% Output:
% - y: reconstructed image [n x n]
%
% Ulugbek Kamilov, Emrah Bostan, BIG @ EPFL, 2011.

% Dimensions of the input signal
[nr, nc] = size(w);

% Initialize wavelet domain output
y = zeros(nr, nc);

% Haar normalization constant
C = 1/realsqrt(2);

% Process
switch(dir)
    case{'h'}
        % Mid-point
        m = nc/2;
        
        y(:, 1:2:end) = C*(w(:, 1:m) - w(:, m+1:end));
        y(:, 2:2:end) = C*(w(:, 1:m) + w(:, m+1:end));
        
        % Shift
        if(shift)
            y = [y(:, end), y(:, 1:end-1)];
        end
        
    case{'v'}
        % Mid-point
        m = nr/2;
        
        y(1:2:end, :) = C*(w(1:m, :) - w(m+1:end, :));
        y(2:2:end, :) = C*(w(1:m, :) + w(m+1:end, :));
        
        % Shift
        if(shift)
            y = [y(end, :); y(1:end-1, :)];
        end
        
    otherwise
        error('Only horizontal and verical shifts supported');
end

