%% Frequency domain kernel calculation
% calculates kernel to be used in frequency domain filtering - ideal
% circular kernel, circular Bessel kernel, Bessel polynomial kernel,
% Gaussian kernel, Butterworth kernel, alpha-stable kernel. Also calculates
% corresponding highpass variants and laplacian frequency kernel

% inputs: X - input image
%       : filter_name - fftLP, fftHP, fftLap
%                     - name of filter used - for lowpass and highpass
%                     - distinction and laplacian filter
%       : kernel_style - idealKer, circBessKer, bessKer, gaussKer, buttKer,
%                      - alpStKer 
%                      - kernel style used
%       : varargin - various kernel parameters
% outputs: K - calculated frequency domain kernel

function K = calcFilterKernel(X, filter_name, kernel_style, varargin)
    dims_X = ndims(X);
    if(dims_X == 3)
        [m, n, p] = size(X);
        M = 2*m; N = 2*n; P = 2*p;
        [u, v, t] = meshgrid(0:N-1, 0:M-1, 0:P-1);
        d = sqrt((u/M - 1/2).^2 + (v/N - 1/2).^2 + (t/P - 1/2).^2);
    end
    d_sq = d.^2;

    if(strcmp(filter_name, 'fftLP'))
        switch kernel_style
            case 'binKer'
                if(numel(varargin{1}) < 1)
                    dist = 0.3;
                else
                    dist = varargin{1}{1};
                end

                K = ones(M, N);
                K(d > dist) = 0;
 
            case 'gaussKer'
                if(numel(varargin{1}) < 1)
                    sigma = 0.25;
                else
                    sigma = varargin{1}{1};
                end
    
                K = exp(-(d_sq)/(2*sigma^2));

            case 'alpStKer'
                if(numel(varargin{1}) < 2)
                    dist = 0.25;
                    alpha = 2;
                else
                    dist = varargin{1}{1};
                    alpha = varargin{1}{2};
                end

                K = exp(-(d / dist).^(alpha));
    
            case 'buttKer'
                if(numel(varargin{1}) < 2)
                    dist = 0.25;
                    r = 1;
                else
                    dist = varargin{1}{1};
                    r = varargin{1}{2};
                end

                K = 1 ./ sqrt((1 + (d / dist).^(2*r)));
                
        end
    elseif(strcmp(filter_name, 'fftHP'))
        K = ones(M, N) - calcFilterKernel(X, 'fftLP', kernel_style, varargin{1});
    elseif(strcmp(filter_name, 'fftLap'))
        K = -4*pi^2*d_sq;
    elseif(strcmp(filter_name, 'fftLoB'))
        K = d_sq .* calcFilterKernel(X, 'fftLP', 'buttKer', varargin);
    elseif(strcmp(filter_name, 'fftLoG'))
        K = d_sq .* calcFilterKernel(X, 'fftLP', 'gaussKer', varargin);
    elseif(strcmp(filter_name, 'fftDoB'))
        if(numel(varargin) < 4)
            dist1 = 0.5;
            r1 = 1;
            dist2 = 0.25;
            r2 = 1;
        else
            dist1 = varargin{1, 1};
            r1 = varargin{1, 2};
            dist2 = varargin{1, 3};
            r2 = varargin{1, 4};
        end
        K = calcFilterKernel(X, 'fftLP','buttKer', {dist1, r1}) ...
            - calcFilterKernel(X, 'fftLP', 'buttKer', {dist2, r2});
    elseif(strcmp(filter_name, 'fftDoG'))
        if(numel(varargin{1}) < 2)
            sigma1 = 5;
            sigma2 = 1;
        else
            sigma1 = varargin{1}{1};
            sigma2 = varargin{1}{2};
        end
        K = calcFilterKernel(X, 'fftLP', 'gaussKer', {sigma1}) ...
            - calcFilterKernel(X, 'fftLP', 'gaussKer', {sigma2});

    end
end