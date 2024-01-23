%% Calculate Binary Mask
% calculates binary mask for normalized image X based on threshold theta

% inputs: X - normalized input image
%       : theta - binary mask threshold
% outputs: bin_mask - binary mask for image X 

function bin_mask = calcBinMask(X, low_val, upp_val, bin_mask_erode, rho, epsilon)
    % binary mask creation
    % volshow(X);
    K = calcFilterKernel(X, 'fftDoB', 0, 0.25, 1);
    X_filt = X .* K;
    [M, N, P] = size(X);
    bin_mask = ones(M, N, P);
    bin_mask(X < low_val) = 0;
    bin_mask(X > upp_val) = 0;

    % volshow(bin_mask);
    % figure;
    % sliceViewer(bin_mask);

    if(bin_mask_erode)
        if nargin == 5
            epsilon = 3/(8*pi*rho^3);
        end

        % binary mask erosion, create a fftn image of a ball with radius rho
        bin_mask = fftn(bin_mask);
    
        [U, V, W] = meshgrid((0:N-1), (0:M-1), (0:P-1));
        omega = sqrt((U/M - 1/2).^2 + (V/N - 1/2).^2 + (W/P - 1/2).^2) + eps;
        rho_om = rho*omega;
        taylor_exp = 1 - (rho^2*omega.^2)/10 + (rho^4*omega.^4)/280;
        rho_om(rho_om <= 1e-3) = taylor_exp(rho_om <= 1e-3);
        B = 3*(sin(rho_om) - rho_om.*cos(rho_om))./(rho_om.^3);
    
    
        B = ifftshift(B);
        Y = B .* bin_mask;

        y = real(ifftn(Y));
        % y_2 = reshape(y, [], 1);
        % max_y = max(y_2);
        % min_y = min(y_2);
        bin_mask = (y > 1 - epsilon);
        % volshow(bin_mask);
        % figure;
        % sliceViewer(bin_mask);
    end

    % binary mask reshaping
    bin_mask = reshape(bin_mask, [], 1);
    % sum_bin = sum(bin_mask(bin_mask == 1));
    % rel_sum = sum_bin / size(bin_mask, 1);
end