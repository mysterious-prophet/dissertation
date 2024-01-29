%% Calculate Binary Mask
% calculates binary mask for normalized image X based on threshold theta

% inputs: X - normalized input image
%       : theta - binary mask threshold
% outputs: bin_mask - binary mask for image X 

function bin_mask = calcBinMask(X, low_val, upp_val, bin_mask_erode, rho, epsilon)
    % binary mask creation
    % volshow(X);

    %% initital image filtering prepro
    K = calcFilterKernel(X, 'fftLP', 'buttKer', {1, 2});
    X_filt = round(real(filterFunctionFFT(X, K, 0)));
    [M, N, P] = size(X);
    max_x = max(max(max(X)));
    min_x = min(min(min(X)));
    max_x_filt = max(max(max(X_filt)));
    min_x_filt = min(min(min(X_filt)));

    figure;
    sliceViewer(X);
    figure;
    sliceViewer(X_filt);
    figure;
    X_diff = X - X_filt;
    sliceViewer(X_diff);

    %% find border, dilate it and get its inside
    border = zeros(M, N, P);
    border(X_filt < 250) = 1;
    figure;
    sliceViewer(border);

    [~, border_dil] = calcEroDil(border, 2);
    figure;
    sliceViewer(border_dil);

    inside = bwlabeln(border_dil, 18);
    % in_count = bwconncomp(border_dil, 6);
    figure;
    sliceViewer(inside);

    [~, inside_dil] = calcEroDil(inside, 1.5);
    figure;
    sliceViewer(inside_dil);

    %% binary mask creation
    bin_mask = ones(M, N, P);
    bin_mask(X_filt < low_val) = 0;
    bin_mask(X_filt > upp_val) = 0;

    % volshow(bin_mask);
    figure;
    sliceViewer(bin_mask);

    %% binary mask erosion
    if(bin_mask_erode)
        [bin_mask, ~] = calcEroDil(bin_mask, rho); 
        figure;
        sliceViewer(bin_mask);
    end

    %% binary mask reshaping
    bin_mask = reshape(bin_mask, [], 1);
    sum_bin = sum(bin_mask(bin_mask == 1));
    rel_sum = sum_bin / size(bin_mask, 1);
end

function [ero_X, dil_X] = calcEroDil(X, rho, epsilon)
    if nargin < 3 
        epsilon = 3/(8*pi*rho^3);
    end

    % erosion, dilation, create a fftn image of a ball with radius rho
    [M, N, P] = size(X);
    X = fftn(X);

    [U, V, W] = meshgrid((0:N-1), (0:M-1), (0:P-1));
    omega = sqrt((U/M - 1/2).^2 + (V/N - 1/2).^2 + (W/P - 1/2).^2) + eps;
    rho_om = rho*omega;
    taylor_exp = 1 - (rho^2*omega.^2)/10 + (rho^4*omega.^4)/280;
    rho_om(rho_om <= 1e-3) = taylor_exp(rho_om <= 1e-3);
    B = 3*(sin(rho_om) - rho_om.*cos(rho_om))./(rho_om.^3);


    B = ifftshift(B);
    Y = B .* X;

    y = real(ifftn(Y));
    ero_X = (y > 1 - epsilon);
    dil_X = (y > epsilon);

end