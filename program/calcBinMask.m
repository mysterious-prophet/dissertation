%% Calculate Binary Mask
% calculates binary mask for normalized image X based on threshold theta

% inputs: X - normalized input image
%       : theta - binary mask threshold
% outputs: bin_mask - binary mask for image X 

function bin_mask = calcBinMask(X, low_val, upp_val, bin_mask_erode, rho, epsilon)
    % binary mask creation
    % volshow(X);

    % figure;
    % montage(X);
    K = calcFilterKernel(X, 'fftDoG', 0, {0.1, 2, 2, 5});
    % X_filt = round(real(filterFunctionFFT(X, K, 0)));
    % X_diff = X - X_filt;
    % 
    % %% edges
    % BW = edge3(X_diff, "approxcanny", 0.3);
    % figure;
    % sliceViewer(BW);
    % 
    % %% grey differences
    % % [M, N, P] = size(X);
    % % W = graydiffweight(X_diff, round(3*M/4), round(N/2), round(P/2), "GrayDifferenceCutoff", 120);
    % % thresh = 0.01;
    % % BW = imsegfmm(W, round(3*M/4), round(N/2), round(P/2), thresh);
    % % figure;
    % % sliceViewer(BW);
    % 
    %% segmentation
    se = strel("sphere", 6);
    % X_open = imopen(X_filt, se);
    % figure;
    % sliceViewer(X_open);
    X_ero = imerode(X, se);
    X_rec = imreconstruct(X_ero, X);
    % figure;
    % sliceViewer(X_rec);
    % X_op_cl = imclose(X_open, se);
    % figure;
    % sliceViewer(X_op_cl);
    X_rec_dil = imdilate(X_rec, se);
    X_rec_dil_rec = imreconstruct(imcomplement(X_rec_dil), imcomplement(X_rec));
    X_rec_dil_rec_com = imcomplement(X_rec_dil_rec);
    % figure;
    % sliceViewer(X_rec_dil_rec);
    % figure;
    % sliceViewer(X_rec_dil_rec_com);
    reg_max = imregionalmax(X_rec_dil_rec_com);
    % figure;
    % sliceViewer(reg_max);
    reg_max_com = imcomplement(reg_max);
    % figure;
    % sliceViewer(reg_max_com);
    se1 = strel("sphere", 10);
    % reg_max_com_cl = imclose(reg_max_com, se1);
    % figure;
    % sliceViewer(reg_max_com_cl);
    reg_max_com_ero = imerode(reg_max_com, se1);
    % figure;
    % sliceViewer(reg_max_com_ero);
    reg_max_final = imfill(imcomplement(bwareaopen(reg_max_com_ero, 1000000)), "holes");
    % figure;
    % sliceViewer(reg_max_final);
    reg_max_final1 = bwareaopen(imopen(reg_max_final, se), 10000);
    % figure;
    % sliceViewer(reg_max_final1);
    reg_max_final2 = imfill(reg_max_final1, "holes");
    % figure;
    % sliceViewer(reg_max_final2);

    bin_mask = reg_max_final2;


    %% result figures
    % figure;
    % montage(reg_max_final2);
    % 
    % fin_vol = volshow(X, OverlayData=reg_max_final2);
    % viewer = fin_vol.Parent;
    % fin_vol.RenderingStyle = "GradientOpacity";
    % fin_vol.Alphamap = linspace(0, 0.2, 256);
    % fin_vol.OverlayAlphamap = 0.6;
    % viewer.BackgroundColor = "white";
    % viewer.BackgroundGradient = "off";


    % this doesn't work
    %% watershed
    % % figure;
    % % sliceViewer(imcomplement(reg_max_com_ero));

    % % bw = imbinarize(X_rec_dil_rec_com, "adaptive");
    % % figure;
    % % sliceViewer(bw);
    % D = bwdist(reg_max_final2);
    % DL = watershed(D, 26);
    % bgm = DL == 0;
    % figure;
    % sliceViewer(bgm);


    % this approach just doesn't work now. Let's try the stuff above and
    % see what happens
    %% initital image filtering prepro
    % K = calcFilterKernel(X, 'fftDoG', 0, {0.1, 2, 2, 5});
    % X_filt = round(real(filterFunctionFFT(X, K, 0)));
    % [M, N, P] = size(X);
    % max_x = max(max(max(X)));
    % min_x = min(min(min(X)));
    % max_x_filt = max(max(max(X_filt)));
    % min_x_filt = min(min(min(X_filt)));
    % 
    % figure;
    % sliceViewer(X);
    % figure;
    % sliceViewer(X_filt);
    % figure;
    % X_diff = X - X_filt;
    % sliceViewer(X_diff)

    %% find border, dilate it and get its inside
    % border = zeros(M, N, P);
    % border(X_filt < 150) = 1;
    % figure;
    % sliceViewer(border);
    % 
    % [~, border_dil] = calcEroDil(border, 1.5);
    % figure;
    % sliceViewer(border_dil);
    % 
    % inside = bwlabeln(border_dil, 6);
    % % in_count = bwconncomp(border_dil, 6);
    % figure;
    % sliceViewer(inside);
    % 
    % [~, inside_dil] = calcEroDil(inside, 1.5);
    % figure;
    % sliceViewer(inside_dil);
    % 
    % % convex hull or concavity/hole filling of inside dil
    % % take bwconnectn of largerst inside cavity and take that and dilate it
    % 
    % %% binary mask creation
    % bin_mask = ones(M, N, P);
    % bin_mask(X_filt < low_val) = 0;
    % bin_mask(X_filt > upp_val) = 0;
    % 
    % % volshow(bin_mask);
    % figure;
    % sliceViewer(bin_mask);
    % 
    % %% binary mask erosion
    % if(bin_mask_erode)
    %     [bin_mask, ~] = calcEroDil(bin_mask, rho); 
    %     figure;
    %     sliceViewer(bin_mask);
    % end
    % 
    % %% binary mask reshaping
    % bin_mask = reshape(bin_mask, [], 1);
    % sum_bin = sum(bin_mask(bin_mask == 1));
    % rel_sum = sum_bin / size(bin_mask, 1);
end

% function [ero_X, dil_X] = calcEroDil(X, rho, epsilon)
%     if nargin < 3 
%         epsilon = 3/(8*pi*rho^3);
%     end
% 
%     % erosion, dilation, create a fftn image of a ball with radius rho
%     [M, N, P] = size(X);
%     X = fftn(X);
% 
%     [U, V, W] = meshgrid((0:N-1), (0:M-1), (0:P-1));
%     omega = sqrt((U/M - 1/2).^2 + (V/N - 1/2).^2 + (W/P - 1/2).^2) + eps;
%     rho_om = rho*omega;
%     taylor_exp = 1 - (rho^2*omega.^2)/10 + (rho^4*omega.^4)/280;
%     rho_om(rho_om <= 1e-3) = taylor_exp(rho_om <= 1e-3);
%     B = 3*(sin(rho_om) - rho_om.*cos(rho_om))./(rho_om.^3);
% 
% 
%     B = ifftshift(B);
%     Y = B .* X;
% 
%     y = real(ifftn(Y));
%     ero_X = (y > 1 - epsilon);
%     dil_X = (y > epsilon);
% 
% end