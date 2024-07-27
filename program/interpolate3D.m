function interpolation = interpolate3D(input_image, dx, dy, dz, spacing)
    % original size
    [m, n, p] = size(input_image);
    
    % original coordinates
    x = (0:n-1)*dx;
    y = (0:m-1)*dy;
    z = (0:p-1)*dz;
    
    % new dimensions diminished by one
    N = floor(dx/spacing*(n-1));
    M = floor(dy/spacing*(m-1));
    P = floor(dz/spacing*(p-1));
    
    % cubic interpolation grid of given positive spacing
    [XQ, YQ, ZQ] = meshgrid(0:N, 0:M, 0:P);
    
    % linear 3D interpolation
    interpolation = interp3(x, y, z, input_image, XQ*spacing, YQ*spacing, ZQ*spacing, 'linear');