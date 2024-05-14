function roundFactor = roundnessFactor(image, somaArea, resolution)
    %%%% Computes the roundness factor of image
    %%%%
    %%%% Inputs:
    %%%%    image: a binary image of a cell
    %%%%    somaArea: area of the soma of the cell in square micrometers
    %%%%    resolution: resolution of the image, in micrometers per pixel
    %%%%        
    %%%% Outputs:
    %%%%    roundFactor: roundness factor of the cell in image

    [y, x] = find(image);                                                   % coordinates of pixels in image
    y = y - mean(y);                                                        % centering coordinates
    x = x - mean(x);
    M = [y x]' * [y x];                                                     % computing the covariance matrix of coordinates
    [V, ~] = eig(M);                                                        % computing the right eigenvectors of the covariance matrix
    pXY = [y x] * V;                                                        % scaling coordinates by their eigenvalues
    axisLength = 2 * std(pXY);                                              % axis length as twice the standard deviation of 
    roundFactor = 4 * somaArea / (pi * (max(axisLength) * resolution) ^ 2); % computing roundness factor
end
