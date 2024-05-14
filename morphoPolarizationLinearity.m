function [polarizationIndex, linearity] = morphoPolarizationLinearity(image, Gx, Gy, resolution) 
    %%%% Computes parameters based 
    %%%%
    %%%% Inputs:
    %%%%    image: a binary image of a cell
    %%%%    soma: a binary image of the cell's soma
    %%%%    Gx: x coordinate of the cell's soma centroid
    %%%%    Gy: y coordinate of the cell's soma centroid
    %%%%    resolution: resolution of the images, in micrometers per pixel
    %%%%        
    %%%% Outputs:
    %%%%    polarizationIndex: ratio of the distance between gravity
    %%%%                       center of the cell and that of its soma on 
    %%%%                       the square root of the area of the cell
    %%%%    linearity: ratio of max on min variances of coordinates of 
    %%%%               pixels in the cell 


    [rows, cols] = size(image);                                             % size of image

    % Compute coordinates of gravity center of the cell
    ux = (1:cols);
    uy = (1:rows);
    Gx1 = sum(image * ux') / sum(image(:));
    Gx1 = round(Gx1);
    Gy1 = sum(uy * image) / sum(image(:));
    Gy1 = round(Gy1);

    % Offset of the center of the cell's soma and gravity center
    deltaX = Gx - Gx1;
    deltaY = Gy - Gy1;
    vectG = [deltaX, deltaY];
    distance = norm(vectG);
      
    % Polarization of the cell
    polarizationIndex = distance * resolution / ...             
                        sqrt(sum(image(:)) * resolution ^ 2);               
    
    % Linearity of the cell
    [listey, listex] = find(image);
    liste(:, 1) = listey; 
    liste(:, 2) = listex;
    
    % for g = 1:length(listex)                                                % revert the x axis
    %     liste(g, 2) = abs(listex(g, :) - cols);
    % end
    
    % PCA on the coordinates of the cell image
    [scoresPCA, ~, ~, ~, ~, ~] = PCA_custom(liste);
    
    variance(:, 1) = var(scoresPCA(:, 1));                                  % variance of PC1
    variance(:, 2) = var(scoresPCA(:, 2));                                  % variance of PC2
    linearity = max(variance(:)) / min(variance(:));                        % compute linearity
end
