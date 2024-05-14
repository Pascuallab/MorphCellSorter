function DisplayImages(path, extension, colors, ranking, classe)
    %%%% Display_Images displays 64 images per figure with ranking or not
    %%%% and cluster colors or not
    %%%%
    %%%% Inputs:
    %%%%    path: path to the images to display
    %%%%    ranking : array of strings of the ranked names of the microglia
    %%%%        (optional, set to [] if not needed)
    %%%%    clusters : array of numbers of the clusters of the microglia
    %%%%        (optional, set to [] if not needed)
    %%%%    colors : colormap for clusters
    %%%%        (optional, set to [] if not needed, needed if clusters is set)

    %% Get binary image files
    if isempty(ranking)
        [~, ~, ~, ~, images] = Open_Microglia_Images(path, extension);
    else
        [~, ~, ~, ~, images] ...
            = Open_Microglia_Images_With_Ranking(path, extension, ranking);
    end
    % Get image sizes
    
    [imWidth, imHeight] = size(images(1).R);
    
    % Zoom on 1/3rd at the center and color
    for n = 1:length(images)
        if ~isempty(classe)
            idxWhite = find(images(n).R);
            redImg = 255 * zeros(imWidth, imHeight, 'uint8');
            greenImg = 255 * zeros(imWidth, imHeight, 'uint8');
            blueImg = 255 * zeros(imWidth, imHeight, 'uint8');
            
            redImg(idxWhite) = colors(classe(n), 1) * 255;
            greenImg(idxWhite) = colors(classe(n), 2) * 255;
            blueImg(idxWhite) = colors(classe(n), 3) * 255;
            images(n).R = cat(3, redImg, greenImg, blueImg);
        end
    end
% %     Display images
    if length(images)<352
        nbCols=22;
        nbRows=round(length(images)/22)+1;
         nbRows=round(347/22)+1;
            figure()
            set(gcf, 'Position', get(0, 'Screensize'));
            montage({images(1:length(images)).R}, 'Size', [nbRows nbCols]);
    else
        for n=1:352:length(images)
            if n+351<=length(images)
                figure()
                montage({images(n:n+351).R},'Size',[16 22]);
            else
                figure;
                montage({images(n:end).R}, 'Size', [11 22]);
            end
        end
    end



end 
