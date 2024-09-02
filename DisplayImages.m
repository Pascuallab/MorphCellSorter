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
    % if isempty(ranking)
    %     [~, ~, ~, ~, images] = Open_Microglia_Images(path, extension);
    % else
    %     [~, ~, ~, ~, images] ...
    %         = Open_Microglia_Images_With_Ranking(path, extension, ranking);
    % end

    files = dir([path '/' extension]);
    nFiles = length(files);
    fileNames = [];

    %% Sort files by number
    fileNo = zeros(1, nFiles);
    
    % Loop: get file numbers and file names
    for n = 1:nFiles
        k = find(files(n).name == '.') - 1;
        k = k(end);
        fileNo(n) = str2double(files(n).name(1:k));
        fileNames = [fileNames convertCharsToStrings(files(n).name(1:k))];
    end
    
    % If no file numbers (strings), set as 1 2 3 4 5...
    if sum(isnan(fileNo(:))) == nFiles
        fileNo = 1:nFiles;
    end

    %% Sort files by ranking
    idxRanking = zeros(1, nFiles);

    for n = 1:size(ranking, 1)
        idxRanking(n) = find(ranking(n) == fileNames);
    end

    fileNo = fileNo(nonzeros(idxRanking));
    files = files(nonzeros(idxRanking));
    fileNames = fileNames(nonzeros(idxRanking));

    % Get image sizes
    
    % [imWidth, imHeight] = size(images(1).R);
    
    % Zoom on 1/3rd at the center and color
    % for n = 1:length(images)
    %     if ~isempty(classe)
    %         idxWhite = find(images(n).R);
    %         redImg = 255 * zeros(imWidth, imHeight, 'uint8');
    %         greenImg = 255 * zeros(imWidth, imHeight, 'uint8');
    %         blueImg = 255 * zeros(imWidth, imHeight, 'uint8');
    % 
    %         redImg(idxWhite) = colors(classe(n), 1) * 255;
    %         greenImg(idxWhite) = colors(classe(n), 2) * 255;
    %         blueImg(idxWhite) = colors(classe(n), 3) * 255;
    %         images(n).R = cat(3, redImg, greenImg, blueImg);
    %     end
    % end

% %     Display images
    % if length(images)<352
    %     nbCols=22;
    %     nbRows=round(length(images)/22)+1;
    %      nbRows=round(347/22)+1;
    %         figure()
    %         set(gcf, 'Position', get(0, 'Screensize'));
    %         montage({images(1:length(images)).R}, 'Size', [nbRows nbCols]);
    % else
    %     for n=1:352:length(images)
    %         if n+351<=length(images)
    %             figure()
    %             montage({images(n:n+351).R},'Size',[16 22]);
    %         else
    %             figure;
    %             montage({images(n:end).R}, 'Size', [11 22]);
    %         end
    %     end
    % end

    % n = length(images);
    size_montage = 1200;
    if nFiles < size_montage
        for k = 1:nFiles
            images(k).R = imread([path '/' files(k).name]);
        end
        [imWidth, imHeight] = size(images(1).R);
        for k = 1:length(images)
            if ~isempty(classe)
                idxWhite = find(images(k).R);
                redImg = 255 * zeros(imWidth, imHeight, 'uint8');
                greenImg = 255 * zeros(imWidth, imHeight, 'uint8');
                blueImg = 255 * zeros(imWidth, imHeight, 'uint8');
    
                redImg(idxWhite) = colors(classe(k), 1) * 255;
                greenImg(idxWhite) = colors(classe(k), 2) * 255;
                blueImg(idxWhite) = colors(classe(k), 3) * 255;
                images(k).R = cat(3, redImg, greenImg, blueImg);
            end
        end
        x = ceil(sqrt(nFiles/12));                                          % desired montage ratio 4/3 -> 4 * 3 = 12
        figure()
        set(gcf, 'Position', get(0, 'Screensize'));
        montage({images(1:nFiles).R}, 'Size', [3*x 4*x])
    else 
        for k=1:size_montage:nFiles
            if k-1+size_montage<nFiles
                for l = 1:size_montage
                    images(l).R = imread([path '/' files(l+k-1).name]);
                end
                [imWidth, imHeight] = size(images(1).R);
                for n = 1:length(images)
                    if ~isempty(classe)
                        idxWhite = find(images(n).R);
                        redImg = 255 * zeros(imWidth, imHeight, 'uint8');
                        greenImg = 255 * zeros(imWidth, imHeight, 'uint8');
                        blueImg = 255 * zeros(imWidth, imHeight, 'uint8');
                
                        redImg(idxWhite) = colors(classe(n+k-1), 1) * 255;
                        greenImg(idxWhite) = colors(classe(n+k-1), 2) * 255;
                        blueImg(idxWhite) = colors(classe(n+k-1), 3) * 255;
                        images(n).R = cat(3, redImg, greenImg, blueImg);
                    end
                end
                x = ceil(sqrt(length(images)/12));                                                   % desired montage ratio 4/3 -> 4 * 3 = 12
                a = figure();
                set(gcf, 'Position', get(0, 'Screensize'));
                montage({images(1:length(images)).R}, 'Size', [3*x 4*x])
                saveas(a, fullfile(path, ['\Montage_' num2str(k) '.png']));
                close(a);
            else
                nn = nFiles;
                for l=1:nn-k
                    images(l).R = imread([path '/' files(l+k-1).name]);
                end
                [imWidth, imHeight] = size(images(1).R);
                for n = 1:length(images)
                    if ~isempty(classe)
                        idxWhite = find(images(n).R);
                        redImg = 255 * zeros(imWidth, imHeight, 'uint8');
                        greenImg = 255 * zeros(imWidth, imHeight, 'uint8');
                        blueImg = 255 * zeros(imWidth, imHeight, 'uint8');
                
                        redImg(idxWhite) = colors(classe(n+k-1), 1) * 255;
                        greenImg(idxWhite) = colors(classe(n+k-1), 2) * 255;
                        blueImg(idxWhite) = colors(classe(n+k-1), 3) * 255;
                        images(n).R = cat(3, redImg, greenImg, blueImg);
                    end
                end
                x = ceil(sqrt(length(images)/12));                                                   % desired montage ratio 4/3 -> 4 * 3 = 12
                a = figure();
                set(gcf, 'Position', get(0, 'Screensize'));
                montage({images(1:length(images)).R}, 'Size', [3*x 4*x])
                saveas(a, fullfile(path, ['\Montage_' num2str(k) '.png']));
                close(a);
            end
        end
    end

end 
