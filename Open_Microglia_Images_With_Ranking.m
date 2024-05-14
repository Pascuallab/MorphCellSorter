function [files, nFiles, fileNo, fileNames, images] = Open_Microglia_Images_With_Ranking(path, extension, ranking)
    %%%% Open_Microglia_Images opens a set of numbered microglia images
    %%%%
    %%%% Inputs:
    %%%%    path: path to the directory where the images are stored
    %%%%    ranking: an array containing image names in ranking order
    %%%%
    %%%% Outputs:
    %%%%    files: informations on the files contained at the path
    %%%%    nFiles: number of image files
    %%%%    fileNo: an array containing the number/id of each image
    %%%%    images: an array containing the image files
    
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
    
    %% Open images
    for n = 1:size(files, 1)
        images(n).R = imread([path '/' files(n).name]);
    end
    
    nFiles = length(fileNames);
end
