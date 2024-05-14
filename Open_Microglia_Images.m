function [files, nFiles, fileNo, fileNames, images] = Open_Microglia_Images(path, extension)
    %%%% Oopens a set of numbered microglia images
    %%%%
    %%%% Inputs:
    %%%%    path: path to the directory where the images are stored
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
        k = find(files(n).name == '.') - 1;                                 % find starting index of end of name (before file extension)
        k = k(end);
        %         fileNo(n) = str2double(files(n).name(1:k));
        fileNames = [fileNames convertCharsToStrings(files(n).name(1:k))];  % concatenate file names in a vector
    end
    
    % If no file numbers (strings), set as 1 2 3 4 5...
    if sum(isnan(fileNo(:))) == nFiles
        fileNo = 1:nFiles;
    end

    % Sort everything by file number
    [fileNo, sortIdx] = sort(fileNo);                                       
    files = files(sortIdx);                                                 % reorder files if numbers in name
    fileNames = fileNames(sortIdx);
    
    %% Open images
    for n = 1:nFiles
        images(n).R = imread([path '/' files(n).name]);                     
    end
end
