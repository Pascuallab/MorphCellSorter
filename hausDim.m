function D = hausDim(I)
    %%%% Returns the Haussdorf fractal dimension of an object represented by
    %%%% a binary image.
    %%%%    Returns the Haussdorf fractal dimension D of an object represented by the
    %%%%    binary image I. Nonzero pixels belong to an object and 0 pixels 
    %%%%    constitute the background.
    %%%%
    %%%%    Algorithm
    %%%%    ---------
    %%%%    1 - Pad the image with background pixels so that its dimensions are a 
    %%%%        power of 2.
    %%%%    2 - Set the box size 'e' to the size of the image.
    %%%%    3 - Compute N(e), which corresponds to the number of boxes of size 'e' 
    %%%%        which contains at least one object pixel.
    %%%%    4 - If e > 1 then e = e / 2 and repeat step 3.
    %%%%    5 - Compute the points log(N(e)) x log(1/e) and use the least squares 
    %%%%        method to fit a line to the points.
    %%%%    6 - The returned Haussdorf fractal dimension D is the slope of the line.
    %%%%
    %%%%    Author
    %%%%    ------
    %%%%    Alceu Ferraz Costa 
    %%%%    email: alceufc [at] icmc [dot] usp [dot] br
    
    maxDim = max(size(I));                                                  % maximum dimension of image I
    newDimSize = 2 ^ ceil(log2(maxDim));                                    % compute the closest bigger dimension that is a power of 2
    
    rowPad = newDimSize - size(I, 1);                                       % Pad the image with background pixels so that its dimensions are a power of 2
    colPad = newDimSize - size(I, 2);
    I = padarray(I, [rowPad, colPad], 'post');                              

    boxCounts = zeros(1, ceil(log2(maxDim)));                               % initializing 
    resolutions = zeros(1, ceil(log2(maxDim)));
    
    boxSize = size(I, 1);                                                   % first dimension of image I
    boxesPerDim = 1;                                                        % initializing number of boxes by dimension to 1
    idx = 0;                                                                % initializing index of dimensions scanned
    while boxSize >= 1
        boxCount = 0;                                                       % initializing number of boxes detected in a dimension
        
        for boxRow = 1:boxesPerDim                                          % for each box
            for boxCol = 1:boxesPerDim                                  
                minRow = (boxRow - 1) * boxSize + 1;                        % compute indices of boxes
                maxRow = boxRow * boxSize;
                minCol = (boxCol - 1) * boxSize + 1;
                maxCol = boxCol * boxSize;
                
                objFound = false;                                           % initializing break variable

                for row = minRow:maxRow                                     % for each coordinates in a box
                    for col = minCol:maxCol
                        if I(row, col)                                      % if a pixel of image is found 
                            boxCount = boxCount + 1;                        % increase box count
                            objFound = true;                                % Break from nested loop.
                        end
                        
                        if objFound
                            break;                                          % Break from nested loop.
                        end
                    end
                    
                    if objFound
                        break;                                              % Break from nested loop.
                    end
                end
            end
        end
        
        idx = idx + 1;                                                      % increase index of dimension scanned
        boxCounts(idx) = boxCount;                                          % store number of boxes in which a pixel was found
        resolutions(idx) = 1 / boxSize;                                     % compute the resolutions of boxes 
        
        boxesPerDim = boxesPerDim * 2;                                      % increment number of boxes of the dimension (multiply by 2)
        boxSize = boxSize / 2;                                              % increment size of boxes in the dimension (divide by 2)
    end
    
    D = polyfit(log(resolutions), log(boxCounts), 1);                       % perform a first order polynomial fit of natural logarithms of resolutions and boxCounts
    D = D(1);                                                               % extract highest order coefficent of fitted polynomial
end
