function [processSkeleton, nbBranchpoints, nbEndpoints, ratioSkeletonProcessArea, ratioEndpointsBranchpoints] = morphoSkeleton(image, soma, processArea, resolution) 
    %%%% Computes parameters based on the detection of the skeleton of the
    %%%% processes of a cell
    %%%% 
    %%%%
    %%%% Inputs:
    %%%%    image: a binary image of a cell
    %%%%    soma: a binary image of the cell's soma
    %%%%    processArea: area of the processes of the cell in square 
    %%%%                 micrometers
    %%%%    resolution: resolution of the image, in micrometers per pixel
    %%%%        
    %%%% Outputs:
    %%%%    processSkeleton: image of the skeleton of the cell's processes
    %%%%    nbBranchpoints: total number of branching points of processes
    %%%%    nbEndpoints: total number of endpooints of the processes
    %%%%    ratioSkeletonProcessArea: ratio between the square length of
    %%%%                              the processes' skeleton and the
    %%%%                              processes total area
    %%%%    ratioEndpointsBranchpoints: ratio between total number of
    %%%%                                branching and end points

    % Skeleton of processes
    skel = bwskel(logical(image));                                          % skeletonize the image
    skelClean = bwareaopen(skel,4);                                         % remove particleq that are less than 4 pixels in area
    skeleton = soma + skelClean;                                            % add images of soma and skelClean (soma INTER slelClean pixels will be equal to 2)
    skeleton = skeleton - (skelClean&soma);                                 % substract the intersection of skelClean and soma to get a binary image back
    processSkeleton = skeleton - soma;                                      % remove the soma to get only the processes
    labelObj = bwlabel(skeleton);                                           % label regions of skeleton
    propsObj = regionprops(labelObj);                                       % compute properties of regions

    % Branchpoints
    bp = bwmorph(skelClean, 'branchpoints');                                % automatically detect branchpoints
    bp = bp - (bp&soma);                                                    % remove the union between bp and soma
    [ybp, xbp] = find(bp);                                                  % coordinates of pixels in bp

    for i = 1:length(ybp)                                                   % for each point in bp
        if bp(ybp(i) - 1, xbp(i) - 1) == 1 || ...                           % if a pixel is 8-connected to another pixel
           bp(ybp(i), xbp(i) - 1) == 1 || ...
           bp(ybp(i) + 1, xbp(i) - 1) == 1 || ...
           bp(ybp(i) - 1, xbp(i)) == 1 || ...
           bp(ybp(i) + 1, xbp(i)) == 1 || ...
           bp(ybp(i) - 1, xbp(i) + 1) == 1 || ...
           bp(ybp(i) + 1, xbp(i) + 1) == 1 || ...
           bp(ybp(i), xbp(i) + 1) == 1

           bp(ybp(i), xbp(i)) = 0;                                          % remove that pixel to get one single pixel per branchpoint
        end
    end

    % Endpoints
    ep = bwmorph(skelClean, 'endpoints');                                   % automatically detect endpoints
    ep = ep - (ep&soma);                                                    % substract the intersection between ep and the soma to get only endpoints pertaining to the processes
    
    maxBout = max([propsObj.Area]);                                         % area of biggest region in skeleton
    refineSkel = bwareaopen(skeleton, maxBout);                             % remove the smaller regions in skeleton
    flyingProc = refineSkel - soma;                                         % substract the soma to get -1 values on the intersection between refineSkel and soma
    liste = find(flyingProc == -1);                                         % coordinates of the pixels in intersection
    flyingProc(liste) = 0;                                                  % set intersection as background
    labelProc = bwlabel(flyingProc);                                        % label regions of flyingProc
    propsProc = regionprops(labelProc);                                     % compute properties of regions
    
    % Variables
    if length(propsProc) > 1                                                % if processes are segmented in several regions
        nbEndpoints = sum(ep(:)) - (length(propsProc) - 1);                 % substract false positive endpoints
    else
        nbEndpoints = sum(ep(:));
    end

    nbBranchpoints = sum(bp(:));                                            % compute number of branching points
    
    lengthSkel = sum(processSkeleton(:)) * resolution;                      % compute length of processes in micrometers
    
    if processArea > 0
        ratioSkeletonProcessArea = lengthSkel ^ 2 / processArea;                % compute ratio between skeleton length and processes area
    else
        ratioSkeletonProcessArea = 0;
    end

    if nbBranchpoints > 0                                                   % if branching points were detected
        ratioEndpointsBranchpoints = nbEndpoints / nbBranchpoints;          % compute ratio between endpoints and branchpoints
    else
        ratioEndpointsBranchpoints = 0;                                     % else, set the ratio to 0
    end
end
