function [totalIntersections, criticalRadius, dendriticMax, branchingIndex] = Sholl2pixels(image, processSkeleton, soma, Gx, Gy, resolution) 
    %%%% Computes parameters based on the Sholl Analysis with circles of
    %%%% radii of 2 pixels increment, centered on the cell's soma centroid
    %%%%
    %%%% Inputs:
    %%%%    image: a binary image of a cell
    %%%%    processSkeleton: a binary image of the skeleton of the cell's
    %%%%                     processes
    %%%%    soma: a binary image of the cell's soma
    %%%%    Gx: x coordinate of the cell's soma centroid
    %%%%    Gy: y coordinate of the cell's soma centroid
    %%%%    resolution: resolution of the images, in micrometers per pixel
    %%%%        
    %%%% Outputs:
    %%%%    totalIntersections: total number of intersections between 
    %%%%                        processes and all Sholl circles
    %%%%    criticalRadius: index of the first Sholl circle for which the 
    %%%%                    number of intersections is maximized
    %%%%    dendriticMax: maximum number of intersections between the 
    %%%%                  cell's processes and a Sholl circle
    %%%%    branchingIndex: Sholl analysis branching index, describing 
    %%%%                    the ramification of a cell's processes
    

    [ybd, xbd] = find(bwmorph(soma, 'remove'));                             % indices of the points of the circumference of the soma
    
    for m = 1:length(ybd)                                                   % distance of each point of the circumference of the soma to its centroid
        distCentroid(m, :) = norm([ybd(m); xbd(m)] - [Gy; Gx]);
    end

    distMaxCentroid = max(distCentroid);                                    % maximum distance to the centroid
    radius = round(distMaxCentroid) + 1;                                    % radius of the circumscribing circle to the soma

    [ysq, xsq] = find(processSkeleton);                                     % indices of the points of the skeleton of the processes
    
    if isempty(ysq) == 0                                                    % if the cell has processes
        
        for m = 1:length(ysq)                                               % distance of each point of the skeleton of the processes to the centroid of the soma
            distSkel(m,:) = round(norm([ysq(m); xsq(m)] - [Gy; Gx]));
        end

        [cy, cx] = size(image);                                             % indices of all the points of the cell
        N2 = max([Gx, Gy, (cx - Gx), (cy - Gy)]);                           % maximum distance in pxels between centroid and the borders of the image
        nbRadius = round((N2 - radius) / 2) - 2;                            % maximum number of Sholl circles that fit between the circumscribing circle and the borders of the image           
        
        for t = 1:nbRadius                                                  % for each successive Sholl circle, starting with the circumscribing circle
            nbIntersections(:, t) = length(find(distSkel == radius));       % number of intersections between the skeleton of the processes and each Sholl circle
            radius = radius + 2;                                            % increase radius by 2 pixel (next Sholl circle)
        end

        totalIntersections = sum(nbIntersections);                          % total number of intersections between processes and all Sholl circles
        BI = 0;

        for t = 2:nbRadius                                                  
            diffIntersections = nbIntersections(t) - nbIntersections(t-1);  % comparing number of intersections for successive Shol circles
            if diffIntersections > 0                                        % if number of intersections increases
                BI = BI + diffIntersections * (t * resolution)...           % add to the branching index 
                                            / sqrt(cy * cx);
            end
        end
        
        maxIntersections = find(nbIntersections == max(nbIntersections));   
        criticalRadius = maxIntersections(1, 1);                            % index of the first Sholl circle for which the number of intersections is maximized
        dendriticMax = max(nbIntersections);                                % maximum number of intersections with a Sholl circle
        branchingIndex = BI;  
    
    else                                                                    % if the cell has no processes, set parameters to 0
        totalIntersections = 0;
        dendriticMax = 0;
        criticalRadius = 0;
        branchingIndex = 0;
    end
end
