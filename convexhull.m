function [convexHullArea,convexHullPerim,convexity, solidity, convHullCircularity, convHullRadii, convHullSpanRatio] = convexhull(image, cellPerimeter, resolution) 
    %%%% Computes parameters based on the convex hull of a cell. The convex
    %%%% hull is the smallest convex polygon that fully encloses the object
    %%%%
    %%%% Inputs:
    %%%%    image: a binary image of a cell
    %%%%    resolution: resolution of the images, in micrometers per pixel
    %%%%    cellPerimeter: perimeter of the image
    %%%%
    %%%% Outputs:
    %%%%    convexity: ratio between the perimeters of the convex hull of
    %%%%               the cell and that of the cell
    %%%%    solidity: ratio between the areas of the cell and that of its
    %%%%              convex hull
    %%%%    convHullCircularity: circularity of the convex hull
    %%%%    convHullRadii: ratio of the longest on the shortest distances
    %%%%                   of points of the convex hull to centroid of the
    %%%%                   convex hull
    %%%%    convHullSpanRatio: ratio between major and minor axis lengths
    %%%%                       of convex hull. The axes are line segments 
    %%%%                       that can be drawn within the region while 
    %%%%                       still being entirely contained within it.


    convexHull = bwconvhull(image);                                         % image of the convex hull of the cell in image
    convexHullArea = sum(convexHull(:));                                    % area of the convex hull
    convexHullPerim = bwperim(convexHull);                                  % perimeter of the convex hull
    solidity = sum(image(:)) / convexHullArea;                              % compute solidity
    convexity = cellPerimeter / (sum(convexHullPerim(:)) * resolution);     % compute convexity

    convHullCircularity = 2 * sqrt(pi * sum(convexHull(:) * ...             % compute circularity of the convex hull
                                   resolution ^ 2)) / ...
                          (sum(convexHullPerim(:)) * resolution);

    centCH = regionprops(convexHull, 'Centroid');
    [yCvhp, xCvhp] = find(convexHullPerim);

    for z = 1:length(yCvhp)                                                 % for each point of the convex hull of the image
        distCH(z, :) = norm([yCvhp(z); xCvhp(z)] - ...                      % distance between the point and the centroid of the convex hull
                            [centCH.Centroid(2); centCH.Centroid(1)]);
    end

    largestRadiusCH = max(distCH) * resolution;                             % longest distance between a point of the convex hull and its centroid, in pixels
    smallestRadiusCH = min(distCH) * resolution;                            % shortest distance between a point of the convex hull and its centroid, in pixels

    lengthAxisCvh = regionprops(convexHull, 'MajorAxisLength', ...          % extract major and minor axes of the convex hull
                                'MinorAxisLength');

    convHullRadii = largestRadiusCH / smallestRadiusCH;                     % compute the ratio of optima radii of the convex hull

    convHullSpanRatio = lengthAxisCvh.MajorAxisLength / ...                 % compute the convex hull span ratio
                        lengthAxisCvh.MinorAxisLength; 
                    
    convexHullArea = convexHullArea*resolution*resolution;                  % in microns
    convexHullPerim = convexHullPerim *resolution;

end
