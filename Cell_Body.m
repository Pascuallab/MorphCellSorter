function [somaArea, Gx, Gy, soma] = Cell_Body(image, resolution)
    %%%% Cell_Body locates the cell body (soma) of a microglia (or cell
    %%%% having a similar structure)
    %%%%
    %%%% Inputs:
    %%%%    image: a segmented and binarized square image of the cell
    %%%%    resolution: number of micrometers in one pixel
    %%%%
    %%%% Outputs:
    %%%%    somaArea: area of the soma, in square micrometers
    %%%%    Gx: x coordinate of the center of the cell body
    %%%%    Gy: y coordinate of the center of the cell body
    %%%%    soma: binary image of the cell body
    
    %% Compute the cell's center
    [rows, cols] = size(image);                                             % dimensions of image
    Z = image;                                                              % copy of image
    h = ones(3, 3);                                                         % "box blur" kernel
    s = Inf;                                                                % initializing s and k
    k = 1;

    while s > 0                                                             % until surface of Z is 0 
        Zold = Z;                                                           % saving Z in Zold
        Z = imerode(Z, h);                                                  % eroding Z by kernel h
        s = sum(Z(:));                                                      
        Surface(k) = s;                                                     % saving the surface of Z in Surface
        k = k + 1;                                                          % incrementing index
    end

    Z = Zold;                                                               % restore Z to one step before its surface is 0

    if max(Z(:)) ~= 1                                                       % if Z is not binary
        indexv = find(Z > 0);                                               
        Z(indexv) = 1;                                                      % binarize Z
    end

    [uy, ux] = find(Z == 1);                                                % indices of pixels in Z
    xo = round(mean(ux));                                                   % coordinates of gravity center of Z, the eroded image
    yo = round(mean(uy));
    
    %% Compute the Inertia Momentum
    [uy, ux] = find(image == 1);                                            % indices of pixels in image
    I = ((xo - ux) .^ 2 + (uy - yo) .^ 2);                                  % compute the inertia momentum of image based on gravity center of Z
    I = sum(I(:));
    
    %% Morphological Growing
    h = [0 1 0; 1 1 1; 0 1 0];                                              % "plus" kernel
    W = zeros(rows, cols);                                                  % initialize growing image W...
    W(yo, xo) = 1;                                                          % ... and seed it at the gravity center of the eroded Z
    sto = 0;                                                                % initialize step surface values sto and stn
    stn = 10000;
    Surf(1) = 1;                                                            % initialize morphological growing surface
    k = 1;                                                                  % initialize step index

    while stn ~= sto || sum(sum(W .* image)) ~= sum(image(:))               % as long as surface changes from step to step
        k = k + 1;                                                          % increment step index
        sto = stn;                                                          % save previous surface value
        W = imdilate(W, h);                                                 % dilate W with kernel h
        Wt = W .* image;                                                    
        stn = sum(Wt(:));                                                   % compute the surface of the intersection of image and growing image W
        Surf(k) = stn;                                                      % store value in Surf
    end
    

    %% Build back the cell body
    NN = find(Surface < Surf(1:length(Surface)), 1, 'first');               % compare erosion and dilation surfaces to find lowest index at which all the soma is isolated
    h1 = [1 0 1; 0 1 0; 1 0 1];                                             % "cross" kernel
    h2 = [0 1 0; 1 1 1; 0 1 0];                                             % "plus" kernel
    soma1 = image;                                                          % copies of image for "cross" and "plus" eroding and dilating
    soma2 = image;

    if I > 50000                                                            % if inertia momentum is high
        for n = 1:NN+2                                                      % erode NN+2 times with kernel h1, or until full erosion
            if sum(sum(imerode(soma1, h1))) > 0
                soma1 = imerode(soma1, h1);                                 
            end

            if sum(sum(imerode(soma2, h2))) > 0                             % erode NN+2 times with kernel h2, or until full erosion
                soma2 = imerode(soma2, h2);
            end
        end

        for n = 1:NN+1                                                      
            soma1 = imdilate(soma1, h1);                                    % regrow soma by dilating with kernel h1 and intersecting with image, NN+1 times
            soma1 = soma1 .* image;
            soma2 = imdilate(soma2, h2);                                    % regrow soma by dilating with kernel h2 and intersecting with image, NN+1 times
            soma2 = soma2 .* image;
        end

    else                                                                    % if inertia momentum is not high, perform like above, one less time
        for n = 1:NN+1
            if sum(sum(imerode(soma1, h1))) > 0
                soma1 = imerode(soma1, h1);
            end

            if sum(sum(imerode(soma2, h2))) > 0
                soma2 = imerode(soma2, h2);
            end
        end

        for n = 1:NN
            soma1 = imdilate(soma1, h1);
            soma1 = soma1 .* image;
            soma2 = imdilate(soma2, h2);
            soma2 = soma2 .* image;
        end
    end

    soma = 0.5 * (soma1 + soma2) > 0.25;                                    % union between soma1 and soma2 and binarization                                             
    [uy, ux] = find(soma > 0);                                              % indices of pixels in soma
    
    xop = round(mean(ux));                                                  % coordinates of the gravity center of the soma
    yop = round(mean(uy));
    
    %% Postprocessing 
    %% Remove possible checkboards-like patterns resulting from dilations
    [uy, ux] = find(soma);                                                  % indices of pixels in soma

    for n = 1:length(uy)
        for m = 1:length(ux)                                                % for each pixel in soma
            if soma(uy(n), ux(m)) ~= 0                                      % if the pixel in soma is not 0
                if soma(uy(n) + 1, ux(m)) == 0 ...                          % if the pixel is orthogonally surrounded by 0
                    && soma(uy(n) - 1, ux(m)) == 0 ...
                    && soma(uy(n), ux(m) + 1) == 0 ...
                    && soma(uy(n), ux(m) - 1) == 0

                    soma(uy(n), ux(m)) = 0;                                 % remove the pixel

                elseif soma(uy(n) + 1, ux(m)) == 0 ...                      % if the pixel is orthogonally surrounded by 0 except by one side
                    && soma(uy(n) - 1, ux(m)) == 0 ...
                    && soma(uy(n), ux(m) + 1) == 0 ...
                    || soma(uy(n) - 1, ux(m)) == 0 ...
                    && soma(uy(n), ux(m) + 1) == 0 ...
                    && soma(uy(n), ux(m) - 1) == 0 ...
                    || soma(uy(n) + 1, ux(m)) == 0 ...
                    && soma(uy(n) - 1, ux(m)) == 0 ...
                    && soma(uy(n), ux(m) - 1) == 0 ...
                    || soma(uy(n) + 1, ux(m)) == 0 ...
                    && soma(uy(n), ux(m) + 1) == 0 ...
                    && soma(uy(n), ux(m) - 1) == 0

                    soma(uy(n), ux(m)) = 0;                                 % remove the pixel
                end
            end
        end
    end

    %% If more than one soma is detected, select the biggest one
    props = regionprops(logical(soma));                                              % extract properties of different regions in soma

    if length(props) > 1                                                    % if more than one region is detected
        best = 0;                                                           % initialize
        score = - Inf;

        for reg = 1:length(props)                                           % for each region
            regScore = props(reg).Area;                                     % exctract area of region

            if regScore > score                                             % detect highest area and save index of region
                best = reg;
                score = regScore; 
            end
        end

        for reg = 1:length(props)                                           % for each region
            if reg ~= best                                                  % if region is not biggest, remove region
                soma(fix(props(reg).BoundingBox(2)):fix(props(reg).BoundingBox(2) + props(reg).BoundingBox(4)), ...
                     fix(props(reg).BoundingBox(1)):fix(props(reg).BoundingBox(1) + props(reg).BoundingBox(3))) = 0;
            end
        end

        xop = round(props(best).Centroid(1));                               % compute coordinates of new gravity center of soma
        yop = round(props(best).Centroid(2));
    end
    
    %% another round of removing potential additional regions in soma
    labCC = bwlabel(soma);                                                  
    propsCC = regionprops(labCC, 'Area'); 
    areasoma = cat(1, propsCC.Area);

    if length(propsCC) > 1
        surfCCmax = max(areasoma);
        Xpp2 = bwareaopen(soma, surfCCmax);
        soma = image&Xpp2;
    end
    
    %% based on the size of the soma, erode then dilate with a disk kernel
    somaArea = sum(soma(:));                                                % compute area of soma 
        
    if somaArea < 100                                                       % erode soma in a disk shape, 1, 3 or 5 times, depending on its area
        XppEr = imerode(soma, strel('disk', 1));
        indice = 1;

    elseif somaArea > 300
        XppEr = imerode(soma, strel('disk', 5));
        indice = 2;

    else
        XppEr = imerode(soma, strel('disk', 3));
        indice = 3;
    end

    labCC2 = bwlabel(XppEr);
    propsCC2 = regionprops(labCC2, 'Area'); 
    areasoma2 = cat(1, propsCC2.Area);

    if length(propsCC2) > 1
        surfCCmax2 = max(areasoma2);
        Xpp22 = bwareaopen(XppEr, surfCCmax2);                              % remove additional regions

        if indice == 1                                                      % dilate soma the same number of times as it was eroded
            Xpp2n = imdilate(Xpp22, strel('disk', 1));

        elseif indice == 2
            Xpp2n = imdilate(Xpp22, strel('disk', 5));

        else
            Xpp2n = imdilate(Xpp22, strel('disk', 3));

        end

        soma = image&Xpp2n;                                                 % intersection between image and erosion dilation image
    end
    
    %% compute parameters
    clear uy ux
    somaArea = sum(soma(:)) * resolution ^ 2 ;                              % area of soma in square micrometers
    [uy, ux] = find(soma > 0);                                              % coordinates of pixels in soma
    Gx = round(mean(ux));                                                   % coordinates of gravity center of final soma
    Gy = round(mean(uy));

%     if indice ==2
%         figure()
%         imagesc(image+2*soma)
%         axis off
%         axis equal
%     end
end
