clear variables
clc
close all

%% Get binary image files
extension = '*.tif*';
path = uigetdir('', ['Select the folder containing the binarized ' ...
                'individualized 300*300pi microglia']);

[files, nFiles, fileNo, fileNames, images] ...
        = Open_Microglia_Images(path, extension);

filePattern = fullfile(path, extension);
theFiles = dir(filePattern);

%% Dectection of bad segmentation images
id_wrong = [];
for n = 1:nFiles
    X = images(n).R(:, :, 1);
    borders = sum(X(:,1)) + sum(X(:,end)) ...
            + sum(X(1,:)) + sum(X(end,:));
    total = sum(X(:));
    if borders ~= 0 || total < 256
        id_wrong = [id_wrong, n];
    end
end

if  ~isempty(id_wrong)
    disp("The following images are not correct !")
    disp(fileNames(id_wrong)')
    return
else
    disp("All images are correct !")
end

%% Image resolution
answerResolution = inputdlg({['Please enter the value of 1 pixel ' ...
                            '(in micrometers):']}, 'Resolution', [1 35]);
resolution = str2num(answerResolution{1});

%% Save Path
savePath = uigetdir('','Select the folder for saving the results tables');

%% Check Cell body & Skeleton
answerCheck = questdlg(['Do you want to check the detection of the ' ...
                        'cell bodies and skeleton ? it will make the ' ...
                        'analysis longer'], 'Parameters', 'Yes', 'No', ...
                        'Yes');

if strcmpi('Yes', answerCheck)
    saving_image_path = uigetdir('', ['Select the folder where you want' ...
                                 ' to save the analyzed images']);
    % Selectionne aléatoirement 30 cellules parmi le jeu de données
    % Affiche les microglies source, cell body par dessus microglie source,
    % squelette aussi
end

%% Pre-allocating the variables
somaArea = zeros(1, nFiles);
perimAreaRatio = zeros(1, nFiles);
circularity = zeros(1, nFiles);
density = zeros(1, nFiles);
roundFactor = zeros(1, nFiles);
ramificationIndex = zeros(1, nFiles);
convexity = zeros(1, nFiles);
solidity = zeros(1, nFiles);
convHullCircularity = zeros(1, nFiles);
convHullRadii = zeros(1, nFiles);
convHullSpanRatio = zeros(1, nFiles);
convexHullArea = zeros(1, nFiles);
convexHullPerim = zeros(1, nFiles);
cellArea = zeros(1, nFiles);
cellPerimeter = zeros(1, nFiles);
processesArea = zeros(1, nFiles);
ratioProcessSomaArea = zeros(1, nFiles);
ratioProcessCellArea = zeros(1, nFiles);
fractalDimension = zeros(1, nFiles);
nbBranchpoints = zeros(1, nFiles);
nbEndpoints = zeros(1, nFiles);
ratioSkeletonProcessArea = zeros(1, nFiles);
ratioEndpointsBranchpoints = zeros(1, nFiles);
totalIntersections = zeros(1, nFiles);
criticalRadius = zeros(1, nFiles);
dendriticMax = zeros(1, nFiles);
branchingIndex = zeros(1, nFiles);
polarizationIndex = zeros(1, nFiles);
linearity = zeros(1, nFiles);
inertiaCell = zeros(1, nFiles);
lacunaritySlope = zeros(1, nFiles);
lacunarityMean = zeros(1, nFiles); 

%% Computation of the morphological parameters
for n = 1:nFiles
    disp(n)
    X = images(n).R(:, :, 1);
    X = double(X);
    X = X / max(X(:));                                                      % binarized microglia image
    
    [somaArea(n), Gx, Gy, soma] = Cell_Body(X, resolution);
     cellArea(n)=sum(X(:))*resolution*resolution;                           % Surface cellule
    processesArea(n) = cellArea(n) - somaArea(n);
    
    [perimAreaRatio(n), cellPerimeter(n)] ... 
        = fftPerimAreaRatio(X, resolution);                                 % log(Perimeter²/Area) 

    circularity(n) = 2 * sqrt(pi * cellArea(n)) / cellPerimeter(n);
    density(n) = cellArea(n) / (size(X,1) * size(X,2) * resolution ^ 2);                  
    roundFactor(n) = roundnessFactor(X, somaArea(n), resolution);

    ramificationIndex(n) ...
        = (cellPerimeter(n) / cellArea(n)) / (2 * sqrt(pi / cellArea(n)));

    [convexHullArea(n), convexHullPerim(n), convexity(n), solidity(n), ...
     convHullCircularity(n), convHullRadii(n), convHullSpanRatio(n)] ...
        = convexhull(X, cellPerimeter(n), resolution);

    processesArea(n) = cellArea(n) - somaArea(n);
    ratioProcessSomaArea(n) = processesArea(n) / somaArea(n);
    ratioProcessCellArea(n) = processesArea(n) / cellArea(n);

    [processSkeleton, nbBranchpoints(n), nbEndpoints(n), ...
        ratioSkeletonProcessArea(n), ratioEndpointsBranchpoints(n)] ...
        = morphoSkeleton(X, soma, processesArea(n), resolution);

    [totalIntersections(n), criticalRadius(n), dendriticMax(n), ...
        branchingIndex(n)] ...
        = Sholl2pixels(X, processSkeleton, soma, Gx, Gy, resolution);

    fractalDimension(n) = hausDim(X);
    [lacunaritySlope(n), lacunarityMean(n)] = lacunarity_glbox(X);

    [polarizationIndex(n), linearity(n)] ...
        = morphoPolarizationLinearity(X, Gx, Gy, resolution);

    inertiaCell(n) = inertia(X);

    % saving an image of the analysis
    if strcmpi('Yes', answerCheck)
        colors = [1      0      0;
                  0.5      1      0.5;
                  0      1      0;
                  1      0.5    0;
                  0.4    0.7608 0.6470; 
                  0.9882 0.5529 0.3843; 
                  0.5529 0.6275 0.7961; 
                  0.9059 0.5412 0.7647; 
                  0.651  0.8471 0.3294; 
                  1      0.851  0.1843; 
                  0.8980 0.7686 0.5804; 
                  1      1      1];


        c1 = 12;
        c2 = 3;
        c3 = 1;

        idx1 = find(X - soma - processSkeleton);
        redImg = 255 * zeros(size(X,1), size(X,2), 'uint8');
        greenImg = 255 * zeros(size(X,1), size(X,2), 'uint8');
        blueImg = 255 * zeros(size(X,1), size(X,2), 'uint8');                    
        redImg(idx1) = colors(c1, 1) * 255;
        greenImg(idx1) = colors(c1, 2) * 255;
        blueImg(idx1) = colors(c1, 3) * 255;
        im1 = cat(3, redImg, greenImg, blueImg);
        
        idx2 = find(soma);
        redImg = 255 * zeros(size(X,1), size(X,2), 'uint8');
        greenImg = 255 * zeros(size(X,1), size(X,2), 'uint8');
        blueImg = 255 * zeros(size(X,1), size(X,2), 'uint8');                    
        redImg(idx2) = colors(c2, 1) * 255;
        greenImg(idx2) = colors(c2, 2) * 255;
        blueImg(idx2) = colors(c2, 3) * 255;
        im2 = cat(3, redImg, greenImg, blueImg);
        
        idx3 = find(processSkeleton);
        redImg = 255 * zeros(size(X,1), size(X,2), 'uint8');
        greenImg = 255 * zeros(size(X,1), size(X,2), 'uint8');
        blueImg = 255 * zeros(size(X,1), size(X,2), 'uint8');                   
        redImg(idx3) = colors(c3, 1) * 255;
        greenImg(idx3) = colors(c3, 2) * 255;
        blueImg(idx3) = colors(c3, 3) * 255;
        im3 = cat(3, redImg, greenImg, blueImg);

        im = im1 + im2 + im3;

        f = figure();
        imshow(im)
        % FigName = ['cell' num2str(n)];
        FigName = fileNames{n};
        saveas(f, fullfile(saving_image_path, [FigName '.png']));
        close(f)
    end
end

%% Save data in xls file
% Set parameter names

header = {'microglia', 'CellArea', 'CellPerimeter', 'ConvexHullPerimeter', ...
        'ConvexHullArea', 'ProcessesArea', 'SomaArea', 'NbBranchpointsSkeleton', ...
        'NbEndpointsSkeleton', 'ShollTotalIntersections', 'ShollcriticalRadius', ...
        'ShollDendriticMax', 'PerimeterAreaRatio', 'Circularity', ...
        'RoundnessFactor', 'RamificationIndex', 'Solidity', 'Convexity', ...
        'ConvexHullRadiiRatio', 'SpanRatio', 'ConvexHullCircularity', ...
        'ProcessesSomaAreasRatio', 'ProcessesCellAreasRatio', 'Density', ...
        'BranchingIndex', 'EndpointsBranchpointsRatio', ...
        'SkeletonProcessesRatio', 'FractalDimension', 'LacunaritySlope', ...
        'PolarizationIndex', 'Linearity', 'Inertia'};

% Set table values
allParameters = table(fileNames', cellArea', cellPerimeter', ...
                convexHullPerim', convexHullArea', processesArea', ...
                somaArea', nbBranchpoints', nbEndpoints', totalIntersections', ...
                criticalRadius', dendriticMax', perimAreaRatio', circularity', ...
                roundFactor', ramificationIndex', solidity', convexity', ...
                convHullRadii', convHullSpanRatio', convHullCircularity', ...
                ratioProcessSomaArea', ratioProcessCellArea', density', ...
                branchingIndex', ratioEndpointsBranchpoints', ...
                ratioSkeletonProcessArea', fractalDimension', ...
                lacunaritySlope', polarizationIndex', linearity', ...
                inertiaCell', 'VariableNames', header);

% Delete old output file if exists and save file
if exist([savePath './Parameters_MorphoCellSorter.xlsx'], 'file') == 2
    delete([savePath './Parameters_MorphoCellSorter.xlsx']);
end

writetable(allParameters, [savePath './Parameters_MorphoCellSorter.xlsx'], ...
           'WriteRowNames', true);


%% TABLE MCS
allParameters_forMCS = table(perimAreaRatio', circularity', ...
                roundFactor', ramificationIndex', solidity', convexity', ...
                convHullRadii', convHullSpanRatio', convHullCircularity', ...
                ratioProcessSomaArea', ratioProcessCellArea', density', ...
                branchingIndex', ratioEndpointsBranchpoints', ...
                ratioSkeletonProcessArea', fractalDimension', ...
                lacunaritySlope', polarizationIndex', linearity', ...
                inertiaCell');
allParameters_forMCS=table2array(allParameters_forMCS);


for r=1:nFiles
    allParameters_forMCS(r,2)=1/ allParameters_forMCS(r,2);
    allParameters_forMCS(r,3)=1/ allParameters_forMCS(r,3);
    allParameters_forMCS(r,5)=1/ allParameters_forMCS(r,5);
    allParameters_forMCS(r,7)=1/ allParameters_forMCS(r,3);
    allParameters_forMCS(r,8)=1/ allParameters_forMCS(r,3);
end
header = {'microglia', 'PerimeterAreaRatio', 'Circularity', ...
        'RoundnessFactor', 'RamificationIndex', 'Solidity', 'Convexity', ...
        'ConvexHullRadiiRatio', 'SpanRatio', 'ConvexHullCircularity', ...
        'ProcessesSomaAreasRatio', 'ProcessesCellAreasRatio', 'Density', ...
        'BranchingIndex', 'EndpointsBranchpointsRatio', ...
        'SkeletonProcessesRatio', 'FractalDimension', 'LacunaritySlope', ...
        'PolarizationIndex', 'Linearity', 'Inertia'};
table_forMCS=table(fileNames',allParameters_forMCS(:,1),allParameters_forMCS(:,2),allParameters_forMCS(:,3),...
    allParameters_forMCS(:,4),allParameters_forMCS(:,5),allParameters_forMCS(:,6),allParameters_forMCS(:,7),...
    allParameters_forMCS(:,8),allParameters_forMCS(:,9),allParameters_forMCS(:,10),allParameters_forMCS(:,11),...
    allParameters_forMCS(:,12),allParameters_forMCS(:,13),allParameters_forMCS(:,14),allParameters_forMCS(:,15),...
    allParameters_forMCS(:,16),allParameters_forMCS(:,17),...
    allParameters_forMCS(:,18),allParameters_forMCS(:,19),allParameters_forMCS(:,20),'VariableNames', header);



if exist([savePath './table_MCS.xlsx'], 'file') == 2
    delete([savePath './table_MCS.xlsx']);
end

writetable(table_forMCS, [savePath './table_MCS.xlsx'], ...
           'WriteRowNames', true);
