clear variables
clc
close all

%% Import data table
extension = '*.tif*';
[xlsFile, xlsPath] = uigetfile('*.xlsx', 'Select a table');                 % Select table containing the morphological parameters
allParameters = readtable(fullfile(xlsPath, xlsFile));  % Importthe table

allParameters = sortrows(allParameters,1,"ascend"); % sort rows in chrono

%  saving_path='E:\Programmes Matlab\figures_article0424\table1\Embryo\0.9';
X = table2array(allParameters(:, 2:end));  % Extract all parameters except image names

fileNames = allParameters{:,1};
%fileNames = allParameters{(X(:,1) + 6.66*X(:,2)-177)>0, 1};
%X = X((X(:,1) + 6.66*X(:,2)-177)>0, 3:end); % to remove CellArea < 100 and CellArea and SomaAre columns


Label_variable = {'Peri2/Area',	'Circ'	'Roundness', 'RamifIndex',...	    
    'Solidity',	'Convexity', 'CHRadiiRatio', 'SpanRatio',...	
    'CHCirc', 'Processes/Soma',...
    'Processes/Cell', 'Density', 'BranchingIndex',...
    'Endp/Branchp', 'Skeleton/Processes',...
    'FractalDim', 'Lacunarity', 'PolarizationIndex',...	
    'Linearity', 'Inertia'};

%% PCA
[rows, cols] = size(X);
[Xp, L, W, Wp, Phi, R] = PCA_custom(X);

figure();
bar(L ./ sum(L));
hold on;
plot(cumsum(L ./ sum(L)));
title('Eigenvalues (Inertia CP_i)');
xlabel('CP_i');
ylabel('Eigenvalues (Normalized)');

%------------------------------------------------------------------------

figure();
quiver(0 * W(:,1), 0 * W(:,2), W(:,1), W(:,2));
hold on;
t = 0 : 0.01 : 6.28;
t = t';
x = cos(t);
y = sin(t);
plot(x, y);
plot(0.707 * x, 0.707 * y);

plot(Xp(:, 1) ./ max(Xp(:, 1)), ...
    (L(2) / L(1)) * Xp(:, 2) ./ max(Xp(:, 2)), 'o');

text(W(:, 1), W(:, 2), string(Label_variable));
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
axis equal
grid minor
title('Correlation Circle');
xlabel('PC1');
ylabel('PC2');



close;



%---------------------------------------------------------------------------

figure();
Wp(:, 2) = abs(Wp(:, 2)) ./ max(abs(Wp(:, 2)));                             % Normalize to PC1 amplitudes
wp = sortrows(Wp, 2, 'descend');                                            % Sort variables based on their projection onto PC1

bar(1:length(wp), wp(:, 2));
hold on;

for k = 1:20
    txt = text(k, 0, Label_variable(wp(k, 1)));
    Label_variable(wp(k, 1));
    txt.Rotation = 90;
    
end
grid on
yyaxis right;
wpn = cumsum(wp(:, 2) ./ sum(wp(:, 2)));

answerThresh = inputdlg(['Choose a threshold for how many parameters ' ...
                         'to include (values between 0.1 and 1)'], ...
                         'Threshold', [1 43], {'0.8'});
if str2num(answerThresh{1}) == 1
    var_max=20;
else
    var_max = find(wpn > str2num(answerThresh{1}), 1, 'first') - 1;
end

plot(1:20, wpn, 'm');

%----------------------------------------------------------------------
 
M1 = X(:, wp(1:var_max, 1));
M1 = (M1 - mean(M1)) ./ std(M1);
M1 = M1 .* (wp(1:var_max, 2)' ./ max(wp(1:var_max, 2)));

Phi_m = Phi(wp(1:var_max));
cos_sin = zeros(length(t), var_max);

for n = 1:var_max
    cos_sin(:, n) = cos(n * t + Phi_m(n));
end

%% Andrews plot
clear S1
S1 = M1 * cos_sin';

figure();
t = (t * 180 / pi);
plot(t, S1');

t_max = find(var(S1) == max(var(S1)));
teta = t(t_max);
hold on;
plot(teta * [1 1], [-1 1] * cols, 'b', 'linewidth', 2);

score = [(1:rows)' S1(:, t_max)];
rank = sortrows(score, 2, 'ascend');
firstMg=find(rank(:,1)==1);
lastMg=find(rank(:,1)==size(allParameters,1));


title('Andrews Plot');
xlabel('Virtual angle');
ylabel('Andrews Score');

text(teta + 5, 3 * cols / 4, ...
    ['var_{max} \rightarrow \theta=' num2str(round(teta)) '^o']);


close;



andrewsValuesMaxVar = S1(:, t_max);
table_rank = table(rank);
vars = {'rank', 'score'};
% table_rank=table_rank(1:end,vars);
% fig = figure();
% uit = uitable(fig, 'Data', rank, 'columnName', vars, 'Units', ...
%               'normalized', 'position', [0 0 1 1]);
% uit.FontSize = 10;

ranking = rank(:, 1);

fileNames = sortrows(fileNames,1,'ascend');
microgliaIdsSorted = fileNames(ranking)';

%% Results of the ranking
answerFigures = menu(['Do you want to display the rankings and/or ' ...
                      'distribution histograms of Andrews Scores?'], ...
                      'Rankings', 'Histograms', 'Both', 'None');

% Number and name of conditions
answerCondition = inputdlg({['Please enter the number of ' ...
                             'conditions:']}, ...
                             'Number of conditions', [1 35]);

nbConditions = str2num(answerCondition{1});
clear answerCondition

for t = 1:nbConditions
    nameCond{1, t} = ['Name Condition', num2str(t)];
end

nameConditions = inputdlg(nameCond, 'Conditions name', [1 35]);
classes = zeros(1, rows);

for t = 1:nbConditions
    findPattern = contains(fileNames, ...
                           nameConditions{t, 1});

    findidx = find(findPattern);
    classes(findidx) = t;
    % clear findPattern findidx
end

classesIdsSorted = classes(ranking);

if answerFigures < 4                                                        % if answer is not 'None'
    answerDisplay = questdlg(['Do you want to display the results ' ...
                              'color-coded according to the ' ...
                              'experimental conditions?'], 'Parameters', ...
                              'Yes', 'No', 'Yes');

    if strcmpi('Yes', answerDisplay)
        % classesIdsSorted = zeros(rows, 1);
        % 
        % for i = 1:rows
        %     numero(i, :) = find(allParameters{:, 1} ...
        %                         == microgliaIdsSorted(:, i));
        % 
        %     classesIdsSorted(i, 1) = classes(:, numero(i, :));
        % end

        % Selection of the display colors
        defaultColors = [0.4    0.7608 0.6470; 
                         0.9882 0.5529 0.3843; 
                         0.5529 0.6275 0.7961; 
                         0.9059 0.5412 0.7647; 
                         0.651  0.8471 0.3294; 
                         1      0.851  0.1843; 
                         0.8980 0.7686 0.5804; 
                         1      1      1];

        if nbConditions <= 8
            answerColors = questdlg(['Do you want to choose the display ' ...
                                     'colors for your conditions?'], ...
                                     'Color display', 'Yes', 'No', 'No');
            
            if strcmpi('Yes', answerColors)
                for t = 1:nbConditions
                    colors(t, :) = uisetcolor(defaultColors(t, :), ...
                                   ['Select color for Condition', ...
                                   num2str(t), ' ', nameConditions{t, 1}]);
                end

            elseif strcmpi('No', answerColors)
                colors = defaultColors;
            end

        else
            for t = 1:nbConditions
                colors(t, :) = uisetcolor([1 1 1], ...
                               ['Select color for Condition', ...
                               num2str(t), ' ', nameConditions{t, 1}]);
            end
        end

        % Display ranking
        if answerFigures == 1 || answerFigures == 3
            path = uigetdir('', ['Select the folder containing ' ...
                            'the binarized individualized ' ...
                            '300*300pi microglia']);
            % [~, ~, ~, ...
            %     fileNames, ~] = Open_Microglia_Images(path, extension);

            %%filePattern = fullfile(path, extension);
            %%theFiles = dir(filePattern);
            % fileNames = sortrows(fileNames,1,'ascend');
            % microgliaIdsSorted = fileNames(ranking);                        % Get microglia ranking with sort
            DisplayImages(path, extension, colors, microgliaIdsSorted', ...
                          classesIdsSorted)
        end

        % Display histograms
        if answerFigures == 2 || answerFigures == 3
            for t = 1:nbConditions
                 groupe(t).R = find(classes(1, :) == t);
    
                 for i = 1:length(groupe(t).R)
                     groupList(t).R(i, 1) ...
                         = andrewsValuesMaxVar((groupe(t).R(1,i)),1);
                 end
            end

            minCombined=zeros(1, nbConditions);
            maxCombined=zeros(1, nbConditions);

            for t = 1:nbConditions
                [minCombined(t), maxCombined(t)] = bounds(groupList(t).R);
            end

            minC = min(minCombined(:));
            maxC = max(maxCombined(:));

            clear minCombined maxCombined
            %nb bins
            nx = numel(andrewsValuesMaxVar);

            % Calculate the 25th and 75th percentiles
            p25 = prctile(andrewsValuesMaxVar, 25);
            p75 = prctile(andrewsValuesMaxVar, 75);

            % Calculate the optimal bin width
            width = 2 * (p75 - p25) / nthroot(nx, 3);

            % Calculate the number of bins
            nbins = ceil((max(andrewsValuesMaxVar) ...
                    - min(andrewsValuesMaxVar)) / width);

            nbins = max(1, nbins);

            BE = linspace(minC, maxC, nbins);

            for t = 1:nbConditions
                bincounts(t).R = histcounts(groupList(t).R, BE, ...
                                 'Normalization','probability');
            end


            % Distributions display
            nbRows=round(nbConditions/2);
            nbCols=round(nbConditions/nbRows);

            figure()
            % set(gcf, 'Position', get(0, 'Screensize'));

            % Compute the max of all distributions
            maxdistrib = 0; 
            for t = 1:nbConditions
                if max(bincounts(t).R) > maxdistrib
                    maxdistrib = max(bincounts(t).R);
                end
            end

            for t = 1:nbConditions
                subplot(nbRows, nbCols, t)
                histogram('BinCounts', bincounts(t).R, 'FaceColor', ...
                          colors(t, :), 'Facealpha', 0.4, 'BinEdges', BE)

                ylim([0, min(1.25*maxdistrib, 1)])
                legend(nameConditions{t, 1})
                xlabel('Andrews score', 'FontSize', 15)
                ylabel('Normalized number of microglia', 'FontSize', 15)
            end

            figure()
            for t = 1:nbConditions
                histogram('BinCounts', bincounts(t).R, 'FaceColor', ...
                          colors(t,:), 'Facealpha', 0.4, 'BinEdges', BE)

                hold on
                legend(nameConditions)
                xlabel('Andrews score', 'FontSize', 15)
                ylabel('Normalized number of microglia', 'FontSize', 15)
                hold on
                legend
            end
        end

    elseif strcmpi('No', answerDisplay)
        if answerFigures == 1 || answerFigures == 3
            path = uigetdir('', ['Select the folder containing the ' ...
                            'binarized individualized 300*300pi' ...
                            ' microglia']);

            % [files, nFiles, fileNo, ...
            %     fileNames, images] = Open_Microglia_Images(path, extension);
% 
            %%filePattern = fullfile(path, extension);
            %%theFiles = dir(filePattern);
            % fileNames = sortrows(fileNames,1,'ascend');
            % microgliaIdsSorted = fileNames(ranking);                          % Get microglia ranking with sort
            DisplayImages(path, extension, [], microgliaIdsSorted', [])
        end
        
        if answerFigures == 2 || answerFigures == 3
        
            minC = min(andrewsValuesMaxVar(:));
            maxC = max(andrewsValuesMaxVar(:));

            clear minCombined maxCombined
            %nb bins
            nx = numel(andrewsValuesMaxVar);
            
            % Calculate the 25th and 75th percentiles
            p25 = prctile(andrewsValuesMaxVar, 25);
            p75 = prctile(andrewsValuesMaxVar, 75);
            
            % Calculate the optimal bin width
            width = 2 * (p75 - p25) / nthroot(nx, 3);
            
            % Calculate the number of bins
            nbins = ceil((max(andrewsValuesMaxVar) ...
                    - min(andrewsValuesMaxVar)) / width);

            nbins = max(1, nbins);
            
            BE = linspace(minC, maxC, nbins);
            
            bincounts = histcounts(andrewsValuesMaxVar, BE, ...
                        'Normalization','probability');
            
            % Distributions display             
            figure()
            histogram('BinCounts', bincounts, 'BinEdges', BE)
            xlabel('Andrews score', 'FontSize', 15)
            ylabel('Normalized number of microglia', 'FontSize', 15)
        end
    end
end

%% Save figures & ranking
% Save the ranking
saving_path = uigetdir('','Select the folder for saving the ranking');
tableHeader = {'rank', 'andrews score', 'groups', 'microglia'};

% Set table values
classement = table((1:size(microgliaIdsSorted, 2))', ...
                    sortrows(andrewsValuesMaxVar, 1, 'ascend'), ...
                    nameConditions(classes(ranking)), ...
                    microgliaIdsSorted', 'VariableNames', tableHeader);

% Delete old output file if exists and save file
if exist([saving_path './ranking.xlsx'], 'file') == 2
    delete([saving_path './ranking.xlsx']);
end

writetable(classement, [saving_path '\ranking.xlsx'], ...
          'WriteRowNames', true);

% Save all figures
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = ['Fig' num2str(iFig)];
    saveas(FigHandle, fullfile(saving_path, [FigName '.png']));
    saveas(FigHandle, fullfile(saving_path, [FigName '.fig']));
end
