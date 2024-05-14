function inertia = inertia(image) 
    %%%% Computes the inertia of the image as the ratio between eigenvalues
    %%%% of covariance matrix of coordinates of image
    %%%%
    %%%% Inputs:
    %%%%    image: a binary image of a cell
    %%%%        
    %%%% Outputs:
    %%%%    inertia: inertia of the image

    [y, x] = find(image);                                                   % coordinates of image
    % N = length(x);                                                        % dimension of image
    y = y - mean(y);                                                        % centering coordinates
    x = x - mean(x);
    M = [y x]' * [y x];                                                     % computing the covariance matrix of coordinates
    [~, D] = eig(M);                                                        % computing the diagonal matrix of eigenvalues
    inertia = max(D(2, 2), D(1, 1)) / min(D(2, 2), D(1, 1));                                            % computing inertia

    if isnan(inertia)                                                       % if inertia is NaN, set inertia to 0
        inertia = 0;
    end   
end
