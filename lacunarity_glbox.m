function [lacunaritySlope, lacunarityMean] = lacunarity_glbox(image) 
    %%%% Gliding box lacunarity algorithm based on the ideas presented in
    %%%% the paper by Tolle et al., Physica D, 237, 306-315, 2008
    %%%% 
    %%%% Returns the lacunarity slope dans mean at various box lengths
    %%%% (edge sizes are multiples of 2)
    %%%%
    %%%% Author: Tegy J. Vadakkan
    %%%% Date: 09/08/2009
    %%%% 
    %%%% Input:
    %%%%    image: is the binary file with 0's and 1's. 0's represent the
    %%%%           holes or the lacunae (must be square)
    %%%% 
    %%%% Output: 
    %%%%    lacunaritySlope: slope of the lacunarity curve
    %%%%    lacunarityMean: mean lacunarity at various box lengths
    
    [rows, ~] = size(image);                                                % first dimension of the image
    image = 1 - image;                                                      % inverting the image
    q = 2;                                                                  % initialization of the index q (power of 2)
    
    while q <= rows                                                         % as long as q is inferior to the image dimension, scan image with boxes of increasing size 
        nn = q - 1;                                                         % box dimension
        rnn = rows - nn;                                                    % number of boxes of dimension nn
        index = uint8(log2(q));                                             % index = 1,2,3,4 etc...                                    
        count(index) = power(rnn, 2);                                       % saving current box dimension
        sigma(index) = 0.0;                                                 % initializing 
        sigma2(index) = 0.0;

        for i = 1:rnn
            for j = 1:rnn                                                   % for each box of dimension nn
                sums = sum(sum(image(i:i + nn, j:j + nn)));                 % compute the sum of pixels (lacunarities)
                sigma(index) = sigma(index) + sums;                         % compute sigma, sum of lacunarities for each box    
                sigma2(index) = sigma2(index) + power(sums,2);              % compute sigma2 sum of squared lacunarites for each box
            end
        end
        q = q * 2;                                                          % increase counter q
    end
    
    for i = 1:index                                                         % for each box dimension
        M(i, 1)= count(i) * sigma2(i) / power(sigma(i), 2);                 % computing lacunarity index M based on sigma and sigma2
    end

    L = polyfit(log(count)', log(M), 1);                                    % regression on log(M) to extract the lacunarity slope and mean
    lacunaritySlope = L(1);
    lacunarityMean = mean(M);

end
