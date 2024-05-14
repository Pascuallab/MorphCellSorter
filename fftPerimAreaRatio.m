function [perimAreaRatio,cellPerimeter] = fftPerimAreaRatio(image, resolution) 
    %%%% Computes the decimal logarithm of the ratio between the squared
    %%%% perimeter and area of the cell in image. The perimeter of the cell
    %%%% is detected using fast fourier transforms and filtering out low
    %%%% frequencies.
    %%%%
    %%%% Inputs:
    %%%%    image: a binary image of a cell
    %%%%    resolution: resolution of the images, in micrometers per pixel
    %%%%        
    %%%% Outputs:
    %%%%    perimAreaRatio: decimal logarithm of the ratio between the
    %%%%                    squared perimeter and area of the the cell
    %%%%    cellPerimeter: perimeter of the cell (microns)

    Xfft = fftshift(fft2(image));                                           % centered 2d fast Fourier transform (fft2) of image
    Xphi = imag(log(Xfft));                                                 % imaginary part of natural logarithm of the fft2 (extract phase or frequencies)
    Xphiinv = abs(ifft2(exp(1i * Xphi)));                                   % inverse fft2 of the phases (
    Q = (imdilate(image, ones(3, 3)) - imerode(image, ones(3, 3)));         % intersection between once dilated and once eroded images of the cell to detect borders
    W = Q .* (Xphiinv > max(Xphiinv(:)) / 4);                               % filtering borders using high frequencies (contour)
    cellPerimeter=sum(W(:)) * resolution;
    perimAreaRatio = (cellPerimeter) .^ 2 ./ ...                            % compute ratio
                     (sum(image(:)) * resolution ^ 2);

    perimAreaRatio = log(perimAreaRatio) / log(10);                         % compute decimal logarithm of ratio
end
