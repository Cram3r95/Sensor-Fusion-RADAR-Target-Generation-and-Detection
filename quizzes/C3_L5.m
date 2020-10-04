% Class 3. Range-Doppler Estimation - Lesson 5. The 2D FFT
% Carlos Gómez Huélamo

% 2-D Transform
% The 2-D Fourier transform is useful for processing 2-D signals and other 2-D data such as images.
% Create and plot 2-D data with repeated blocks.

P = peaks(20); % peaks is a function of two variables, obtained by translating and scaling Gaussian distributions,
% which is useful for demonstrating mesh, surf, pcolor, contour, and so on.
X = repmat(P,[5 10]); % repmat repeat copies of array (in this case, 5 rows and 10 columns)
% imagesc(X) % Display image with scaled colors

% TODO : Compute the 2-D Fourier transform of the data.  
% Shift the zero-frequency component to the center of the output, and 
% plot the resulting 100-by-200 matrix, which is the same size as X.

% Compute the 2-D Fourier transform of the data. Shift the zero-frequency component to the center of the output, 
% and plot the resulting 100-by-200 matrix, which is the same size as X

% Y = fft2(X);
% imagesc(abs(fftshift(Y)))

% Pad X with zeros to compute a 128-by-256 transform

Y = fft2(X,2^nextpow2(100),2^nextpow2(200));
imagesc(abs(fftshift(Y)))