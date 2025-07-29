function [] = filtering()

% RGB = imread('glories.png');
RGB = imread('airbos_f7_p5.jpg');
G = double(rgb2gray(RGB));

clf
figure(1)
imshow(G / max(G(:)))
colormap gray
axis equal
title('Initial')

figure(2)
r = 20;
[B, PSF] = blur(G, 2);

B = medfilt2(B,[r r]);

% B = blur(B, 2);

imshow(B / max(B(:)))
colormap gray
axis equal
title(['Blurred, radisu: ', num2str(r)])

B = double(B);
B(:) = 255 * (B(:) - min(B(:))) / (max(B(:)) - min(B(:)));

% imwrite(B,'med.gif')
% imwrite(B,'schlieren.gif')
save('schlieren.mat', 'B');
stop





% figure(3)
% maxDeconvIterations = 10;
% sharper = deconvlucy(B,PSF,maxDeconvIterations);
% image(sharper)
% colormap gray
% axis equal
% title(['Deblurred, maxIterations: ' num2str(maxDeconvIterations)])

figure(4)
E = edge(B,'canny');
imshow(E)
% colormap gray
axis equal
title('Edge detection')


% W = fft2(G);
% Y = fftshift(W);
% 
% figure(2)
% colormap gray
% image(log(abs(Y)))
% 
% Y2 = lowPass(Y,70);
% figure(3)
% colormap gray
% image(log(abs(Y2)))
% 
% W2 = ifftshift(Y2);
% G2 = ifft2(W2);
% figure(4)
% colormap gray
% image(abs(G2))

end

function Y = lowPass(Y0, radius)
    [m,n] = size(Y0);
    
    [i, j] = ind2sub([m n], 1:(m*n));
    
    w1 = m / 2;
    w2 = n / 2;
    
    select = ((i-w1).^2 + (j-w2).^2) < (radius^2);
    
    Y = Y0 * 0;
    
    Y(select) = Y0(select);
end

function [blurred, PSF] = blur(I, radius)
    PSF = fspecial('gaussian',radius, radius);
    blurred = imfilter(I,PSF,'symmetric','conv');
end

