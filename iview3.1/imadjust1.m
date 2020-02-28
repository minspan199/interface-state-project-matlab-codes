% close all
clc
for i = 1:1:1
I = imread('t6.tif',i);
% I = Image14;
J = imadjust(I,[0.02 0.2],[0 0.5]);

% imshow(Image53)
figure
imshow(J)
% imwrite(J,'test1_intensified.tif','WriteMode','Append')
imwrite(J,'waveguide1.tif');
% T = imrotate(J,0.867);
% figure;imshow(T)
end