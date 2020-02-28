% srcFiles = dir('D:\Google Drive Buffalo\Matlab Codes\170426\C\Filtered\*.tif');
% for i = 1 : length(srcFiles)
%   TFilename = strcat('D:\Google Drive Buffalo\Matlab Codes\170426\C\Filtered\',srcFiles(i).name);
%   I = imread(strcat('Image1_',num2str(i),'.tif'));
%   % eval(['Image',num2str(i),'=','imread(TFilename);']);
%   imwrite(I,strcat('Image1_',num2str(i+200),'.tif'))
% %   F{i} = I;
% end
% % clear i TFilename srcFiles
% srcFiles = dir('D:\Google Drive Buffalo\Matlab Codes\170426\C\Filtered\*.tif');
for i = 1 : 400
%   TFilename = strcat('D:\Google Drive Buffalo\Matlab Codes\170426\C\Filtered\',srcFiles(i).name);
%   I = imread(strcat('Image5_',num2str(i),'.tif'));
  % eval(['Image',num2str(i),'=','imread(TFilename);']);
  
  I = imread('test1.tif',i);
%   imwrite(I,strcat('Image4_',num2str(i+200),'.tif'))
I = imadjust(I,[0 0.3],[]);
imwrite(I,'test1_intesified.tif','WriteMode','append');
%   F{i} = I;
end
% clear i TFilename srcFiles