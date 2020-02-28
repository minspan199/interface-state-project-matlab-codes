% v = VideoReader('3darraydefect2_monitor.mpg');
% 
% 
% for img = 1:v.NumberOfFrames
%     filename=strcat('frame.tif'); 
%     b{img} = readFrame(v); 
%     imwrite(b{img},filename,'WriteMode','append'); 
% end

clc
close all
clear all
source='3darrayloss2_monitor.mpg';
vidobj=VideoReader(source);
frames=vidobj.Numberofframes;
for f=1:frames
    thisframe=read(vidobj,f);
    figure(1);imagesc(thisframe);
    for K=1:size(thisframe,3)
        itframe=thisframe(:,:,K);
        thisfile = sprintf('frame_%04d.jpg', K);
        imwrite( thisframe, thisfile );
    end
end

