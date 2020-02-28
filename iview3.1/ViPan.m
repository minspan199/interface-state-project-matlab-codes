srcFiles = dir('D:\Camera\170514\H-SLD\Filtered\*.tif');
k = 0;
mov = VideoWriter('Movie1.avi');
mov.FrameRate = 20;
open(mov)
for i = 1 : length(srcFiles)
    TFilename = strcat('D:\Camera\170514\H-SLD\Filtered\Filtered_',num2str(i+10),'.tif');
    writeVideo(mov,imread(TFilename))
end
close(mov);



function Movie_from_frames(name,filetype,number_of_frames,display_time_of_frame)
% mov = avifile('MOVIE.avi');
count=0;
name1=strcat(name,num2str(1),filetype);
    mov=imread(name1);
for i=1:number_of_frames
    name1=strcat(name,num2str(i),filetype);
    a=imread(name1);
    while count<display_time_of_frame
        count=count+1;
        imshow(a);
        F=getframe(gca);
        mov=addframe(mov,F);
    end
    count=0;
end
close all
mov=close(mov);
end
    