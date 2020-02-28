
function iViewMovieCbacks(cback) %% input callback switch

debug = 0;

global F Props Dat

datNames      = get(Props.view.tag.datNames,'String'); 
Props.isMovie = 1;

% Set gui init properties %
set(Props.movie.tag.path,'String',Props.movie.path)   

switch cback

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'writeImages' % % writes processed images (including FFT, etc.) to .bmp files (may change format later)

cd(Props.movie.path)
	
Props.isWriteImages   = 1;
for ii                = 1:length(Props.movie.files)

Props.movie.currFrame = Props.movie.files(ii);
iName                 = datNames{2+Props.movie.files(ii)};
             
% call iviewer to change to relevant data %
iViewerCbacks('datname');

cd(Props.movie.path) % go to relevant directory

if Props.isAx1
    if strcmp(Props.view.ax1Type,'ift')
        if get(Props.view.tag.isPhase,'value')
            datWrite = angle(Dat.ift.image);   
        else
            datWrite = abs(Dat.ift.image);
        end
    else
            datWrite = Dat.proc.image ;
     end
 else
            datWrite = Dat.proc.image ;
 end
 
%  datWrite = 255* imadjust(datWrite/max(datWrite(:)), stretchlim(datWrite,0), [0 1]);
% datWrite = datWrite - min(datWrite(:));
% datWrite = datWrite/max(datWrite(:)) * 64; % why 63-64 ?? why does the image turn black if saturated?

datWrite = datWrite - min(datWrite(:));
datWrite = datWrite/max(datWrite(:)) * 64; % why 64 ?? why does the image turn black if saturated?

if ~strcmp(iName(end-2:end),'bmp'), iName = [iName,'.bmp']; end

isJpg     = 0; % add this or other compressed format ... %
fileFmt   = 'bmp';
imwrite( ...
        datWrite ,...
        eval(Props.cmap),...
        [get(Props.movie.tag.movName,'string') ,' ', iName], ...
        fileFmt...
        )
    
end
 
Props.isWriteImages = 0;    

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'create' % % makes movie 

cd(Props.movie.path)	
	
Props.movie.figMov = figure;
set(gcf,'PaperPositionMode','auto')

if get(Props.movie.tag.isTight,'value')
    set(Props.movie.figMov,...
        'Position',[50 50 Props.roiInds(2) Props.roiInds(4)],...
        'Menubar','none'...
        );
end

 if get(Props.movie.tag.isFigSize,'value')
    figX = str2num(get(Props.movie.tag.figX,'string'));
    figY = str2num(get(Props.movie.tag.figY,'string'));   
    set(Props.movie.figMov,'Position',...
        [100 100 figX figY]...
        );
end

% background figure color %
set(Props.movie.figMov,'Color',Props.movie.fig.Color)
  
for ii = 1:length(Props.movie.files)

    
    Props.movie.currFrame = Props.movie.files(ii);
    iName                 = datNames{2+Props.movie.files(ii)};
             
    % call iviewer to plot relevant data %
    iViewerCbacks('datname');

 
    if get(Props.movie.tag.what,'value') == 2;
        figAx = Props.view.ax2;  %% the axes to make movie with
    else
        figAx = Props.view.ax1;
    end

    newPlot   = copyobj(figAx,Props.movie.figMov);
    axes(newPlot)

    set(newPlot,'Position',[0.1 0.1 0.8 0.8])

    if get(Props.movie.tag.isTight,'value')
        set(newPlot,'Position',[0 0 1 1]), end
    
    if ~strcmp(Props.view.ax2Type,{'profile','hist'})
        feval(F.aspectRatios)
        if Props.isCbar
            cb    = colorbar; 
            cbPos = get(cb,'position'); % set colorbar positions
            set(cb,'position',[cbPos(1:2) 0.01 cbPos(4)]);
        end
        colormap(Props.cmap)
    end
       
plotText(iName) % add text to plot


if get(Props.movie.tag.axOn,'value') == 1
    set(get(gca,'title'),'string',''), axis off    
end    

colormap(Props.cmap)

if get(Props.movie.tag.isFigset,'value')               % set axes with figset
    if exist('figset'), figset,  end
end

drawnow

          
if get(Props.movie.tag.isWriteImages,'value')	% write images from figure %
	set(gcf,'inverthardcopy','off')
	saveas(gcf,[iName,'.jpg'],'jpg')
else
	Frame(ii) = getframe(Props.movie.figMov);
end

clf
end % image loop

% compress movie %
if ~get(Props.movie.tag.isWriteImages,'value')	% write images from figure %
if ~debug 
    codecStr = get(Props.movie.tag.codec,'string');
    codec    = codecStr{get(Props.movie.tag.codec,'value')};
    movie2avi(Frame,get(Props.movie.tag.movName,'string'),...
             'compression',codec,...
             'fps',str2num(get(Props.movie.tag.fps,'string'))...
             );
end
end

close(Props.movie.figMov)    
figure(Props.movie.figH)
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'close' % % closes figure 

    close(gcbf);
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'view' % % test case   
    cd(        Props.movie.path)
    fileInfo = aviinfo(get(Props.movie.tag.movName,'string'));
    winopen(   fileInfo.Filename);
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'dirChange' % % change movie figure color
    fooPath     = uigetdir(Props.movie.path) ;
    if fooPath ~= 0 
        Props.movie.path = fooPath;
    end
    set(Props.movie.tag.path,'String',Props.movie.path)    
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'figColor' % % change movie figure color
    Props.movie.fig.Color = uisetcolor([0 0 0],'Figure Color');
    set(Props.movie.tag.figColor,'BackgroundColor',Props.movie.fig.Color)
        
case 'txtColor' % % change movie figure color
    Props.movie.txt.Color = uisetcolor([0 0 0],'Figure Color');
    set(Props.movie.tag.Color,'BackgroundColor',Props.movie.txt.Color)    
	Props.movie.txt.FontSize	= str2num(get(Props.movie.tag.FontSize,  'string')); 
% 	Props.movie.txt.FontWeight	=         get(Props.movie.tag.FontWeight,'string') ; 
	
end

% set data based on previous settings %

% Set figure properties %
% Props.movie = movie;

Props.isMovie = 0;


function plotText(iName)

global F Props Dat

% create plot text %
[txt1,txt2,txt3] = deal('');
if get(Props.movie.tag.isText1,'value')
    txt1 = get(Props.movie.tag.text1,'string');
end
        
if get(Props.movie.tag.isText2,'value')
    ind1     = str2num(get(Props.movie.tag.tokInd1,'string'));
    ind2     = str2num(get(Props.movie.tag.tokInd2,'string'));
    tokenStr = get(Props.movie.tag.token,'string');
    tokenInd = findstr(iName,tokenStr);
    tokenInd = tokenInd + length(tokenStr);
    txt2     = iName(tokenInd+ind1:tokenInd+ind2);
    
        %%%%%%%%%%%%%%%%%%%
% Insert code to operate on Token here %        





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
        % % % % % % % % % % % % %
 end
    
 if get(Props.movie.tag.isText3,'value')
    txt3 = get(Props.movie.tag.text3,'string');
 end
            
txt      = sprintf([txt1,txt2,txt3]);    
txtHndl  = text(str2num(get(Props.movie.tag.textX,'string')),...
               str2num(get(Props.movie.tag.textY,'string')),...
               txt);
set(txtHndl,Props.movie.txt)    












