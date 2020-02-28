

function iVFuns

global F     % declares the structure of function handles as global
              % in scope, used by any function 

F = struct('createTags',@createTags,...
	        'refreshTags',@refreshTags,...
	        'findHndl',@findHndl,...
	        'localPlot',@localPlot,...
            'axset',@axset,...
            'centroid', @centroid,...
            'clearIndicators', @clearIndicators,...
            'localFWHM', @localFWHM,...
            'convertUnits', @convertUnits,...
            'pix2mm', @pix2mm,...
            'mm2pix', @mm2pix,...
            'pix2nm', @pix2nm,...
            'nm2pix', @nm2pix,...
            'nextPlot', @nextPlot,...
            'drawProfLines', @drawProfLines,...
            'catProfile', @catProfile,...
            'colord', @colord,...
            'colortable', @colortable,...
            'nearest',@nearest,...
            'aspectRatios',@aspectRatios,...
            'cbarsSet',@cbarsSet,...
            'newfig',@newfig,...
            'DesignFilter',@DesignFilter,...
            'shift',@shift,...
            'computeFFT',@computeFFT,...
            'computeIFFT',@computeIFFT,...
			'background_subtract',@background_subtract,...
			'prpl',@prpl,...
			'xsection',@xsection...
			  );

		  
		  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % subfunctions
% % % %

% Create structure array of handles accessed by tags %
function Out = createTags(hndls,tags)

% *** same as built in function guihandles !!! 

for ii = 1:length(tags)
	iTag = tags{ii};
	if ~isempty(iTag)	
		Out.(iTag) = hndls(ii);
	end		
end

% updates all the indicators with data in relevant structures %
function refreshTags(gui,Fields)

global Props

for jj = 1:length(Fields)
    jField = Fields{jj};
    nm = fieldnames(Props.(gui).(jField));
    dummy = nm; % ?? weird bug if not use this
    for ii = 1:length(nm)
        iNm = dummy{ii};

        hndl = Props.(gui).tag.(iNm);
        prop = Props.(gui).(jField).(iNm);
        if isstr(prop)
            propType = 'string'; 
        else
            propType = 'value';
        end
        set(hndl,propType,prop)
    end
    
end

% find and return handle to specific tag %
function Hndl = findHndl(gui,tag)
global Props

Hndl = Props.(gui).Hndls(find(strcmp(Props.(gui).tags,tag)));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Plotting Function 										  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function localPlot(plotType) %% use this so that you can replot if
                                                           % if chng images
global Props Dat

switch plotType
case 'orig'
    axes(Props.view.ax2), cla, set(gca,'visible','off'), title('')
    axes(Props.view.ax1); set(gca,'visible','on')
    imagesc(Dat.orig.x,Dat.orig.y,Dat.orig.image); 
    title('Original Image'); 
        
    % % set user data and Hndls
    Props.view.ax1     = gca;
    Dat.proc           = Dat.orig;
    Props.roi          = [min(Dat.orig.x) max(Dat.orig.x) ...
                          min(Dat.orig.y) max(Dat.orig.y)]; % ROI in current units
    Props.view.ax1Type = 'orig'; 
    Props.isAx1        = 1;
    axset  
        
case 'roi'
	
    if Props.isAx1
        axes(Props.view.ax1)
        imagesc(Dat.proc.x,Dat.proc.y,Dat.proc.image); 
        Props.view.ax1 = gca;
        Props.view.ax1Type = 'roi'; 
    else
        axes(Props.view.ax2)
        if strcmp(Props.view.ax2Type,'ft')
            datPlot = log10(abs(Dat.proc.image));
            
        else
            datPlot = Dat.proc.image;
        end
        imagesc(Dat.proc.x,Dat.proc.y,datPlot);
        Props.view.ax2 = gca;
    end
    
    title('ROI of Original Image') 
       
    % % set user data and Handls
    axset 
    
case 'contour'
    axes(Props.view.ax2);  set(gca,'Visible','on')
    contour(Dat.proc.x,Dat.proc.y,double(Dat.proc.image),Props.ncontours)
    set(gca, 'YDir', 'reverse'); 
    title('Contour Plot'); 

    % % set data and handles
    Props.view.ax2Type = 'contour';
    Props.iview.ax2    = gca;
    Props.isAx1        = 0;
    axset 
    
case 'plot3'
	%% fix this to open in new window
	%% matlab renders really slowly in iView window, not so bad
	% in new one... just need to keep careful track of the axes to 
	% keep from getting lost on the axes handles
	
	%if ~get(Props.view.tag.is3DnewFig,'value')
	%figure
	%Props.isAx1        = 1; 		
	
	%else
	axes(Props.view.ax2), set(gca,'Visible','on')
	% % set data and handles
	Props.view.ax2Type = 'plot3'; 
    Props.view.ax2     = gca;
    Props.isAx1        = 0; 
    axset 
				
	%end
	
	if get(Props.view.tag.surf,'value')
        surf(Dat.proc.x,Dat.proc.y,double(Dat.proc.image)), shading interp
    else
        mesh(Dat.proc.x,Dat.proc.y,double(Dat.proc.image))
    end
    
	zlabel('Intensity'), title('3-D Intensity Profile'); axis tight
	

case 'ft'
    axes(Props.view.ax2); set(gca,'Visible','on')
    imagesc(Dat.proc.x,Dat.proc.y,abs(log10(Dat.proc.image))); title('Fourier Transform [log]'); 
		  
    % % set data and Hndls
    Props.view.ax2Type = 'ft'; 
    Props.view.ax2     = gca ;
    Props.isAx1        = 0   ;   
    axset  
    
case 'ift'
    axes(Props.view.ax1);
    if get(Props.view.tag.isPhase,'value')
        imagesc(Dat.proc.x,Dat.proc.y,angle(Dat.ift.image)); title('Phase of  Inverse Fourier Transform [rad]'); 
    else
        imagesc(Dat.proc.x,Dat.proc.y,abs(Dat.ift.image));   title('Magnitude of Inverse Fourier Transform'); 
    end
		  
    % % set data and Hndls
    Props.view.ax1Type = 'ift'; 
    Props.view.ax1     = gca;
    Props.isAx1        = 1; 
    axset 
    
case 'inten'
    axes(Props.view.ax2), set(gca,'Visible','on')
    imagesc(Dat.proc.x,Dat.proc.y,Dat.proc.image,Dat.proc.Irng)
    title(['Intensity Range [',num2str(Dat.proc.Irng),']']);
    
    % % set dat and handles
    Props.view.ax2Type = 'inten';
    Props.view.ax2     = gca;
    Props.isAx1        = 0;
    axset 

case 'profile'
    
    axes(Props.view.ax2)
    nLines    = size(Dat.profile.dat,2); 
    lProfHndl = plot(Dat.profile.xDat,Dat.profile.dat);  
    xLabStr   = 'Position in profile';
    if Props.isMm
        unitLab = ', mm'; else unitLab = ', pixels'; end
    if Props.isNm, xLabStr = 'Wavelength';
        unitLab = ', nm'; end
        
    ylabel('Intensity')       , title('Image Profile');     
    xlabel([xLabStr,unitLab]) , Props.view.ax2Type = 'profile';    
    Props.view.ax2 = gca;
		  
    axes(Props.view.ax1) 
    tit = get(get(gca,'title'),'string'); 
    
    if strcmp(Props.view.ax1Type,'ift'), 
		if get(Props.view.tag.isPhase,'value')
			datFoo = angle(Dat.ift.image);
		else
			datFoo = abs(Dat.ift.image);
		end
    else
        datFoo = Dat.proc.image;
    end
    
    Props.iview.ax1 = imagesc(Dat.proc.x,Dat.proc.y,datFoo); 
    
    title(tit), axset 
        
    drawProfLines(Dat.profile);
		  
    for ii = 1:length(lProfHndl)
        set(lProfHndl(ii),'color',colord(ii)), end
        
case 'hist'
        
    % % get data, plot and label    
    dat 		= double(Dat.proc.image);
    axes		(Props.view.ax2);
	set			(gca,'Visible','on')
    [counts,x] 	= hist(dat(:),256); %% only 256 bins here ... can alter
    histHnd		= plot(x,counts);
% 	set(histHnd,'FaceColor',colord(1))
    title('Histogram of Current Image'); ylabel('Occurrences'); xlabel('Intensity')
    axis([min(x) max(x) min(counts) max(counts)])
    
    % % set data and handles
    Dat.proc.hist.x    = x;      Dat.proc.hist.counts = counts;
    Props.view.ax2Type = 'hist'; Props.view.ax2       = gca; 
    
    set(Props.view.tag.lin,'Value',1), set(Props.view.tag.log,'Value',0) 
    axes(Props.view.ax1)
    
otherwise
    set(Props.cmdHndl,'string','Error: No other kind of plot available'), return
end

aspectRatios %% set plot aspect ratios
    
if strcmp(plotType,'ft')
    set(Props.view.tag.lin,'Value',0), set(Props.view.tag.log,'Value',1)
else
    if ~Props.isAx1
        if strcmp(Props.view.ax2Type,'ft')
            set(Props.view.tag.lin,'Value',0), set(Props.view.tag.log,'Value',1)
        end
    end
    set(Props.view.tag.lin,'Value',1), set(Props.view.tag.log,'Value',0)
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Axse Function 											  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function axset

global Props

% % define labels depending on units % %
if get(Props.view.tag.isMm,'value')
    xLab = 'position, mm'; yLab = xLab;
else, xLab = 'pixels'; yLab = xLab; end

if get(Props.view.tag.isNm,'value')
	xLab = 'wavelength, nm'; end

if ~strcmp(Props.view.ax2Type,'ft') %% set axis labels
    set(get(gca,'xlabel'),'string',xLab)
    set(get(gca,'ylabel'),'string',yLab)
end

set(Props.view.tag.isAx1,'Value',Props.isAx1) % % set dots and ax
set(Props.view.tag.isAx2,'Value',~Props.isAx1)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% newfig													  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function newfig %

global Props

if Props.isAx1 
  if strcmp(Props.view.ax2Type,'profile') | strcmp(Props.view.ax2Type,'hist')
      figAx = Props.view.ax2; 
  else figAx = Props.view.ax1; 
  end
else figAx = Props.view.ax2; end
    
newFig = figure; 
if get(Props.view.tag.colordef,'value')
	set		(gcf,'color','k'), end

newPlot = copyobj(figAx,newFig);
axes(newPlot)

set(newPlot,'Position',[0.1 0.1 0.8 0.8])

if ~strcmp(Props.view.ax2Type,{'profile','hist'})
    aspectRatios
    if Props.isCbar
        cb = colorbar; 
        cbPos = get(cb,'position'); % set colorbar positions
        set(cb,'position',[cbPos(1:2) 0.01 cbPos(4)]);
    end
    colormap(Props.cmap)
end

% % % %
function [xCom,yCom] = centroid(x,y,dat) %% conventional algorithm
dat  = double(dat);
mass = sum(sum(dat)); 

xArray = repmat(x,size(dat,1),1);
yArray = repmat(y',1,size(dat,2));

xMom = sum(sum(dat.*xArray));            %% x and y moments
yMom = sum(sum(dat.*yArray));

xCom = xMom/mass; yCom = yMom/mass;

% % % %
function clearIndicators(varargin) % % clears text boxes

global Props

if nargin == 0, type = 'all';
elseif nargin == 1, type = deal(varargin{:}); end

% % gets handles for text boxes % %
txtTags = {'max','yFWHM','xFWHM','xCOM','yCOM',...
           'minOrig','sumOrig','avgOrig',... % % original props
           'minProc','sumProc','avgProc',...
           'maxProc','yFWHMProc','xFWHMProc','xCOMProc','yCOMProc',...
           'profFWHM','profMax','xCoord','yCoord','ICoord'};

switch type
    case 'all' , iStart = 1;
    case 'orig', iStart = 9;
end

iEnd = size(txtTags,2);

for ii = iStart:iEnd
    txtHndl.(txtTags{ii}) = findobj(Props.view.fig,'Tag',txtTags{ii});
    iHndl = findobj(Props.view.fig,'Tag',txtTags{ii});
    set(iHndl,'String','-');
end

% % computes full width % %
function [sFWHM,sMax,d1,d2] = localFWHM(imSlice,x)

if nargin < 2, x = 1:length(imSlice); end

[sMax,iM] = max(imSlice);

srng = x(imSlice>=sMax/2);
if ~isempty(srng)
    d2 = srng(end);  d1 = srng(1);
    sFWHM = abs( srng(end) - srng(1) );
else
    sFWHM = NaN;
end

if sFWHM == length(imSlice(~isnan(imSlice))) - 1
    sFWHM = NaN;
end


% % % %
function [x,y,units] = convertUnits(xPix,yPix) %% convert to proper units (use pixels as input)

global Props

if Props.isMm %% should end in mm
    x = pix2mm('x',xPix,Props.pixSize);
    y = pix2mm('y',yPix,Props.pixSize);
    units.x = 'mm'; units.y = 'mm';
else
    x = xPix; y = yPix;
    units.x = 'pixels'; units.y = 'pixels';

end

if Props.isNm
     x = pix2nm(xPix);
     units.x = 'nm';
end

% % % %
function mm = pix2mm(dir,pix,pixSize)

global Props

mmX = str2num(get(Props.view.tag.mmX,'string'));
mmY = str2num(get(Props.view.tag.mmY,'string'));

if isnan(mmX) | isnan(mmY)
    set(Props.view.tag.cmd,'string','Invalid image size value'), return
end

switch dir
    case 'x'
        mm = mmX * pix / pixSize(2);  
    case 'y'
        mm = mmY * pix / pixSize(1);
end

% % % %
function pix = mm2pix(dir,mm,pixSize)

global Props

mmX = str2num(get(Props.view.tag.mmX,'string'));
mmY = str2num(get(Props.view.tag.mmY,'string'));

if isnan(mmX) | isnan(mmY)
    set(Props.view.tag.cmd,'string','Invalid image size value'), return
end

switch dir
    case 'x'
        pix = mm / mmX * pixSize(2);  
    case 'y'
        pix = mm / mmY * pixSize(1);
end


% % % %
function nm = pix2nm(xCoord)

global Props

xPix = str2num(get(Props.view.tag.nmPix,'string'));
xCal = str2num(get(Props.view.tag.nmCal,'string'));
xDisp = str2num(get(Props.view.tag.nmDisp,'string'));

if isnan(xPix) | isnan(xCal) | isnan(xDisp)
    set(Props.view.tag.cmd,'string','Invalid wavelength calibration value'), return
end

nm = xCal + (xCoord - xPix) * xDisp;

% % % %
function pix = nm2pix(nm)

global Props

xPix  = str2num( get(Props.view.tag.nmPix, 'string') );
xCal  = str2num( get(Props.view.tag.nmCal, 'string') );
xDisp = str2num( get(Props.view.tag.nmDisp,'string') );

if isnan(xPix) | isnan(xCal) | isnan(xDisp)
    set(Props.view.tag.cmd,'string','Invalid wavelength calibration value'), return
end
pix =  xPix + (nm - xCal) / xDisp;

% % To plot the next image in similar fashion to that currently displayed % %
function nextPlot(DatPrev)

global Props Dat

clearIndicators('all')    

[Dat.orig.x,Dat.orig.y,Dat.units] = convertUnits(Dat.orig.xPixels,Dat.orig.yPixels);
Props.pixSize                     = size(Dat.orig.image); %% 

if (Props.roiInds(2) < size(Dat.orig.image,2) ) ...
 & (Props.roiInds(4) < size(Dat.orig.image,1) )

    Dat.proc.image = Dat.orig.image(Props.roiInds(3):Props.roiInds(4),... % Set new range ... could
                                    Props.roiInds(1):Props.roiInds(2)...  % also change to IMCROP          
                                   );
    Dat.proc.x = Dat.orig.x(Props.roiInds(1):Props.roiInds(2));
    Dat.proc.y = Dat.orig.y(Props.roiInds(3):Props.roiInds(4));
    
else
    Dat.proc = Dat.orig;                % case when new image too small for ROI
end
    
% change up this order... redundant %%%%%%%%%%%%%%%
fooAx1Type         = Props.view.ax1Type;

axes(Props.view.ax1)
Props.isAx1 = 1;

if ~Props.isWriteImages % don't do this if writing images
    if strcmp(Props.view.ax1Type,'ift'),
% % % %         shift 
% % % %         computeIFFT
% % % %         localPlot(Props.view.ax1Type); 
    else   
        localPlot('roi');                                  % plot axes 1 as roi
    end %% plot axes 1 end  % new version- check out %

end

Props.view.ax1Type = fooAx1Type;
Props.view.ax1     = gca;
cbarsSet('view')

if strcmp(Props.view.ax2Type,'ft'),  computeFFT, end   % new version- check out 
  
if strcmp(Props.view.ax1Type,'ift'),
    shift 
    computeIFFT
    if ~Props.isWriteImages % don't do this if writing images
        localPlot(Props.view.ax1Type); 
		cbarsSet('view')
    end
end

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(DatPrev.profile)
    Dat.profile.dat      = [];  Dat.profile.coords.lims = [];
    Dat.profile.fwhm     = [];  Dat.profile.max         = [];
    Dat.profile.coords.x = [];  Dat.profile.coords.y    = [];
    Dat.profile.dist     = [];  Dat.profile.xDat        = [];
    
    for ii = 1 : size(DatPrev.profile.coords.lims,1)
        lim(1:2) = DatPrev.profile.coords.lims(ii,1:2);
        lim(3:4) = DatPrev.profile.coords.lims(ii,3:4);
        
        [xCoords,yCoords,imProf,xLim,yLim] = improfile(Dat.proc.x,Dat.proc.y,Dat.proc.image,...
                                                       lim(1:2),lim(3:4),...
                                                       Props.nProfPoints,'bilinear'...
                                                       );     

        Dat.profile = catProfile(Dat.profile,xCoords,yCoords,imProf,xLim,yLim);
    end
    localPlot('profile'); 
   
    set(Props.view.tag.profMax, 'string',...
        num2str(Dat.profile.max(end)))    
    set(Props.view.tag.profFWHM,'string',...
        num2str(Dat.profile.fwhm(end)))  
	
else
	
	if ~isempty(DatPrev.Irng) %% set range Intensity range if specified
     Dat.proc.Irng = DatPrev.Irng; 
	else Dat.proc.Irng = [min(Dat.proc.image(:)) max(Dat.proc.image(:))];
	end

	if ~strcmp(Props.view.ax2Type,'none') %% plot axes 2
 	   if ~Props.isWriteImages % don't do this if writing images
  	      axes(Props.view.ax2)
   	      Props.isAx1    = 0;
    	  localPlot(Props.view.ax2Type);
     	  cbarsSet('view')
		  Props.view.ax2 = gca;
	   end
   end
end

set(Props.view.fig,'Name',['Image= ',Dat.fname])
set(Props.view.tag.datNames,'Value',Props.currImage + 2)    

% % % %
function drawProfLines(profile)

global Props
nLines = size(profile.dat,2);

hold on
for iProf = 1:nLines % % plot line colors one at a time
	xc              = profile.coords.lims(iProf,1:2); 
    yc              = profile.coords.lims(iProf,3:4);
    l1              = plot(xc,yc,'o-'); 
	l2              = plot(xc,yc,'wo:');  
    l3              = plot(xc,yc,'wx');
	Props.profHndl  = cat(1,Props.profHndl,[l1;l2;l3]);
	iCol            = colord(iProf);
    
	set(l1, 'markerfacecolor',iCol,        'color',iCol )
	set(Props.view.tag.profMax,  'foregroundcolor',iCol )    
	set(Props.view.tag.profFWHM, 'foregroundcolor',iCol )
end
hold off

% % % %
function profile = catProfile(profile,xCoords,yCoords,imProf,xLim,yLim)

global Props 

profile.xUnits = 'pixels'; profile.yUnits = 'pixels';

% % total position vector
xDat = sqrt((xCoords-xLim(1)).^2 + (yCoords-yLim(1)).^2);

if Props.isMm, profile.xUnits = 'mm'; profile.yUnits = 'mm'; end

if Props.isNm
    xDat = xCoords; %% total position vector ... dy is 0 in this case, want 
                     % absolute x position, not relative
    profile.xUnits = 'nm';
end

% % FWHM for improfile 
[profFWHM,profMax] = localFWHM(imProf,xDat);

lims = [xLim(:)' yLim(:)'];

profile.coords.lims = cat(1, profile.coords.lims , lims    ); 
profile.dat         = cat(2, profile.dat         , imProf  );    
profile.max         = cat(2, profile.max         , profMax );
profile.fwhm        = cat(2, profile.fwhm        , profFWHM);
profile.xDat        = cat(2, profile.xDat        , xDat    );
profile.coords.x    = cat(2, profile.coords.x    , xCoords );
profile.coords.y    = cat(2, profile.coords.y    , yCoords );

% % % %
function [col,name,colOrder] = colord(in) %% sets the colororder (see colord.m)

if ~isstr(in)
   colOrder = 'iauvqkfdnberljzsptmcghy';
   num = rem(in,length(colOrder));
   if ~num, num = num+length(colOrder); end	%% remainder zero set to last value
   n = colOrder(num);
else
   n = in;
end

[col,name] = colortable(n);

% % % %
function [col,name]=colortable(n) %% called by colord
vc='kymcrgbadefhijlnpqstuvz';       % valid color codes
%   k               y               m               c
cn=[0.00,0.00,0.00; 1.00,1.00,0.00; 1.00,0.00,1.00; 0.00,1.00,1.00;
%   r               g               b               a
    1.00,0.00,0.00; 0.00,1.00,0.00; 0.00,0.00,1.00; 0.00,0.70,0.00;
%   d               e               f               h
    0.40,0.40,0.40; 0.00,0.40,0.00; 0.90,0.00,0.40; 1.00,0.80,0.00;
%   i               j               l               n
    0.00,0.00,0.70; 0.20,0.80,0.50; 0.80,0.40,0.80; 0.50,0.20,0.00;
%   p               q               s               t
    1.00,0.40,0.60; 1.00,0.40,0.00; 0.00,0.80,1.00; 0.80,0.40,0.00;
%   u               v               z              w no white...
    0.70,0.00,0.00; 0.60,0.00,1.00; 0.60,0.60,0.60; 1.00 1.00 1.00];
                     % plot color table

name={'black','yellow','magenta','cyan',...
      'red','green','blue','apple green',...
      'dark gray','evergreen','fuchsia','honey',...
      'indigo','jade','lilac','nutbrown',...
      'pink','kumquat','sky blue','tan',...
      'umber','violet','zinc','white'};

ind = findstr(vc,n);
col = cn(ind,:);
name = name(ind);

% % finds the nearest value in a vector to the given points (see .m NEAREST)
function NearestValue = nearest(Numbers,PinPoints)

if nargin < 2,NearestValue = []; return, end

Numbers = real(Numbers); Pins = length(PinPoints);
Max = max(Numbers); Min = min(Numbers);

for iPin = 1:Pins
   if PinPoints(iPin) > Max, NearestValue(iPin) = Max;
   elseif PinPoints(iPin) < Min, NearestValue(iPin) = Min;
   else
      nul = abs(Numbers - PinPoints(iPin));
      iValue = find(nul==min(nul));
      NearestValue(iPin) = Numbers(iValue);
   end
end

% % Sets data aspect ratios % %
function aspectRatios

global Dat Props

dx = Dat.proc.x(end) - Dat.proc.x(1); %% values used in all
dy = Dat.proc.y(end) - Dat.proc.y(1); 

if ~Props.isNm
    if Props.isAx1
        daspect([1,1,1]) %% good for both pixels and mms ... 1:1
    else
        if ~strcmp(Props.view.ax2Type,{'hist','profile','plot3'})
            daspect([1,1,1])
        elseif strcmp(Props.view.ax2Type,'plot3')
            xmax = max([dy,dx]); zmax = max(double(Dat.proc.image(:)));
            daspect(Props.view.ax2,[1 1 zmax/xmax])

        end
    end
else %% setting aspect ratios with wavelength present- need pix, mm support
      
    %% need to adjust for the mm values
    if Props.isMm
%           dydxPix = dy/length(Dat.proc.xPixels);       
        dydxPix = length(Dat.proc.y)/length(Dat.proc.xPixels);        
    else        
        dydxPix = length(Dat.proc.yPixels)/length(Dat.proc.xPixels);        
    end
    dydx = dy/dx;
    
    if Props.isAx1 %% set to pixel aspect ratio in cases w/ wavelength
        daspect([dx*dydxPix,dy,1])
    else
        if ~strcmp(Props.view.ax2Type,{'hist','profile','plot3'})
            daspect([dx*dydxPix,dy,1])
        elseif strcmp(Props.view.ax2Type,'plot3')
            xmax = max([dy,dx*dydxPix]); zmax = max(double(Dat.proc.image(:)));
            daspect(Props.view.ax2,[dx*dydxPix,dy, zmax])
        end
    end
end

% % set colorbars widths
function cbarsSet(gui)
% gui = 'view';

global Props

Props.isCaxis = get(Props.view.tag.isCaxis,'value');
if Props.isCaxis, 
    Props.cAx = [ str2num( get(Props.view.tag.cAxLow, 'string') ) ,...
                  str2num( get(Props.view.tag.cAxHigh,'string') ) ] ;
    caxis( Props.cAx )
end
   
if ~Props.isCbar, Props.(gui).cb1 = []; Props.(gui).cb2 = []; 
else
    if Props.isAx1
        Props.(gui).cb1 = colorbar;
    else
        Props.(gui).cb2 = colorbar;
    end
    if ishandle(Props.(gui).cb1)
        cb1pos = get(Props.(gui).cb1,'position'); % set colorbar positions
        set(Props.(gui).cb1,'position',[cb1pos(1:2) 0.01 cb1pos(4)]);
    end
    if ishandle(Props.(gui).cb2)
        cb2pos = get(Props.(gui).cb2,'position');
        set(Props.(gui).cb2,'position',[cb2pos(1:2) 0.01 cb2pos(4)]);
    end
end

% % Design Filter % %
function  [h, r] = DesignFilter(ImageSize)

global Props

xyratio	  = Props.filt.xyratio;
order     = Props.filt.order;
cutoff    = Props.filt.cutoff;
windStr   = get(Props.filt.tag.window,'string');
wind      = windStr{Props.filt.window};
methodStr = get(Props.filt.tag.method,'string');
method      = methodStr{Props.filt.method};
Highpass  = Props.filt.highpass;

%Create desired frequency responce
[f1,f2]     =       freqspace(order, 'meshgrid');
d           =       find(f1.^2 + f2.^2 < cutoff^2);
Hd          =       zeros(order);
Hd(d)       =       1;
if Highpass
    Hd          =       1-Hd;   %Highpass filter
end
%Design filter
switch method
case 'fsamp2',
    h       =       fsamp2(Hd);
case 'fwind1',
    w       =       feval(lower(wind), order);
    h       =       fwind1(Hd,w);
case 'fwind2',
    w       =       feval(lower(wind), order);
    h       =       fwind2(Hd, w*w');
case 'ftrans2',
    F       =       [0 cutoff-0.05 cutoff+0.05 1];
    M       =       [1 1 0 0];
    if Highpass
        M   =       1-M;   %Highpass filter
    end
    h       =       ftrans2(remez(order-1,F,M));
end

%Compute frequency responce
a = 1:1:order;
a = ones(order,1) * a;

% [r, f1, f2] =       freqz2(h, ImageSize(1,1), ImageSize(1,2));

if xyratio < 1,
    [r, f1, f2] =       freqz2(h, ImageSize(1,1), fix(ImageSize(1,2) / xyratio));
    r           =       r((size(r, 1) - ImageSize(1,2)) / 2 + 1: (size(r, 1) - ImageSize(1,2)) / 2 + ImageSize(1,2), :);
else
    [r, f1, f2] =       freqz2(h, fix(ImageSize(1,1) * xyratio), ImageSize(1,2));
    r           =       r(:, (size(r, 2) - ImageSize(1,1)) / 2 + 1: (size(r, 2) - ImageSize(1,1)) / 2 + ImageSize(1,1));
end


% performs cyclic shift of the 1D or 2D matrix %
%  Dx amount of shift
% function [y] = shift(x, Dx);
function shift

global Dat Props

imageCenter      = [(Dat.proc.y(end) - Dat.proc.y(1)) /2,...    % [Y, X]
                    (Dat.proc.x(end) - Dat.proc.x(1)) /2   ] + ...
                   [Dat.proc.y(1) Dat.proc.x(1) ] ;             % offset by roi       

if get(Props.view.tag.isComShift,'value')
	if isfield(Props,'fft')
		scCOM        = [Props.fft.yCom, Props.fft.xCom]; % location of the spatial carrier center of mass [Y, X]
		Dx           = round( imageCenter - scCOM );
	else
		set(Props.view.tag.cmd,'string','Need fft for shift, jerky...')
		return
	end
else
	if isfield(Props,'coords')
		Dx			 = round(imageCenter - [Props.coords.y, Props.coords.x]);
	else
		set(Props.view.tag.cmd,'string','Need coords for shift, first...')
		return
	end
end
x = Dat.proc.image;                                              % convert to slava vars

arraySize   =   ndims(x);
shiftSize   =   ndims(Dx);

%Detect a vector
if length(x) == prod(size(x))
    arraySize = 1;
end
y           =       x;

switch arraySize
case    1,
    SizeX  = prod(size(x));
    switch sign(Dx)
        case    1,
            y(1+Dx:SizeX) =   x(1:SizeX-Dx);
            y(1:Dx)     =   x(SizeX-Dx+1:SizeX);
        case    0,
            ;
        case   -1,
            y(1:SizeX+Dx)     =   x(1-Dx:SizeX);
            y(SizeX+Dx+1:SizeX)       =     x(1:-Dx);
        end
case    2,
    [SizeY, SizeX]  = size(x);
    % Shift the rows
    switch sign(Dx(1))
        case    1,
            y(1+Dx(1):SizeY, :) =   x(1:SizeY-Dx(1), :);
            y(1:Dx(1), :)     =   x(SizeY-Dx(1)+1:SizeY, :);
        case    0,
            ;
        case   -1,
            y(1:SizeY+Dx(1), :)     =   x(1-Dx(1):SizeY, :);
            y(SizeY+Dx(1)+1:SizeY, :)       =     x(1:-Dx(1), :);
        end
    % shift the columns
    x       =       y;
    switch sign(Dx(2))
        case    1,
            y(:, 1+Dx(2):SizeX) =   x(:, 1:SizeX-Dx(2));
            y(:, 1:Dx(2))     =   x(:, SizeX-Dx(2)+1:SizeX);
        case    0,
            ;
        case   -1,
            y(:, 1:SizeX+Dx(2))     =   x(:, 1-Dx(2):SizeX);
            y(:, SizeX+Dx(2)+1:SizeX)       =     x(:, 1:-Dx(2));
        end
end

Dat.proc.image = y; % new
return


function computeFFT    % simple version ... may make a new variable Dat.fft.image in future ...

global Dat

Dat.proc.image = fftshift(fft2(fftshift(Dat.proc.image))); 

dx				= diff(Dat.orig.x);
dx				= dx(1);
dy				= diff(Dat.orig.y);
dy				= dy(1);

% Dat.proc.x		= linspace(-(1./dx)/2,(1./dx)/2,length(Dat.orig.x));
% Dat.proc.y		= linspace(-(1./dy)/2,(1./dy)/2,length(Dat.orig.y));

function computeIFFT % not yet used

global Dat Props

Dat.ift = Dat.proc;

if Props.view.tag.isFilter
    datFoo = Dat.proc.image .* abs(Props.filt.H);
else
    datFoo = Dat.proc.image;
end

 Dat.ift.image = fftshift(ifft2(fftshift(datFoo)));
 
function background_subtract 
%- subtracts background and references as determined by iViewDrkGui -%
global Dat Props
 
if Props.isDrk
if Props.drk.isPrpl
	im	= prpl(Dat.orig.imageRaw); 
else
	im	= Dat.orig.imageRaw;
end

if Props.drk.isMDrk
	 if Props.drk.isPrpl
		 drk	= prpl(Props.drk.drk); 
	 else
		 drk	= Props.drk.drk;
	 end
	im	= im - drk;
end

if Props.drk.isMRef1
	if Props.drk.isPrpl
		ref1	= prpl(Props.drk.ref1); 
	else
		ref1= Props.drk.ref1;
	end
	if Props.drk.isMDrk
		ref1	= ref1	- drk; end	
	im	= im	- ref1;
end

if Props.drk.isMRef2
	if Props.drk.isPrpl
		ref2	= prpl(Props.drk.ref2); 
	else
		ref2= Props.drk.ref2;
	end
	if Props.drk.isMDrk
		ref2	= ref2	- drk; end		
	im	= im	- ref2;
end
else % no drk subtract, check off
im	= Dat.orig.imageRaw;    

end
Dat.orig.image	= im;
 
function imOut = prpl(im)
 
global Dat Props
 
if ~isfield(Props.drk,'drk'), return, end
drk			= Props.drk.drk;	% get variables
thresh			= Props.drk.thresh;
 
 % brute force 1D interpolation - can make 2D is desired %
inds1D			= find(drk >= thresh);
imTmp			= im;
imTmp(inds1D)	= NaN;
 
 for ii	= 1:length(inds1D)
	iInd	= inds1D(ii);
	% wraps around image (top to bottom) like this 
	% - not good with data at edges %
	
	% find first good pixel before bad %
	ind1	= iInd;
	pix1	= imTmp(ind1);
	while isnan(pix1)
		ind1 	= ind1 - 1;
		pix1	= imTmp(ind1);
	end
	
	% find first good pixel after bad %
	ind2	= iInd;
	pix2	= imTmp(ind2);
	while isnan(pix2)
		ind2 	= ind2 + 1;
		pix2	= imTmp(ind2);
	end
	
	imTmp(iInd) = mean([pix1, pix2]);
end

imOut			= imTmp;


%% Add cross section plots to the image, created by imagesc
% Left click - activate/deactivate trace mode

% Author R. Rokitski, 2006
function    [] = xsection(option)
if ~exist('option'),
    if isempty(get(gcf,'WindowButtonUpFcn')),
        option = 'initialize';
    else
        option = '';
    end
end

switch option,
    case 'initialize',
        % save current image data
        % rearrange the figure to add subplots
        % plot xsections and xhairs

        % Identify active axes with image data
        H.figure  =   gcf;
        H.children      =   get(H.figure, 'Children');
        switch length(H.children),
            case    0,
                error('Figure is empty'); return
            case    1,
                if ~strfind(get(get(H.figure, 'Children'), 'Type'), 'axes'),
                    error('No valid axes in current figure'); return
                end
                H.axes      =   gca;
            otherwise
                axesInd =   find(strcmp(get(H.children, 'Type'), 'axes') == 1);
                if length(axesInd) ~= 1,
                    error('Figure should have only one subplot'); return
                end
                H.axes      =  H.children(axesInd);
        end

        H.lines  =   get(H.axes, 'Children');
        if length(H.lines) ~= 1,
            error('Axes should have only one object');return
        end

        H.image  =   get(H.lines);
        if ~sum(strcmp(fieldnames(H.image), 'CData')),
            error('No image in current axes'); return
        end

        %% Create GUI in the figure
        x = H.image.CData;
        figure(H.figure)
        AXim = subplot(4,4,[1 2 3 5 6 7 9 10 11]);
        imagesc(x); axis xy; %colorbar;

        AXvp = subplot(4,4,[4 8 12]);
        plot(x(:, fix(size(x, 2)/2)), 1:1:size(x, 1))
        set(AXvp, 'YLim', [1 size(x, 1)])
        set(AXvp, 'XAxisLocation', 'top')
        set(AXvp, 'YAxisLocation', 'right')
        AXhp = subplot(4,4,[13 14 15]);
        plot(x(fix(size(x, 1)/2),:))
        set(AXhp, 'XLim', [1 size(x, 2)])
        
        AXvpTI =   get(AXvp, 'Position');
        AXhpTI =   get(AXhp, 'Position');

        Font.FontName   = 'Helvetica';
        Font.FontUnits  = 'Pixels';
        Font.FontSize   = 18;
        Font.FontWeight = 'bold';
        Font.FontAngle  = 'normal';
        
        H.xvalue      = uicontrol(H.figure,'Style','edit','Units','Normalized', Font, ...
                                  'Position',[AXvpTI(1) AXhpTI(1)+AXhpTI(1)-0.05  AXvpTI(3), 0.05],...
                                  'BackGroundColor',[ 0 0 1],'ForeGroundColor',[ 1 1 1],...
                                  'Tag','XVALUE',...
                                  'TooltipString','X value (Read Only)',...
                                  'String',num2str(fix(size(x, 2)/2)));
        H.yvalue      = uicontrol(H.figure,'Style','edit','Units','Normalized', Font, ...
                                  'Position',[AXvpTI(1) AXhpTI(1)+AXhpTI(1)-0.05 - 0.075  AXvpTI(3), 0.05],...
                                  'BackGroundColor',[ 0 0 1],'ForeGroundColor',[ 1 1 1],...
                                  'Tag','YVALUE',...
                                  'TooltipString','Y value (Read Only)',...
                                  'String', num2str(fix(size(x, 1)/2)));
                          
        set(H.figure,'WindowButtonUpFcn','xsection(''toggle'');');
        set(H.figure, 'Toolbar', 'Figure');
    case 'trace',
        axesHandles =   get(gcf, 'Children');
        controlHandles =  axesHandles(find(strcmp(get(axesHandles, 'Type'), 'uicontrol') == 1));
        axesHandles = axesHandles(find(strcmp(get(axesHandles, 'Type'), 'axes') == 1));
        axisHandle  =   axesHandles(end);

        x_rng = get(axisHandle,'Xlim');
        y_rng = get(axisHandle,'Ylim');
        temp = get(axisHandle,'currentpoint');
        set(controlHandles(1), 'String', num2str(fix(temp(1,1))));
        set(controlHandles(2), 'String', num2str(fix(temp(1,2))));

        if temp(1,1) >= x_rng(1) && temp(1,1) <= x_rng(2) && ...
                temp(1,2) >= y_rng(1) && temp(1,2) <= y_rng(2),
            AllLines    =   get(axisHandle, 'Children');
            delete(AllLines(find(strcmp(get(AllLines, 'Type'), 'line')==1)));
            x = get(get(axisHandle, 'Children'), 'CData');
            line(x_rng, temp(:,2));
            line(temp(:,1), y_rng);

           % plot(axesHandles(2), x(:, fix(temp(1,1))), 1:1:size(x, 1))
           % slavas line
		   axes(axesHandles(2)) 
		   plot(x(:, fix(temp(1,1))), 1:1:size(x, 1))
            
			set(axesHandles(2), 'YLim', y_rng);
            set(axesHandles(2), 'XAxisLocation', 'top')
            set(axesHandles(2), 'YAxisLocation', 'right')
          
			axes(axesHandles(1)) 
			plot( x(fix(temp(1,2)), :))
            %% plot(axesHandles(1), x(fix(temp(1,2)), :)) 
			% slavas line
			
            set(axesHandles(1), 'XLim', x_rng);

        end

    case 'toggle',
        if isempty(get(gcf,'WindowButtonMotionFcn')),
            set(gcf,'WindowButtonMotionFcn','xsection(''trace'');');
        else
            set(gcf,'WindowButtonMotionFcn','');            
        end
    otherwise
end
