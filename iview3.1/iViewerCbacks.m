% Callbacks for iView v2.1
%
% GUI for image visualization and manipulation- particularly suited
% to gaussian beam profiling, scectroscopy, etc.
%
% Basic framework here- could be used as base to add many more
% file import/export, comparison+multiple images, multiple peaks
% processing/analysis features (image processing toolbox features, etc.)
%
% I'll call the length and wavelength mm and nm, respectively, thought 
% they don't necessarily have to be so ...
%
% User .M files: WINSPECREAD, iViewInitProps (sets GUI data initially)
%                Colormaps: CMAPRAINBOW2, JET2, JET3
% User MEX, C Files, and other: UIGETFILES.DLL, WHYUNO.DLL, Calibration_Lines.pdf
% User .MAT files: none 
% Subfunctions: A few
% See also: iViewer1280x1024.fig, iViiewerXP.fig
%
% * Note: Requires Image Processing Toolbox to be installed *


% Author: Kevin Tetz (as inspired by Jessica Yuen)
% Last revision: 06-Apr-2004
% Revision History:
% 21-Apr-2005: v2.75: updated drk subtract (fixed bugs)
%                     Added units to FFT 
% 06-Apr-2004: v1.2:  Added movie export, made Props and Dat globals
%              - made separate .m file for subfunctions, to be called from elsewhere               
%              - made Dat a single 1-D structure array 
% 06-Jan-2004: Finished zooming/aspect ratios w/ units
% 23-Nov-2003: change x- and ydata to plot instead of just labels
% 11-Nov-2003: Changed max value (was for COM slice)
%              Support multiple images
%              Figures in new window
%              COM/FWHM butons combined into 'STATS'
%              Used IMCROP for ROI
% 07-Nov-2003: added export to workspace, support for 'WinSpec Read files
% 10-Jun-2002: Created

% To add: - support for multiple peaks (finding, COM, etc.)
%         - ellips eccentricity , etc.
%         - zoom memory (zoom out/prev)
%         - export images
%         - smoothing (data)/image contrast, etc.
%           add the same plot in ax2 if change the image
%         - check to make sure contour, and 3D plot not crash b.c. too big
%         - hist axes scale to max/min w/ new plot- could have memory
%         - maybe save the FFT and delete it... right now recomputes...
%         - sticthing of multiple images, simple addition/subtration, etc.
%         of images
%         - more easily compare properties (profile, especially) btw
%         images... display in fig.. say: select: image1, profile 1&2, etc.,
%          (have a profile... open a dialog to select the profile for
%          arbitrary images ... can easily use the same functions already here)

% Known Bugs: - 3-D data aspect ratios w/ different units... z sometimes
%               curious with mixed nm/mm (also code somewhat redundant)
%             - Bkg and Scaling need tweaking (image specific??)
%               get slightly off FWHM when scaled ...                   
%             - some axes tick labels obscured on thin plots
%             - "uigetfiles.dll" seems to only be able to "see" 250 files/dir


function varargout = iViewerCbacks(cback) %% input callback switch

global F Props Dat Dats                     %  globals to be used by all
                                            % F is struct of function
                                            % handles 
                                            % Props is struct of relevant
                                            % gui properties
                                            % Dat is current data
                                            % Dats is the read image data

warning off                                 % turns off warnings (log of 0, etc.)

if ~exist('whyUNO.dll') | ~exist('uigetfiles.dll') %% make sure dlls are found
    set(cmdHndl,'string','Error: Need whyUNO.dll & uigetfiles.dll in matlab path')
end

% well worn handles %
cmdHndl      = Props.view.tag.cmd;

% check a few values %
Props.isMm   = get(Props.view.tag.isMm,'value');
Props.isNm   = get(Props.view.tag.isNm,'value');
Props.isCbar = get(Props.view.tag.isCbar,'value');
datNames     = get(Props.view.tag.datNames,'String'); 

% get colormap %
cmaps        = get(Props.view.tag.cmaps,'string');
Props.cmap   = cmaps{get(Props.view.tag.cmaps,'value')};

%% error checking for is changed units... needs to be orig, open, next,
%% prev, datNames, hist, cmap, or the like
%% or better, error for roi, coords, inten, 

if ~strcmp(cback,{'orig','open'}) %% can only change units then open or orig
    errmsg = 'Error: have to display original if you change units';
    if Props.isMm %% check to see that each profile has same units
        if Props.isNm
            if ~strcmp(Props.prevUnits.x,'nm')
                set(cmdHndl,'string',errmsg)
                return, end   
        else
            if ~strcmp(Props.prevUnits.x,'mm')
                set(cmdHndl,'string',errmsg)
                return, end   
        end
    else %% mm not checked
        if Props.isNm
            if ~strcmp(Props.prevUnits.x,'nm')
                set(Props.view.tag.cmd.Hndl,'string',errmsg)
                return, end
        else
            if ~strcmp(Props.prevUnits.x,'pixels')
                set(Props.view.tag.cmd.Hndl,'string',errmsg)
                return, end       
        end
    end
end

if ~ishandle(Props.view.ax1)|~ishandle(Props.view.ax2) %% bug check
    set(cmdHndl,'string','Axes got mixed up... restart&kill bug'), return
end

if ~strcmp(cback,'exit') % can exit before opening a file
if ~isstruct(Dat), cback = 'open'; %% open data if none 
	set(cmdHndl,'string','You need some data first, jerky...'), end 
end
if (isfield(Props,'coordsHndl')&ishandle(Props.coordsHndl)), delete(Props.coordsHndl), end

% % Clear variables and plot features if needed %% (change these var names)
axChk    = strcmp(cback,{'profile','coords','cmd','setcmap','newfig','export'});
axChkDat = strcmp(cback,{'profile','coords','linear','log','cmd','export','newfig','setcmap','next','prev','datname'});

% % remove profile lines on ax1
if (isfield(Props,'profHndl')&ishandle(Props.profHndl)&isempty(find(axChk==true)))
    delete(Props.profHndl), Props.profHndl = [];  end

% % remove profile data ... leave empty
if isempty(find(axChkDat==true)),...
    coords.lims = []; coords.y = []; coords.x = [];    
    profile     = struct('dat' ,[],'fwhm',[],'max',[],...
                         'xDat',[],'coords',coords);
    Dat.profile = profile;
end

if Props.isMm %% check that mm range is valid if checked
    if (str2num(get(Props.view.tag.mmX,'string')) <= 0) |...
            (str2num(get(Props.view.tag.mmY,'string')) <= 0) 
        set(cmdHndl,'string','Please enter a number > 0 if using mm...')
        return
    end
end

if Props.isNm %% check that nm range is valid if checked
    if (str2num(get(Props.view.tag.nmPix,'string')) <= 0) |...
       (str2num(get(Props.view.tag.nmCal,'string')) <= 0) |...
       (str2num(get(Props.view.tag.nmDisp,'string')) <= 0) 
       set(cmdHndl,'string','Please enter a number > 0 if using nm...')
       return
   end
end

if strcmp(Props.view.ax2Type,'profile') %% saves coords to use for profile
    DatPrev.profile = Dat.profile;   else DatPrev.profile = []; end
if strcmp(Props.view.ax2Type,'inten'), 
    DatPrev.Irng    = Dat.proc.Irng; else DatPrev.Irng    = []; end

switch cback  % % switchyard for GUI callbacks       
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'open' % % opens new image, clears indicators

feval(F.clearIndicators,'all')      

wd = cd; 
if get(Props.view.tag.isOpenVariable,'value')
	% opens variable from workspace
	isOpenVar	= 1;
	
	vars = evalin('base','whos');
	zzstr  = {vars.name};
	assignin('base','zzstr',zzstr)	
	[varNo,isOk] = evalin('base','listdlg(''PromptString'',''Select a variable to import:'',''ListString'',zzstr)');		
	if isOk	== 0, set(cmdHndl,'string','Couldn''t find anything good??'), return, end
	fname		= zzstr(varNo);
	assignin('base','varNo',varNo)
	pname		= 'workspace';
		
else
	isOpenVar	= 0;
% opens file 
if isfield(Props,'pname') %% give start path as last directory
	 [fname, pname] = uigetfiles('*.hdf;*.pcx;*.xwd;*.ico;*.cur;*.gif;*.tif;*.jpg;*.bmp;*.png;*.dat;*.spe',...
        'Open Image File',Props.pname);   
else
    [fname, pname] = uigetfiles('*.hdf;*.pcx;*.xwd;*.ico;*.cur;*.gif;*.tif;*.jpg;*.bmp;*.png;*.dat;*.spe', 'Open Image File');   
end
end % workspace or file

if pname	== 0, set(cmdHndl,'string','Couldn''t find anything good??'), return, end
Props.pname = pname;

nFiles		= size(fname,2);

filenameChkWin 	= fullfile(pname, fname{1});
if findstr		(lower(filenameChkWin),'.spe') %% for WINSPEC32 data files
     datWinspec = winspecread('filename',filenameChkWin); 
     nDatSets 	= size(datWinspec,3);
     fname 		= repmat(fname,nDatSets,1);   
else
     nDatSets = nFiles;
end

% % reads multiple files, and adds each concurrent to gui via setappdata
iFile = 1;
for iData    = Props.nImages + 1 :Props.nImages + nDatSets 

	filename = fullfile(pname, fname{iFile});

	if isequal(fname,0)|isequal(pname,0)
		disp('File not found'), return, end
if ~isOpenVar,	cd(pname); end
	if findstr(filename,'.dat') %% for ascii data files
		dat    = load(filename);
		iFname = fname{iFile};
	elseif findstr(lower(filename),'.spe') %% for WINSPEC32 data files
		if~exist('winspecread')
			set(cmdHndl,'string','You need winspecread for this ...'), return
		else 
			if ~exist('datWinspec'), set(cmdHndl,'String','Make the .spe file first in the list'), return, end
			if nFiles > 1 & ndims(datWinspec) > 2
				set(cmdHndl,'string','Please multiframe WINSPEC file by intelf ...'), cd(wd), return
			else    
				dat = uint16(datWinspec(:,:,iFile));
			end
		end
		iFname = [fname{iFile}, ', Frame',num2str(iFile)];
	elseif strcmp(filename(1:9),'workspace')
		
		% 	var		= d(wVar).name;
		assignin	('base','iFile',iFile)
		dat			= evalin('base','eval(zzstr{varNo(iFile)})');
		iFname		= fname{iFile};
								
	else
		dat    = (imread(filename)); 
% 		iFname = fname{iFile - 1};
		iFname = fname{iFile};
	end

	newFiles{iFile} = iFname;

	% % change color images to intensity only
	if ndims(dat) >2
		dat = mean(dat,3); 
	end

	% % (re)set user variable data  
	roi          = [1 size(dat,2) 1 size(dat,1)]; %% opt out-> Props.roi in future... set to Dat.proc.roi as wells
	orig.image   = dat; units.x = 'pixels';
    proc.image   = dat; units.y = 'pixels';
       
    [orig.xPixels,orig.x,proc.x] = deal(linspace(1,roi(2),roi(2)));
	[orig.yPixels,orig.y,proc.y] = deal(linspace(1,roi(4),roi(4)));

	[orig.x,orig.y,units]        = feval(F.convertUnits,orig.xPixels,orig.yPixels);
    
	coords  = struct('lims',[], 'x',[]  , 'y',[] );
	profile = struct('dat',[] , 'coords',coords, 'max',[], 'fwhm',[], 'xDat',[]);

	Dat = struct('orig',    orig,      'proc',proc   ,...
                 'profile', profile,   'fname',iFname,...
                 'units',   units);
   if iData == 1;  
		Dats = Dat;
	else
	   Dats(iData) = Dat;
	end

	if iFile == 1, datInit = Dat; end
	iFile = iFile + 1;
end %% load new files
cd(wd)

Props.currImage 	= Props.nImages + 1;
Props.nImages   	= nDatSets + Props.nImages;

Dat             	= Dats(Props.currImage);         %% make first file read current %
Dat.orig.image  	= double(Dat.orig.image);        %% less trouble later- may keep uint8, etc., for memory concerns
Dat.orig.imageRaw	= Dat.orig.image;            %% if needed
set					(Props.view.fig,...
					'Name',['Image= ',fname{1}]) %% figure title
    
% % plot and label image and axes
feval(F.clearIndicators,'orig')                 %% clear indicators and ax2
axes(Props.view.ax2), cla, set(gca,'visible','off'), title('')
axes(Props.view.ax1)

Props.view.ax2Type = 'none';
Props.roiInds      = [1 size(datInit.orig.image,2) 1 size(datInit.orig.image,1)];
Props.pixSize      = size(datInit.orig.image);
Props.view.ax1     = gca;

feval(F.localPlot,'orig');

newList 			= cat(1,datNames,newFiles(:));  

set(Props.view.tag.datNames,'string',newList,...
	                        'value',Props.currImage + 2)
set(Props.view.tag.lin,'Value',1), 
set(Props.view.tag.log,'Value',0) 

feval			(F.axset)  %% set current axes 
Props.view.ax1 	= gca;
dispWhy			(cmdHndl)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'orig'  % % (re)plots original image and makes current

feval(F.clearIndicators,'all') %% clear indicators and ax2

if ishandle(Props.view.cb2), set([gca,Props.view.cb2],'Visible','off'), end

Props.view.ax2Type        = 'none';
dat                       = Dat.orig; 

if get(Props.view.tag.isDrk,'value')
	image	= dat.image;
else
	image	= dat.imageRaw;
end

orig                      = struct('image'  	,image  	  , ...
								   'imageRaw'	,dat.imageRaw 	, ...
                                   'xPixels'	,dat.xPixels, ...
                                   'yPixels'	,dat.yPixels  ...
                                  );                                     %% only keep pixels
Props.roiInds             = [1 size(orig.image,2) 1 size(orig.image,1)]; %% ROI in indices (pixels)
[orig.x,orig.y,Dat.units] = feval(F.convertUnits,...
                                  Dat.orig.xPixels,Dat.orig.yPixels...
                                  );
Dat.orig                  = orig;
Dat.orig.image            = double(Dat.orig.image);
Dat.proc                  = orig;                                        %% overwrite old data
Props.view.ax2Type        = 'none'; 
Props.view.isFFT          = 0;

feval					(F.background_subtract);
feval					(F.localPlot,'orig');                                               %% replot new
dispWhy					(cmdHndl)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'origRoi'  % % (re)plots original image and makes current

if ~Props.isAx1, axes(Props.view.ax2), cla, set(gca,'visible','off'), title(''),
    if ishandle(Props.view.cb2), set([gca,Props.view.cb2],'Visible','off'), end,  end

Props.view.isFFT   = 0;
Props.view.ax2Type = 'none';    
Dat.proc           = Dat.orig;     
                                                            %% Set new range
Dat.proc.x         = Dat.proc.x(Props.roiInds(1):Props.roiInds(2));
Dat.proc.y         = Dat.proc.y(Props.roiInds(3):Props.roiInds(4));
Dat.proc.image     = Dat.proc.image(...                    
                                    Props.roiInds(3):Props.roiInds(4),...
                                    Props.roiInds(1):Props.roiInds(2)...
                                    ) ;
Props.isAx1        = 1;
axes(Props.view.ax1), cla
feval(F.localPlot,'roi');                                   %% subfun to plot
dispWhy(cmdHndl)                                    

 % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'roi' % % ROI button

feval(F.clearIndicators,'orig') % % set panel and get data
	
roiIndsPrev	= Props.roiInds;

% % determine manual or line entry, then check for valid ranges    
if ~get(Props.view.tag.manZoom,'value')
    set(cmdHndl,'string','Tip: use right mouse button or hold shift for square...')
%     axes(Props.view.ax1) %% needed this now ... why not before? could use ax2 ...
    roi = getrect;
    roi = [roi(1) roi(1)+roi(3) roi(2) roi(2)+roi(4)]; %% use this notation to be consistent with old code
                                                        % and avoid other changes... can remove if everything                                                        % changes... can
                                                        % else reconciled                                             
else %% manual (prompt) entry

    limits = Props.roi;
        
    if Props.isMm,limUnits = 'mm';
    else, limUnits = 'pixels'; end

    if Props.isNm, limUnits = ['nm; ',limUnits]; end
    
    % % Creates dialog box to prompt user for ROI
    promptStr = {['X-Y extents [',limUnits,']: [xmin xmax ymin ymax]']}; 
    titleStr = 'Image Region of Interest'; nLines = 1; 
    iniStr = {num2str(limits)}; %% init values
    set(0,'defaultfiguretoolbar','none') %% Must turn it off for dialog to work properly
    Roi = inputdlg(promptStr,titleStr,nLines,iniStr); %% Sets range equal to input dialog that allows user to type in range
    set(0,'defaultfiguretoolbar','figure') % turn def back on
    if isempty(Roi), return, end
    roi = str2num(char(Roi)); 
    roi(1:2) = sort(roi(1:2)); roi(3:4) = sort(roi(3:4));  %% make sure right order
    
end %% manual or line ... ROI set above, both plot below

% % maybe add a new error check on the manual entries % %
% if xPix(1) < 1 | yPix(2) > size(Dat.orig.image,2) | xPix(2) < 1 | yPix(2) > size(Dat.orig.image,1)
%     set(cmdHndl,'String','Your box isn''t in the image, jerky... '), return, end

% % need to select sub portion of image ... index with x y data % %
xInds = [find(Dat.proc.x == feval(F.nearest,Dat.proc.x,roi(1))),find(Dat.proc.x == feval(F.nearest,Dat.proc.x,roi(2)))];
yInds = [find(Dat.proc.y == feval(F.nearest,Dat.proc.y,roi(3))),find(Dat.proc.y == feval(F.nearest,Dat.proc.y,roi(4)))];
         
Dat.proc.x = Dat.proc.x(xInds(1):xInds(2));
Dat.proc.y = Dat.proc.y(yInds(1):yInds(2));
Dat.proc.image = Dat.proc.image(yInds(1):yInds(2),xInds(1):xInds(2)); %% Set new range

if ~Props.view.isFFT
    Props.roi = roi; %% always in CURRENT UNITS
    Props.roiInds	= [xInds,yInds];
% 	Props.roiInds	= Props.roiInds + roiIndsPrev - 1;	% adds previoius values 
end

% % use these for aspect ratio scaling only % %
Dat.proc.xPixels = Dat.proc.x; Dat.proc.yPixels = Dat.proc.y;
if Props.isMm
    Dat.proc.xPixels = feval(F.mm2pix,'x',Dat.proc.x,Props.pixSize);
    Dat.proc.yPixels = feval(F.mm2pix,'x',Dat.proc.x,Props.pixSize);
end
if Props.isNm, Dat.proc.xPixels = feval(F.nm2pix,Dat.proc.x); end

feval(F.localPlot,'roi'); %% subfun to plot

dispWhy(cmdHndl)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'hist' % % histogram

feval(F.localPlot,'hist'); %% subfun
dispWhy(cmdHndl)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'ft' % % Fourier Transform button

if ~Props.view.isFFT
    Props.view.isFFT = 1;
    feval(F.computeFFT)
end

feval(F.localPlot,'ft'); %% subfun
dispWhy(cmdHndl)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'inten' % % intensity button 
   
dat = double(Dat.proc.image);
lowIrng = Props.lowIrng; highIrng = Props.highIrng;

if ~get(Props.view.tag.isLimRng,'value') & ...
   ~get(Props.view.tag.isLimRng2,'value') %% auto scale off
   % % prompt for caxis range
   promptstr = {'Enter the range of intensity values: '}; 
   titlestr = 'Range of Intensity Values to be plotted'; nlines = 1; 
   limits = num2str(([min(dat(:)), max(dat(:))])); % Sets the initial range limits
   inistr = {limits};
   set(0,'defaultfiguretoolbar','none') % Must turn it off for dialog to work properly
   Range = inputdlg(promptstr,titlestr,nlines,inistr); 
   set(0,'defaultfiguretoolbar','figure') %% turn def back on
   if isempty(Range), return, end
   range = str2num(char(Range)); 
   dispWhy(cmdHndl)
elseif ~get(Props.view.tag.isLimRng,'value') & ...
        get(Props.view.tag.isLimRng2,'value') %% auto 5%-95%, or whatever the default range is ... 256 bins
      
    %% this uses a percentage of the maximum (comment to use the above)
    range = [lowIrng*max(dat(:)),highIrng*max(dat(:))]; %% or this
    dispWhy(cmdHndl)
elseif ~get(Props.view.tag.isLimRng2,'value') & ...
        get(Props.view.tag.isLimRng,'value')
        %% this uses a percentage of the pixels
        [counts,x] = hist(dat(:),256);
        sumCounts = sum(counts);
        numCounts = 0; ii = 1;
        while numCounts < lowIrng * sumCounts
            numCounts = numCounts + counts(ii);
            ii = ii+1;
        end
        range(1) = x(ii-1);
        numCounts = 0; ii = length(counts);
        while numCounts < (1-highIrng) * sumCounts
            numCounts = numCounts + counts(ii);
            ii = ii - 1;
        end
        range(2) = x(ii+1);
        if range(2) == range(1), range(2) = x(ii+2); end
        dispWhy(cmdHndl)
else %% error box boxes checked (could make a callback to turn off if get around to it)
    set(cmdHndl,'string','Have to choose one or the other...'), return
         
end
Dat.proc.Irng = range;

feval(F.localPlot,'inten');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'contour' % % contour button

feval(F.localPlot,'contour');

dispWhy(cmdHndl)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'plot3' % % 3-D plot button

feval(F.localPlot,'plot3');

dispWhy(cmdHndl)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'stats' % % stats button ... might be nice to consilidate varible names a touch...

datOrig = double(Dat.orig.image); %% get appropriate data
datProc = abs(double(Dat.proc.image));
  
% % Computes centroid -- with appropriate units % % 
[Dat.orig.xCom,Dat.orig.yCom] = feval(F.centroid,Dat.orig.x,Dat.orig.y,datOrig);
[Dat.proc.xCom,Dat.proc.yCom] = feval(F.centroid,Dat.proc.x,Dat.proc.y,datProc);
 
% % Finds stats % %
maxOrig = (max(datOrig(:))); maxProc = (max(datProc(:)));
minOrig = (min(datOrig(:))); minProc = (min(datProc(:)));
sumOrig = sum(datOrig(:));   sumProc = sum(datProc(:));
avgOrig = sumOrig/length(datOrig(:)); 
avgProc = sumProc/length(datProc(:));

% % gets column(row) to define slice % %
ixOrig = find(Dat.orig.y == feval(F.nearest,Dat.orig.y,Dat.orig.yCom));
jy0rig = find(Dat.orig.x == feval(F.nearest,Dat.orig.x,Dat.orig.xCom));
ixProc = find(Dat.proc.y == feval(F.nearest,Dat.proc.y,Dat.proc.yCom));
jyProc = find(Dat.proc.x == feval(F.nearest,Dat.proc.x,Dat.proc.xCom));

xOrig = datOrig(ixOrig,:); yOrig = datOrig(:,jy0rig);
xProc = datProc(ixProc,:); yProc = datProc(:,jyProc);

%% call to sub, returns FWHM, max for slice @ COM
[xFw,xM,xd1,xd2]     = feval(F.localFWHM,xOrig,Dat.orig.x);
[yFw,yM,yd1,yd2]     = feval(F.localFWHM,yOrig,Dat.orig.y); 
[xPFw,xPM,xPd1,xPd2] = feval(F.localFWHM,xProc,Dat.proc.x);
[yPFw,yPM,yPd1,yPd2] = feval(F.localFWHM,yProc,Dat.proc.y);

% % some funny variable names, i guess ... what are each of these? % %
xfwx = [xd1;xd2] ; xfwy = Dat.orig.yCom*[1;1] ;
yfwy = [yd1;yd2] ; yfwx = Dat.orig.xCom*[1;1] ;

xPfwx = [xPd1;xPd2] ; xPfwy = Dat.proc.yCom*[1;1];
yPfwy = [yPd1;yPd2] ; yPfwx = Dat.proc.xCom*[1;1] ;

% % Set data and indicators to appropriate values
Dat.orig.max   = maxOrig;  Dat.proc.max   = maxProc;
Dat.orig.avg   = avgOrig;  Dat.proc.avg   = avgProc;
Dat.orig.sum   = sumOrig;  Dat.proc.sum   = sumProc; 
Dat.orig.xFwhm = xFw;      Dat.orig.yFwhm = yFw;
Dat.proc.xFwhm = xPFw;     Dat.proc.yFwhm = yPFw;

% mark processed data (see comments below for orig) %
hold on
pC1 = line(xPfwx,xPfwy);
pC2 = line(yPfwx,yPfwy);
p4  = plot(Dat.proc.xCom,Dat.proc.yCom,'ko','markerfacecolor',[0 1 0]);
p3  = plot(Dat.proc.xCom ,Dat.proc.yCom,'kx');

set([p3,p4],'Markersize',6)
set([pC1,pC2],'MarkerFaceColor',[0 1 0],'markersize',4,...
              'Marker','o','MarkerEdgeColor',[0 0 0],'Color',[0 0 0]);
hold off

% % marks FWHM with lines % % 
if ~Props.view.isFFT %
    if ~Props.isAx1, axes(Props.view.ax1), end
    
    hold on
    p1 = line(xfwx,xfwy); p2 = line(yfwx,yfwy);

    % % set colors appropriately % %
    set([p1,p2],'MarkerFaceColor',[1 0 1],'markersize',4,...
                'Marker','o','MarkerEdgeColor',[0 0 0],'Color',[0 0 0]);

     % % Marks the center of mass on the figure with an x
     p2 = plot(Dat.orig.xCom,Dat.orig.yCom,'ko','markerfacecolor',[1 0 1]);
     p1 = plot(Dat.orig.xCom ,Dat.orig.yCom,'kx');
     hold off
     set([p1,p2],'Markersize',6)
end

if ~Props.isAx1, axes(Props.view.ax2), end %% want figure to stay current

if strcmp(Props.view.ax2Type,'ft')
    Props.fft.xCom = Dat.proc.xCom;
    Props.fft.yCom = Dat.proc.yCom; end

% % set ALL indicators % %
set(Props.view.tag.xCOM,      'string', num2str(Dat.orig.xCom,  6)  )    
set(Props.view.tag.yCOM,      'string', num2str(Dat.orig.yCom,  6)  )   
set(Props.view.tag.xCOMProc,  'string', num2str(Dat.proc.xCom,  6)  )    
set(Props.view.tag.yCOMProc,  'string', num2str(Dat.proc.yCom,  6)) 
set(Props.view.tag.yFWHM,     'string', num2str(Dat.orig.yFwhm, 5))    
set(Props.view.tag.xFWHM,     'string', num2str(Dat.orig.xFwhm, 5))    
set(Props.view.tag.xFWHMProc, 'string', num2str(Dat.proc.xFwhm, 5))    
set(Props.view.tag.yFWHMProc, 'string', num2str(Dat.proc.yFwhm, 5))    
set(Props.view.tag.max,       'string', num2str(maxOrig,        5) )  
set(Props.view.tag.maxProc,   'string', num2str(maxProc,        5))
set(Props.view.tag.sumOrig,   'string', num2str(sumOrig,        '%1.3g'))
set(Props.view.tag.sumProc,   'string', num2str(sumProc,        '%1.3g'))
set(Props.view.tag.avgOrig,   'string', num2str(avgOrig,        5))
set(Props.view.tag.avgProc,   'string', num2str(avgProc,        5))
set(Props.view.tag.minOrig,   'string', num2str(minOrig,        5))
set(Props.view.tag.minProc,   'string', num2str(minProc,        5))

dispWhy(cmdHndl)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
case 'coords' % % get coords button
  
dat = double(Dat.proc.image);

[xCoord, yCoord] = ginput(1);               %% allows user to pinpoint point
xCoord           = feval(F.nearest,Dat.proc.x,xCoord);
yCoord           = feval(F.nearest,Dat.proc.y,yCoord);
xInd             = [find(Dat.proc.x == xCoord)];
yInd             = [find(Dat.proc.y == yCoord)];   
datVal           = dat(yInd, xInd);

hold on                                                   %% plot the point 
c1               = plot(xCoord,yCoord,'ko','markerFaceColor',[1 0.5 0]); 
c2               = plot(xCoord,yCoord,'k+'); hold off
set([c1, c2], 'markersize', 6)

% % Set Indicators % %
set(Props.view.tag.xCoord, 'String', num2str(xCoord));
set(Props.view.tag.yCoord, 'String', num2str(yCoord));
set(Props.view.tag.ICoord, 'string', num2str(datVal))

% % Set data values % %
Props.coords.x      = xCoord;  Props.coords.y     = yCoord;
Props.coords.I      = datVal;  Props.coordsHndl = [c1; c2];    

dispWhy(cmdHndl)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'profile' % % Image Profile button
      
 % % check to make sure the ROI is correctly handled    
if ~Props.isAx1
    set(cmdHndl,'String','Please profile the left axis only ...'), return, end

if get(Props.view.tag.profXsec,'value')
% opens new window and profiles using Slavas two dimensional xsection
	figure
	imagesc(Dat.proc.image)
	feval(F.xsection)
	return
end

axes(Props.view.ax1)     
[Xa,Ya,A,state] = getimage;

if ~state, set(cmdHndl,'string','Can''t... need image (hit orig)')
    return, end

if get(Props.view.tag.profMan,'value') %% mouse input or dialog

% % Creates dialog box to prompt user for profile coords
    if isfield(Dat.profile,'coords') & ~isempty(Dat.profile.coords.lims)
        limitsInit = Dat.profile.coords.lims(end,:); %% Sets the initial values to the current axis range
    else limitsInit = Props.roi; end
    
    if Props.isMm, limUnits = 'mm';
    else limUnits = 'pixels'; end
     
    if Props.isNm, limUnits = ['[nm; ',limUnits,']']; end
		  
	promptStr = {['Profile coordinates, in ',limUnits,': [x1 x2 y1 y2]']}; 
    titleStr  = 'Image Profile'; nLines = 1; 
    iniStr    = {num2str(limitsInit)};
    set(0,'defaultfiguretoolbar','none')                    %% Must turn it off for dialog to work properly
    Lim       = inputdlg(promptStr,titleStr,nLines,iniStr); %% Sets range equal to input dialog that allows user to type in range
    set(0,'defaultfiguretoolbar','figure')                  %% set def
    if isempty(Lim), return, end
    lim       = str2num(char(Lim)); 
    
    if Props.isNm
        if lim(3)~=lim(4), lim(4)=lim(3); 
            set(cmdHndl,'string','Note: wavelength profile sets ''y'' to horizontal')
        end
    end
    limPlot(1:2) = lim(1:2); limPlot(3:4) = lim(3:4);
	
	imDatProf		= Props.view.ax1Type;
	if strcmp(imDatProf,'ift')
	% select the data to profie - ROI, FFT or PROC, etc.
		imToProfile	= Dat.ift.image;
	else
		imToProfile	= Dat.proc.image;
	end	
	
	
%     [xCoords, yCoords, imProf, xLim, yLim] = improfile(Dat.proc.x,Dat.proc.y,Dat.(Props.view.ax1Type).image,limPlot(1:2),limPlot(3:4),...
%                                                        Props.nProfPoints,'bilinear'); 
    [xCoords, yCoords, imProf, xLim, yLim] = improfile(Dat.proc.x,Dat.proc.y,imToProfile,limPlot(1:2),limPlot(3:4),...
                                                       Props.nProfPoints,'bilinear'); 
    if size(xLim(:),1) ~= 2
        set(cmdHndl, 'string', 'Please use only two points in each profile...'), return
    end
else %% user mouse input  
     [xCoords,yCoords,imProf,xLim,yLim] = improfile(Props.nProfPoints,'bilinear'); 
    if size(xLim(:),1) ~= 2
        set(cmdHndl,'string','Please use only two points in each profile...'), return
    end
    if Props.isNm
        if yLim(1)~=yLim(2), yLim(2) = yLim(1); 
            [xCoords,yCoords,imProf,xLim,yLim] = improfile(Dat.proc.x,Dat.proc.y,Dat.(Props.view.ax1Type).image,xLim,yLim,...
                Props.nProfPoints,'bilinear'); 
            set(cmdHndl,'string','Note: wavelength profile sets ''y'' to horizontal')
        else
            dispWhy(cmdHndl)
        end
    end
end

% % add profile and plot all % %
Dat.profile = feval(F.catProfile,Dat.profile,xCoords,yCoords,imProf,xLim,yLim);
feval(F.localPlot,'profile')

% % set indicators % %
set(Props.view.tag.profMax,'string',...
    num2str(Dat.profile.max(end)))    
set(Props.view.tag.profFWHM,'string',...
    num2str(Dat.profile.fwhm(end)))

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'setcmap' % % cmap pull down menu

if ~exist(Props.cmap) 
    if any(strcmp(Props.cmap,{'cmaprainbow2','jet2','jet3'}))
        set(cmdHndl,'string','Ask Kevin for this map if you want it...')
    else, return
    end
else
    colormap(Props.cmap);
    set(cmdHndl,'string',whyUNO)
end
        

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'linear' % % linear radio button

axType = 'linear';

set(Props.view.tag.log,'Value',0); %% radiobuttons
set(Props.view.tag.lin,'Value',1)

if Props.isAx1  %% cases for left and right plots (left only imagesc)
    axes(Props.view.ax1), tit = get(get(gca,'Title'),'String');
    logInd = findstr(tit,'log');
    if strcmp(Props.view.ax2Type,'profile')
       axes  (Props.view.ax2),  plot(Dat.profile.dat)
       ylabel('Intensity'),     title('Image Profile')
       xlabel('Position in profile')
       feval(F.cbarsSet,'view')
% % keep the figure current so you can slice multiple times
       axes(Props.view.ax1) 
       
    end
	if strcmp(Props.view.ax2Type,'hist')
		set(Props.view.ax2,'yscale','linear'), end
    if strcmp(Props.view.ax1Type,'ift')
        if get(Props.view.tag.isPhase,'value')
            datFoo = angle(Dat.ift.image);
        else
            datFoo = abs  (Dat.ift.image);
        end
    else  
        datFoo = Dat.proc.image;
    end
    imagesc(Dat.proc.x,Dat.proc.y,datFoo)
    feval  (F.aspectRatios)
    feval  (F.drawProfLines,Dat.profile); 
    feval  (F.axset)
    Props.view.ax1 = gca;
    
else %% all other types of images
    axes(Props.view.ax2), tit = get(get(gca,'Title'),'String'); logInd = findstr(tit,'log');
    if strcmp(Props.view.ax2Type,'plot3')
        dat3D = double(Dat.proc.image);
        set(gca,'ZScale',axType)
        xmax = max(size(dat3D)); zmax = max(dat3D(:));
        feval(F.aspectRatios)
    elseif strcmp(Props.view.ax2Type,'contour')         
        contour(Dat.proc.x,Dat.proc.y,Dat.proc.image,Props.ncontours), feval(F.aspectRatios)
        set(gca, 'YDir', 'reverse');
    elseif strcmp(Props.view.ax2Type,'inten')
        imagesc(Dat.proc.x,Dat.proc.y,Dat.proc.image,Dat.proc.Irng), feval(F.aspectRatios)
    elseif strcmp(Props.view.ax2Type,'ft')
        ft = abs(fftshift(fft2(fftshift(double(Dat.proc.image))))); 
        imagesc(ft), feval(F.aspectRatios)
    else
        set(findobj(cmdHndl,'string','Can''t, you ninny...')), return
    end
    
    feval(F.axset), Props.view.ax2 = gca;
    if ~strcmp(Props.view.ax2Type,'ft'), , end
end

feval(F.cbarsSet,'view')
if ~isempty(logInd)
    title(tit(1:logInd-2))
else, title(tit), end

dispWhy(cmdHndl)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'log' % log radio button
  
axType = 'log';

datLog = Dat.proc.image;                        %% make sure double
ind0 = find(datLog <= 0);
datLog(ind0) = NaN;
if ind0, set(cmdHndl,'string','Note: ignoring zero and negative numbers'), end
set(Props.view.tag.lin,'Value',0); %% radiobuttons
set(Props.view.tag.log,'Value',1)

if Props.isAx1
    axes(Props.view.ax1), tit = get(get(gca,'Title'),'String');
    if strcmp(Props.view.ax2Type,'profile')
        datProfLog = Dat.profile.dat;
        if ind0,datProfLog(find(datProfLog <= 0)) = NaN; end
		  axes(Props.view.ax2), semilogy(Dat.profile.xDat,datProfLog),
		  ylabel('Intensity'), xlabel('Position in profile'), title('Image Profile [log]'); 
% % keep the figure current so you can slice multiple times
    axes(Props.view.ax1)
    end
	if strcmp(Props.view.ax2Type,'hist')
		set(Props.view.ax2,'yscale','log'), end
    if strcmp(Props.view.ax1Type,'ift')
        if get(Props.view.tag.isPhase,'value')
            set(cmdHndl,'string','Can''t plot log of negative numbers... '), return
        else
            datFoo = log10(abs  (Dat.ift.image));
        end
    else  
            datFoo = log10(datLog);
    end
    imagesc(Dat.proc.x,Dat.proc.y,datFoo)
%     imagesc(Dat.proc.x,Dat.proc.y,log(datLog)),
    feval(F.aspectRatios)
    if Props.isCbar, Props.view.cb1 = colorbar; else Props.view.cb1 = []; end
    feval(F.drawProfLines,Dat.profile); 
    feval(F.axset) %% set current axes
    Props.view.ax1 = gca;
else
    axes(Props.view.ax2), tit = get(get(gca,'Title'),'String');
    if strcmp(Props.view.ax2Type,'plot3')
        dat3D = datLog;
        set(gca,'ZScale',axType)
        xmax = max(size(dat3D)); zmax = max(dat3D(:));
        feval(F.aspectRatios)
    elseif strcmp(Props.view.ax2Type,'contour')         
        contour(Dat.proc.x,Dat.proc.y,log10(datLog),Props.ncontours), feval(F.aspectRatios) 
        set(gca, 'YDir', 'reverse');
    elseif strcmp(Props.view.ax2Type,'inten')
        imagesc(Dat.proc.x,Dat.proc.y,log10(datLog),log10(Dat.proc.Irng)), feval(F.aspectRatios)
    elseif strcmp(Props.view.ax2Type,'ft')
        imagesc(Dat.proc.x,Dat.proc.y,log10(datLog)), feval(F.aspectRatios)
    else
        set(cmdHndl,'string','Can''t, you ninny...'), return
    end
    
    if Props.isCbar, Props.view.cb2 = colorbar; else Props.view.cb2 = []; end
    feval(F.axset) %% set current axes
    Props.view.ax2 = gca;
    if ~strcmp(Props.view.ax2Type,'ft'), , end
end %% checking ax1 or ax2

title([tit,' [log]'])

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'cmd' % % enter text in command line
    str     = get(gcbo,'String');
    if iscell(str)
        str = str{:}; end
    isError = 'NaN';
    
    if get(Props.view.tag.isCmdOutput,'value')
        ans = 0;
        eval(str,isError);
        set(cmdHndl,'string','Done: did it work?  (Checkbox off for numeric return)')
        if isnan((ans)) 
            set(cmdHndl,'string','Error: Invalid Expression')
        end
    else
        out = eval(str,isError);

        if isnan(out)
            set(cmdHndl,'string','Error: Invalid Expression (try the checkbox ...)')
        else
            set(cmdHndl,'string',['ans= ', num2str(out,5)])
        end
    end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'cla2' % % clear axis 2
    axes(Props.view.ax2), cla, set(gca,'visible','off'), title('')
    if ishandle(Props.view.cb2), set([gca,Props.view.cb2],'Visible','off'), end
    axes(Props.view.ax1), dispWhy(cmdHndl)
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'export' % % writes variable to workspace

mess = sprintf('If Variable ''Dat'' exists in workspace\n          it will be overwritten ...\nDo you want to continue?');
button = questdlg(mess,...
	'Continue Operation','Yes','No','Help','Yes');
switch button
    case 'Yes'
        set(cmdHndl,'string','Wrote ''Dat'' to workspace')
        assignin('base','Dat',Dat);
    case 'No'
        set(cmdHndl,'string','Canceled export ...')
    case 'Help'
	set(cmdHndl,'string','Sorry- you''re on your own, poncho...')
end
  

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'newfig' % % opens a new figure with current plot
    
feval(F.newfig)       
                     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'close' % % closes any number of images

datNamesHeader  = datNames(1:2);

[closeFiles,ok] = listdlg('PromptString','Choose Files to Close',...
                          'liststring',datNames(3:end));

if ok %% picked some files to close
    lenCloseFiles          = length(closeFiles);
    Props.nImages          = Props.nImages - lenCloseFiles;
    Dats(closeFiles)       = [];
    datNames(closeFiles+2) = [];
   
    feval(F.clearIndicators,'all')
    
else %% no files selected
    return
end

if Props.nImages %% there are some images left

    Props.currImage = 1;
    Dat             = Dats(Props.currImage);
   
    set(Props.view.fig,         'Name'  ,['Image= ',Dat.fname])
    set(Props.view.tag.datNames,'value' ,Props.currImage + 2,...
                                'string',datNames)
                            
    feval(F.localPlot,'orig');                            
    dispWhy(cmdHndl)
else %% all are removed
    Props.currImage = 0;
    axes(Props.view.ax2), cla, set(gca,'visible','off'), title('')
    axes(Props.view.ax1), cla, set(gca,'visible','off'), title('')
    set(Props.view.tag.datNames,'value' ,1,...
                                'string',datNames)
	set(cmdHndl,                'string','Load more images to continue...')
    clear global Dat Dats
    
    if ishandle(Props.view.cb1)&~isempty(Props.view.cb1), delete(Props.view.cb1), end
    if ishandle(Props.view.cb2)&~isempty(Props.view.cb2), delete(Props.view.cb2), end
    return % return if all files closed
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'next' % % goes to next image
      
if Props.currImage == Props.nImages
    set(cmdHndl,'string','You''re pushing the limits already...'), return
end

Props.currImage 	= Props.currImage + 1;
Dat             	= Dats(Props.currImage);
Dat.orig.image  	= double(Dat.orig.image);
Dat.orig.imageRaw	= Dat.orig.image;
feval				(F.background_subtract)
feval				(F.nextPlot,DatPrev);
dispWhy				(cmdHndl)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'prev' % % goes to previous image

if Props.currImage == 1
    set(cmdHndl,'string','You''re pushing the limits already...'), return
end

Props.currImage 	= Props.currImage - 1;
Dat             	= Dats(Props.currImage); 
Dat.orig.image  	= double(Dat.orig.image);
Dat.orig.imageRaw	= Dat.orig.image;
feval				(F.background_subtract)
feval				(F.nextPlot,DatPrev);
dispWhy				(cmdHndl)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'datname' % %  goes to selected image

if Props.isMovie
    Props.currImage = Props.movie.currFrame;
else
    val = get(gcbo,'value');
    if val < 3, return, end
    Props.currImage = val - 2;
    dispWhy(cmdHndl)
end

Dat            		= Dats(Props.currImage); 
Dat.orig.image 		= double(Dat.orig.image);
Dat.orig.imageRaw	= Dat.orig.image;
feval				(F.background_subtract)
feval				(F.nextPlot,DatPrev);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'calLines' %%  opens Calibration_Lines.pdf ,...
                %% or gives a prompt to help compute wavelength dispersion
                %% depends on 'isNm' checkbox
                
if ~get(Props.view.tag.isNm,'value') %% opens cal lines file
    calName = 'Calibration_Lines.pdf';
    if exist(calName)
        wd  = cd;
        foo = which(calName);
        foo = fileparts(foo);
        cd(foo);
        winopen(calName);
        cd(wd);
        set(cmdHndl,'string','Opening Calibration_Lines.pdf ...')
    else
        set(cmdHndl,'string','You need to have Calibration_Lines.pdf in the matlab path for this feature...' )
    end
else %% open cal dialog
    promptStr = {'[Pixel1, Wavelength1 : Pixel2 Wavelength2]'}; 
    titleStr  = 'Calibration of Linear Dispersion'; nLines = 1; 
    xPix      = str2num(get(findobj(gcf,'tag','nmPix'),'string'));
    xCal      = str2num(get(findobj(gcf,'tag','nmCal'),'string'));
    iniStr    = {num2str([xPix,xCal])};
    set(0,'defaultfiguretoolbar','none')                              %% Must turn it off for dialog to work properly
    CalVals   = inputdlg(promptStr,titleStr,nLines,iniStr);           %% Sets range equal to input dialog that allows user to type in range
    set(0,'defaultfiguretoolbar','figure') 
    if isempty(CalVals), return, end
    set(0,'defaultfiguretoolbar','figure')                            %% set def
    calVals   = abs(str2num(char(CalVals)));
    dispVal   = (calVals(2) - calVals(4)) / (calVals(1) - calVals(3));
    set(Props.view.tag.nmDisp,'string',dispVal);
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'shift' % %  performs cyclic shift on processed image
        
feval(F.shift);
feval(F.localPlot,'roi')                                    %% plot results

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'movie' % % creates movie with selected images
    
datNames = get(Props.view.tag.datNames,'string');    
datNamesHeader = datNames(1:2);

movie            = Props.movie;
[movie.files,ok] = listdlg('PromptString','Choose images for movie ',...
                          'liststring',datNames(3:end));    
movie.path       = Props.pname;

if ok %% picked some files to close                    
    movie.figH = open('iViewMovieGui.fig');
    set(movie.figH,'color',[0.8 0.8 0.85])
else
    return
end

% find children and relevent handles and tags, save to global for easy access %
hndls       = allchild(gcf);
tags        = get(hndls,'tag');
movie.tag   = feval(F.createTags,hndls,tags);
Props.movie = movie;   
feval(F.refreshTags,'movie',{'str','txt','val'})

set(Props.movie.tag.fileName,'string',datNames{2+Props.movie.files(1)})
                      

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'colordef' % % changes colordef property (default figure colors)
    
    if get(gcbo,'value')
        colordef black
        set(gcbf,'color',[0 0 0])
    else
        colordef white
        set(gcbf,'color',[0.8 0.8 0.85])

    end        

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'isDrk' % %background/scale button
%- -%
Props.isDrk = ~Props.isDrk;
feval(F.background_subtract)
dispWhy(cmdHndl)
	
	
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'drk' % %background/scale button
  
set				(0,'defaultfiguretoolbar','figure')
Props.drk.fig 	= open('iViewDrkGui.fig');
Props.drk.tag 	= guihandles(Props.drk.fig);
set				(gcf,'color',[0.6 0.6 0.6])
 

% if Props.isDrk
% check if the actual data exists, not whether the checks are on
if isfield(Props.drk,'drk')
	axes(Props.drk.tag.axes1)
	imagesc(Props.drk.drk), axis off, colorbar
	set(Props.drk.tag.drkDisp,'string', Props.drk.drkname)
	set(Props.drk.tag.isMDrk, 'value', 1)
end
if Props.drk.isPrpl
	axes(Props.drk.tag.axes2)
	plot(Props.drk.bins,Props.drk.counts);	 
	set(Props.drk.tag.isPrpl, 'value', 1)
	set(Props.drk.tag.threshDisp,'string',num2str(Props.drk.thresh))
end
if isfield(Props.drk,'ref1')
	axes(Props.drk.tag.axes3)
	imagesc(Props.drk.ref1), axis off, colorbar
	set(Props.drk.tag.ref1Disp,'string', Props.drk.ref1name)
	set(Props.drk.tag.isMRef1, 'value', 1)
end
if isfield(Props.drk,'ref2')
	axes(Props.drk.tag.axes4)
	imagesc(Props.drk.ref2), axis off, colorbar
	set(Props.drk.tag.ref2Disp,'string', Props.drk.ref2name)
	set(Props.drk.tag.isMRef2, 'value', 1)

end
% end

axes(Props.drk.tag.axes5), axis off

dispWhy(cmdHndl)	

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'rotate' % %  
	
Props.rotate.fig 		= open('iViewRotate.fig');
Props.rotate.tag 	 	= guihandles(Props.rotate.fig);	
Props.rotate.caxis		= get(Props.view.ax1,'clim');

axes(Props.rotate.tag.ax1)
imagesc(Dat.proc.image)

dispWhy(cmdHndl)	


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'designFilter' % %  
    	
Props.filt.fig = open('iViewFilter.fig');
Props.filt.tag = guihandles(Props.filt.fig);


if get(Props.view.tag.colordef,'value')
    colordef black
    set(Props.filt.fig,'color',[0 0 0],...
						'menubar', 'figure'...
   )
else
    colordef white
	set(Props.filt.fig, 'menubar', 'figure',...
 	                   'Color'  , [0.8 0.8 0.85]...
   )
end        

set(Props.filt.tag.order, 'string',num2str(Props.filt.order));
set(Props.filt.tag.cutoff,'string',num2str(Props.filt.cutoff));   
set(Props.filt.tag.xyratio,'string',num2str(Props.filt.xyratio));   

axes    ( Props.filt.tag.axes1)
imagesc ( log10( abs( double( Dat.proc.image))))
colormap( Props.cmap )
feval   ( F.aspectRatios)
Props.filt.cb2 = [];

if Props.isCbar
    Props.filt.cb1 = colorbar; end

feval(F.cbarsSet, 'filt')
title('Shifted FFT [log]')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'ift' % %  

Dat.ift       = Dat.proc;
% Dat.ift.image = double(Dat.ift.image);                % make sure not uint8

if get(Props.view.tag.isFilter,'value')
     if ~isfield(Props.filt,'H')
        set(cmdHndl,'string', 'Please design a filter first...'), return
     end
    if size(Dat.ift.image) ~= size(Props.filt.H)
        set(cmdHndl,'string', 'You need a new filter, jerky ...'), return
    else
        foo = Dat.ift.image .* abs(Props.filt.H);
    end
else
    foo = Dat.ift.image;
end

Dat.ift.image = fftshift(ifft2(fftshift(foo)));

feval(F.localPlot,'ift')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'exit' % % exit viewer and clear variables
    
mess = sprintf('Do you want to exit (clears all globals)?');
button = questdlg(mess,...
	'Continue Operation','Yes','No','Help','Yes');
    
switch button
    case 'Yes'

		rmpath(Props.iViewDir)
		
        if exist('Props.movie.figH')
            close(Props.movie.figH), end
        if exist('Props.filt.fig')
            close(Props.filt.fig), end
        close(gcbf)        
        clear global F Props Dat Dats
		
		colordef white
        return
    case 'No'
        set(cmdHndl,'string','Haven''t had enough YET?'), return
    case 'Help'
        set(cmdHndl,'string','Sorry- you''re on your own, poncho...'), return
end
        


end %% switch
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
% % finish program, set colorbars, userdata

if ~strcmp(cback,'movie') 
    if ~Props.isWriteImages %don't want colobar on movie gui !
        feval(F.cbarsSet,'view') %% set colorbar widths
    end
end

if exist('Dat')
    Props.prevUnits   = Dat.units;
else
    Props.prevUnits.x = 'pixels'; 
    Props.prevUnits.y = 'pixels'; 
end

if Props.isMm %% change to mm labels
    set(Props.view.tag.xAxLab,'string','x (mm)')
    set(Props.view.tag.yAxLab,'string','y (mm)')
	 set(Props.view.tag.profFWHMLab,'string','FWHM (mm)')
else 
    set(Props.view.tag.xAxLab,'string','x (pixels)')
    set(Props.view.tag.yAxLab,'string','y (pixels)')
	 set(Props.view.tag.profFWHMLab,'string','FWHM (pixels)')
end

if Props.isNm %% change to nm labels
    set(Props.view.tag.xAxLab,'string','x (nm)')
    set(Props.view.tag.profFWHMLab,'string','FWHM (nm)')
end

warning on                                    % turns warnings back on

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% local ONLY subfunctions (you sure? err to iViewFunctions.m)%
% % % %
function dispWhy(cmdHndl)
set(cmdHndl,'string',whyUNO)
