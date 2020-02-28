% iVIEW image display and processing gui
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
% User .M files: Many... WINSPECREAD, iView (sets GUI data initially)
%                Colormaps: CMAPRAINBOW2, JET2, JET3
% User MEX, C Files, and other: UIGETFILES.DLL, WHYUNO.DLL, Calibration_Lines.pdf
% User .MAT files: none 
% Subfunctions: A few
% See also: iViewer1280x1024.fig, iViiewerXP.fig
%
% * Note: Requires Image Processing Toolbox to be installed *


% Author: Kevin Tetz 
% Last revision: 03-Feb-2005
% Revision History:
% 11-Jul-2006: Add xsection in 2D simultaneous in new window
% 03-Feb-2005: Import from workspace
% 26-Oct-2004: v2.5: Added Drk Gui, fixed phase profiles, added shift via
% 				coords
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

% To add: - Make 3D plots render in new window- much faster for some reason
% 			(need to track axes carefully with this) 
%		  - units on FT
%         - support for multiple peaks (finding, COM, etc.)
%         - ellips eccentricity , etc.
%         - zoom memory (zoom out/prev)
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
%		 _ MakE LIN AND LOG A HELLA LOT EASIER!  JUST SWITCH AXES SCALES
%		 ONLY (SHOULD WORK WITH SOME PRECAUTIONS)
%			- add log scale histogram 
%		- make coords a Prop not Dat, so stays with a new image...
% Known Bugs: - 3-D data aspect ratios w/ different units... z sometimes
%               curious with mixed nm/mm (also code somewhat redundant)
%             - Bkg and Scaling need tweaking (image specific??)
%               get slightly off FWHM when scaled ...                   
%             - some axes tick labels obscured on thin plots
%             - "uigetfiles.dll" seems to only be able to "see" 250 files/dir
%			  - check out FFT values when performing a shift (change
%			  slightly??)
%				- remove colorbars from new figs like hist, prof

function iView(varargin)

vers = '3.1'; 

clear global

global F Props Dat Dats

% add iviewer directory to path for uninterrupted use %
iviewPath	= which('iview');
iviewPath	= iviewPath(1:end-7);

iViewDir = [iviewPath,'iview',vers];
addpath(iViewDir)

if nargin == 1
	viewername = deal(varargin{:});
else
	viewername = 'iviewer.fig';	
end

open(viewername) % open gui
set (gcf,'color',[0.8 0.8 0.85])

iViewFunctions        % sets global function handles for subfunctions %

ax = findobj(gcf,'type','axes');
view = struct('ax1',ax(2),...
              'ax2',ax(1),...
              'cb1',[],...
              'cb2',[],...
              'ax1Type','none',...
              'ax2Type','none',...
              'isFFT',0,...
              'fig',gcf...
              );

% find children and relevent handles and tags, save to global for easy access %
hndls    = allchild (gcf);
tags     = get      (hndls,'tag');
view.tag = feval    (F.createTags,hndls,tags);

% Movie Properties %
fig.Color = [0.8 0.8 0.8];
txt.Color = [000 000 000]/256;
% txt.FontAngle = 'normal';
txt.FontSize = 30;                   % Label (Text, Title) font size
txt.FontWeight = 'normal';           % Label (Text, Title) weight


val = struct('isTight',0,...
              'isFigset',0,...
              'isFigSize',0,...
              'isText1',0,...
              'isText2',0,...
              'isText3',0,...
              'codec',1,...
              'axOn',1,...
              'what',1);

str = struct(  'movName','Mov1',...
               'path',pwd,...
               'tokInd1',1,...
               'tokInd2',2,...
               'token','',...
               'fps',1,...
               'quality',75,...
               'figX',1000,...
               'figY',1000,...
               'textX',0,...
               'textY',0);
           
movie = struct('txt',txt,'val',val,'str',str,'fig',fig);
           
           
% set filter properties %
filt = struct('order'   , 50,...
              'cutoff'  , 0.15,...
              'window'  , 1,...
              'method'  , 1,...
			  'xyratio' , 1,...
              'highpass', 0);
		  
% set filter properties %
drk = struct('thresh'  	, [],...
			 'isPrpl'	,0,...
			 'isDrk'	,0,...
			 'isRef1'	,0,...
			 'isRef2'	,0);		  

% set global properties %
units.x     = 'pixels';
units.y     = 'pixels';
datnames    = cell(2,1);
datnames{1} = 'Loaded Images';
datnames{2} = '-------------------------';

Props = struct('view',view,...
               'movie',movie,...
               'filt',filt,...
			   'drk',drk,...
               'ncontours',10',...
               'profHndl',[],...
               'nImages',0,...
               'currImage',0,...
               'roi',[0 0 0 0],...
               'nProfPoints',5000,...
               'lowIrng',0.1,...
               'highIrng',0.98,...
               'prevUnits',units,...
               'isMovie',0,...
               'isWriteImages',0,...
			   'isDrk',0,...
			   'iViewDir',iViewDir...
			   );






