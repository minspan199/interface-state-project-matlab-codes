


function iViewFilterDesignCbacks(cback) %% input callback switch

debug = 0;

global F Props Dat

% set the values % 
Props.drk.isLog		= get(Props.drk.tag.log,	'value');
Props.drk.isMDrk 	= get(Props.drk.tag.isMDrk, 'value');
Props.drk.isMRef1 	= get(Props.drk.tag.isMRef1,'value');
Props.drk.isMRef2 	= get(Props.drk.tag.isMRef2,'value');
Props.drk.isPrpl 	= get(Props.drk.tag.isPrpl, 'value');

switch cback

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'drk' % % shows filtered plot

[dat,fname]	= load_file;
	
set(Props.drk.tag.drkDisp,'string', fname)

axes			(Props.drk.tag.axes1)
imagesc			(dat), axis off, colorbar
axes			(Props.drk.tag.axes2)
[counts,bins]	= hist(double(dat(:)),256); %% only 256 bins here ... can alter
histHnd			= plot(bins,counts);

set(Props.drk.tag.isMDrk,'value',1)

Props.drk.counts= counts;
Props.drk.bins	= bins;
Props.drk.drk	= double(dat);
Props.drk.isDrk	= 1;
Props.drk.drkname	= fname;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'thresh' % % overlays filter over plot

axes(Props.drk.tag.axes2)
plot(Props.drk.bins,Props.drk.counts);	
if Props.drk.isLog
	set(gca,'YScale','log'), end

[foo, foo2]	= ginput(1);
Props.drk.thresh	= foo;
hold on

l1 = line([foo foo], [1 ,max(Props.drk.counts)]);
set(l1,'color','r','linestyle','--','linewidth',2)
hold off

set(Props.drk.tag.threshDisp,'string',num2str(foo))
set(Props.drk.tag.isPrpl,'value',1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'threshEdit' % % overlays filter over plot

thresh	= str2num(get(Props.drk.tag.threshDisp,'string'));
axes(Props.drk.tag.axes2)
plot(Props.drk.bins,Props.drk.counts);	
if Props.drk.isLog
	set(gca,'YScale','log'), end

hold on
l1 = line([thresh thresh], [1 ,max(Props.drk.counts)]);
set(l1,'color','r','linestyle','--','linewidth',2)
hold off

Props.drk.thresh	= thresh;
set					(Props.drk.tag.isPrpl,'value',1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'ref1'
	
[dat,fname]	= load_file;
set(Props.drk.tag.ref1Disp,'string', fname)
axes			(Props.drk.tag.axes3)
imagesc			(dat), axis off, colorbar
Props.drk.ref1	= double(dat);
Props.drk.ref1name	= fname;
Props.drk.isRef1= 1;
set(Props.drk.tag.isMRef1,'value',1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'ref2'
	
[dat,fname]	= load_file;	
set(Props.drk.tag.ref2Disp,'string', fname)
axes			(Props.drk.tag.axes4)
imagesc			(dat), axis off, colorbar
Props.drk.ref2	= double(dat);
Props.drk.ref2name	= fname;
Props.drk.isRef2= 1;
set(Props.drk.tag.isMRef2,'value',1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'prplChk'
	
% Props.drk.isPrpl = 1;
% set(Props.drk.tag.isMRef2,'value',1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'drkChk'
	
% Props.drk.isRef2= 1;
% set(Props.drk.tag.isMRef2,'value',1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'ref1Chk'
	
% Props.drk.isRef2= 1;
% set(Props.drk.tag.isMRef2,'value',1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'ref2Chk'
	
% Props.drk.isRef2= 1;
% set(Props.drk.tag.isMRef2,'value',1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'lin' % % opens a new figure with current plot
	axes(Props.drk.tag.axes2)
	set(gca,'YScale','lin')
	set(Props.drk.tag.log,'value',0)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'log' % % opens a new figure with current plot	
	axes(Props.drk.tag.axes2)
	set(gca,'YScale','log')
	set(Props.drk.tag.lin,'value',0)
	
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'isLinAll' % % opens a new figure with current plot
	set(Props.drk.tag.isLogAll,'value',0)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'isLogAll' % % opens a new figure with current plot	
	set(Props.drk.tag.isLinAll,'value',0) 
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'display' % % opens a new figure with current plot	

% feval  (F.background_subtract)	
% axes   (Props.drk.tag.axes5)
% if get(Props.drk.tag.isLogAll,	'value')
%     imagesc(log10(Dat.orig.image))
% else
%     imagesc(Dat.orig.image)
% end
% axis off, colorbar
				
	
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'close' % % closes figure 
	if Props.drk.isPrpl | Props.drk.isDrk | ...
			Props.drk.isRef1 | Props.drk.isRef2
		Props.isDrk	= 1;
		set(Props.view.tag.isDrk,'value',1)
	end
	
    close(gcbf);
    iViewerCbacks('orig')
    return
end

%- plot data in display window -%
Props.isDrk = 1;
feval  (F.background_subtract)	
axes   (Props.drk.tag.axes5)
if get(Props.drk.tag.isLogAll,	'value')
    datLog = Dat.orig.image;                        %% make sure double
    ind0 = find(datLog <= 0);
    datLog(ind0) = NaN;
    imagesc(log10(datLog))
else
    imagesc(Dat.orig.image)
end
axis off, colorbar

%- Subfunctions -%
function [dat,fname] = load_file

global Props

wd	= cd;	
if isfield(Props,'pname') %% give start path as last directory
    [fname, pname] = uigetfile('*.hdf;*.pcx;*.xwd;*.ico;*.cur;*.gif;*.tif;*.jpg;*.bmp;*.png;*.dat;*.spe',...
        'Open Image File',Props.pname);   
else
    [fname, pname] = uigetfile('*.hdf;*.pcx;*.xwd;*.ico;*.cur;*.gif;*.tif;*.jpg;*.bmp;*.png;*.dat;*.spe', 'Open Image File');   
end

filename = fullfile(pname, fname);

if isequal(fname,0)|isequal(pname,0)
		disp('File not found'), return, end

cd(pname);
if findstr(filename,'.dat') %% for ascii data files
	dat    = load(filename);
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
	else
		dat    = (imread(filename)); 
	end
cd(wd)


% % change color images to intensity only
if ndims(dat) >2
	dat = mean(dat,3); 
end










