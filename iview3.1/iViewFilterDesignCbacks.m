


function iViewFilterDesignCbacks(cback) %% input callback switch

debug = 0;

global F Props Dat

% set the values % 
Props.filt.order    = str2double(get(Props.filt.tag.order,'string'));
Props.filt.xyratio  = str2double(get(Props.filt.tag.xyratio,'string'));   
Props.filt.cutoff   = str2double(get(Props.filt.tag.cutoff,'string')); 
Props.filt.window   = get(Props.filt.tag.window,'value'); 
Props.filt.method   = get(Props.filt.tag.method,'value'); 
Props.filt.highpass = get(Props.filt.tag.highpass,'value');     

switch cback

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'apply' % % shows filtered plot
      
foo  = abs(Props.filt.H) .* double(Dat.proc.image);      % make sure double (not uint8)
iFFT = fftshift(ifft2(fftshift(foo)));
axes ( Props.filt.tag.axes2)
    
if get( Props.filt.tag.isPhase, 'value')
    imagesc(angle(iFFT))
else
    imagesc(abs(iFFT))
end

title('IFFT')
colormap(Props.cmap), Props.filt.cb2 = [];

if Props.isCbar
   Props.filt.cb1 = colorbar; end

feval(F.cbarsSet,'filt')
feval(F.aspectRatios)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'show' % % overlays filter over plot
    
    
axes    ( Props.filt.tag.axes1)
cla
imagesc ( log10( abs( double( Dat.proc.image))))
hold on
colormap( Props.cmap )
feval   ( F.aspectRatios)
Props.filt.cb2 = [];

if Props.isCbar
    Props.filt.cb1 = colorbar; end

feval(F.cbarsSet, 'filt')
title('Shifted FFT [log]')


imageSize           = size(Dat.proc.image);
[h, Props.filt.H]   = feval(F.DesignFilter, imageSize);
    
Props.filt.H        = Props.filt.H' ; % Switched the X and Y before ... ?
contour(log10(abs(Props.filt.H)),'k')
	
axes(Props.filt.tag.axes2)
imagesc(log10(abs(Props.filt.H)))
colormap(Props.cmap)
Props.filt.cb2      = [];
feval(F.aspectRatios)

if Props.isCbar
   Props.filt.cb1   = colorbar; end
feval(F.cbarsSet,'filt'), title('Filter Response [log]')
    
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
case 'newfig' % % opens a new figure with current plot

newFig = figure; 

newPlot = copyobj(Props.filt.tag.axes2,newFig);
axes(newPlot)
set(newPlot,'Position',[0.1 0. 0.8 0.8])

feval(F.aspectRatios)
if Props.isCbar
    cb = colorbar; 
    cbPos = get(cb,'position'); % set colorbar positions
    set(cb,'position',[cbPos(1:2) 0.01 cbPos(4)]);
end
    colormap(Props.cmap)
    

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'close' % % closes figure 

    close(gcbf);
    
end














