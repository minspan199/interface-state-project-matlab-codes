


function iViewRotate(cback) %% input callback switch

debug = 0;

global F Props Dat

% set the values % 
isReplot	= 0;

wRng	= get(Props.rotate.tag.range,'value');
rng		= [-5 5 ; -90 90 ; -180 180];

rMin	= rng(wRng,1);
rMax	= rng(wRng,2);

switch cback
	
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'range' % % shows filtered plot	
	
set(Props.rotate.tag.angleSlider,'min',rMin,'max',rMax)
set(Props.rotate.tag.max,'string',num2str(rMax))
set(Props.rotate.tag.min,'string',num2str(rMin))

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'angle' % % shows filtered plot
	
isReplot	= 1;

sAngle		= str2num(get(Props.rotate.tag.setAngle,'string'));

if sAngle <= rMin
	set(Props.rotate.tag.angleSlider,'value',rMin)
elseif sAngle >= rMax 
	set(Props.rotate.tag.angleSlider,'value',rMax)
else
	set(Props.rotate.tag.angleSlider,'value',sAngle)
end

set(Props.rotate.tag.currAngle,'string',num2str(sAngle))

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'angleSlider' % % shows filtered plot
	
isReplot	= 1;
	
set(Props.rotate.tag.currAngle,'string',...
	num2str(get(Props.rotate.tag.angleSlider,'value')))



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'onoff' % % shows filtered plot	
	
isReplot	= 1;	
	

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'close' % % closes figure 

    close(gcbf);
    break
	
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
case 'accept' % % closes figure 
	
% angle	= str2num(get(Props.rotate.tag.currAngle,'string'));
% Dat.proc.image = imrotate(Dat.proc.image,angle,'bilinear','crop');
close(gcbf);
	break	
	
	
end % case


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

angle	= str2num(get(Props.rotate.tag.currAngle,'string'));

set(Props.rotate.tag.currAngle,'string',num2str(angle))

if isReplot

if 	get(Props.rotate.tag.onoff,'value')
	foo		= imrotate(Dat.proc.image,angle,'bilinear','crop');
else
	foo		= Dat.proc.image;
end

imagesc(foo)
title('Rotated Image')
colormap(Props.cmap), Props.filt.cb2 = [];

if Props.isCbar
    cb = colorbar; 
    cbPos = get(cb,'position'); % set colorbar positions
    set(cb,'position',[cbPos(1:2) 0.01 cbPos(4)]);
end
	
end

caxis(Props.rotate.caxis)

% if Props.isCaxis, 
%     caxis( Props.cAx )
% end














