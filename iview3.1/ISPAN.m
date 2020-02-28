clc
close all
clear all
vers = '3.1'; 
clear global
global F Props Dat imageSize
% add iviewer directory to path for uninterrupted use %
iviewPath	= which('iview');
iviewPath	= iviewPath(1:end-7);
ImageA = imread('Image1.tif');
ImageB = imread('Image2.tif');

kkinin = 1;
N = 400;
for kk = 1:1:N

    filename = strcat('Image4_',num2str(kk),'.tif');
    iViewDir = [iviewPath,'iview',vers];
    addpath(iViewDir)
	viewername = 'iviewer.fig';	
    iViewFunctions        % sets global function handles for subfunctions %
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
           iFile = 1;
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
		
		assignin	('base','iFile',iFile)
		dat			= evalin('base','eval(zzstr{varNo(iFile)})');
		iFname		= fname{iFile};
								
	else
		dat    = (imread(filename)); 
	end

	if ndims(dat) >2
		dat = mean(dat,3); 
	end
	roi          = [1 size(dat,2) 1 size(dat,1)]; %% opt out-> Props.roi in future... set to Dat.proc.roi as wells
	orig.image   = dat; units.x = 'pixels';
    proc.image   = dat; units.y = 'pixels';
    % find children and relevent handles and tags, save to global for easy access %
    hndls    = allchild (gcf);
    tags     = get      (hndls,'tag');
    view.tag = feval    (F.createTags,hndls,tags);
    
    Dat.orig.imageRaw = double(dat);
    Props.drk.ref1	= double(ImageA);
    Props.drk.isRef1= 1;
    Props.drk.ref2	= double(ImageB);
    Props.drk.isRef2= 1;
    ref1 = Props.drk.ref1;
    ref2 = Props.drk.ref2;
    im	= Dat.orig.imageRaw;
    if Props.drk.isRef1
        im	= im	- ref1;
    end
    if Props.drk.isRef2
        im	= im	- Props.drk.ref2;
    end
    Dat.orig.image	= im;
    
    
    roi = [1 size(dat,2) 1 size(dat,1)]; %% opt out-> Props.roi in future... set to Dat.proc.roi as wells
    units.x = 'pixels';
    units.y = 'pixels';
    [Dat.orig.xPixels,orig.x,proc.x] = deal(linspace(1,roi(2),roi(2)));
    [Dat.orig.yPixels,orig.y,proc.y] = deal(linspace(1,roi(4),roi(4)));
    Dat.proc = Dat.orig;
    Dat.proc.x = proc.x;
    Dat.proc.y = proc.y;
    Dat.orig.x = orig.x;
    Dat.orig.y = orig.y;
    datInit = Dat;   
    Props.roiInds = [1 size(datInit.orig.image,2) 1 size(datInit.orig.image,1)];

                                                            %% Set new range
Dat.proc.x         = Dat.proc.x(Props.roiInds(1):Props.roiInds(2));
Dat.proc.y         = Dat.proc.y(Props.roiInds(3):Props.roiInds(4));
Dat.proc.image     = Dat.proc.image(...                    
                                    Props.roiInds(3):Props.roiInds(4),...
                                    Props.roiInds(1):Props.roiInds(2)...
                                    ) ;
Dat.proc.image = fftshift(fft2(fftshift(Dat.orig.image)));
% figure; imagesc(Dat.orig.image); axis off; colorbar;
% figure; imagesc(abs(Dat.proc.image)); set(gca,'CLim',[10000 1000000]); axis off; colorbar; colormap jet;
%                                                                             shift([91 106]);%for waveguide array in D:\Camera\170503\C2
%                                                                               shift([90 140]);
%                                                                                shift([93 140])%for single waveguide in D:\Camera\170505\D
%                                                                                  shift([93,177]);
%                                                                                  shift([94,138]);
% shift([93,111]);%D:\Camera\170510\B
shift([94,139]);%D:\Camera\170510\A
Props.filt.order    = 200;
Props.filt.xyratio  = 4;   
Props.filt.cutoff   = .05; 
Props.filt.window   = 1; 
Props.filt.method   = 1; 
Props.filt.highpass = 0;  

imageSize           = size(Dat.proc.image);
[h, Props.filt.H]   = feval(@DesignFilter, imageSize);
Dat.ift.image = fftshift(ifft2(fftshift(Dat.proc.image.*(abs(Props.filt.H))')));
% J = imagesc(((abs(Dat.ift.image))),[100 1250]); %for single waveguide in D:\Camera\170505\D
J = imagesc(((abs(Dat.ift.image))),[200 1300]); %for single waveguide in D:\Camera\170505\
I = ColorRangeAdujust(J.CData,200);

hold on
% Intensity = mean(I.CData);
% plot(Intensity);
% for TT = Waveguide(1):1:(Waveguide(1)+Waveguide(3))
%     J = I.CData((Waveguide(2)-Waveguide(4)):(Waveguide(2)),Waveguide(1):(Waveguide(1)+Waveguide(4)));
% figure
% imshow(J)
% tim(kk) = ((kk)-kkinin-7)*3.33;
tim(kk) = ((kk)-kkinin-130)*3.33;
if tim(kk) > 0
    text(80,220,strcat('+',num2str(tim(kk)),'  fs'),'Color','white','FontSize',16)
else
    text(80,220,strcat(num2str(tim(kk)),'  fs'),'Color','white','FontSize',16)
end
axis off
colorbar
colormap jet
% set(gca, [100 1200]);

Intensity(kk,:) = flipud(mean(I)')';
plot(mean(I),'w-');
hold off

% imwrite(tmp2.cdata,strcat('Filtered_image1_',num2str(kk+10),'.tif'));
saveas(gca,strcat('Filtered_image13_',num2str(kk+10),'.tif'));

end


figure
imagesc(Intensity)
h = get(gca,'xtick');
set(gca,'xticklabel',(h-153)*0.29776);
h = get(gca,'ytick');
set(gca,'yticklabel',(h-131)*3.33);
colorbar
colormap jet
xlabel('Position of Waveguide(um)')
ylabel('Time Elapse(fs)')
title('Position and Time of Distribution of the Pulse')



x = (0:0.01:320)*0.29776;
xx = ((1:1:320))*0.29776;
for kk = 1:1:N
    
    g{kk} = fit(xx.',Intensity(kk,:).','gauss2');
    ff(kk,:) = g{kk}.a1*exp(-((x-g{kk}.b1)/g{kk}.c1).^2) + g{kk}.a2*exp(-((x-g{kk}.b2)/g{kk}.c2).^2);
    [m(kk),Zc(kk)] = max(ff(kk,:));
    Zc(kk) = x(Zc(kk));

end


figure
plot(xx,(Intensity(200,:).^2)/max(Intensity(200,:).^2),'red')
hold on;plot(xx,Intensity(210,:).^2/max(Intensity(210,:).^2),'yellow')
% hold on;plot(xx-26.2029,Intensity(55,:).^2/max(Intensity(55,:).^2),'red')
hold on;plot(xx,Intensity(220,:).^2/max(Intensity(220,:).^2))
% hold on;plot(xx-26.2029,Intensity(65,:).^2/max(Intensity(65,:).^2),'red')
% hold on;plot(xx-26.2029,Intensity(70,:).^2/max(Intensity(70,:).^2),'red')
% hold on;plot(xx-26.2029,Intensity(75,:).^2/max(Intensity(75,:).^2),'red')
hold on;plot(xx,Intensity(190,:).^2/max(Intensity(190,:).^2),'green')
xlabel('waveguide position from 0um')
ylabel('normalized intensity')
% Create textarrow
annotation('textarrow',[0.374359518154055 0.452930946725483],...
    [0.723809523809524 0.722809523809525],'String',{'+123.2fs'});

% Create textarrow
annotation('textarrow',[0.388662198846284 0.467233627417712],...
    [0.616666666666668 0.615666666666669],'String',{'+139.86fs'});

% Create textarrow
annotation('textarrow',[0.6833729216152 0.603015778758058],...
    [0.590476190476192 0.59047619047619],'String',{'+239.76fs'});

% Create textarrow
annotation('textarrow',[0.431485408890397 0.510056837461825],...
    [0.555050922102947 0.554050922102948],'String',{'+189.81fs'});


    figure
    plot(xx,Intensity(210,:))
    hold on
    plot(g{210})
    xlabel('Position of Waveguide(um)')
    ylabel('Amplitude')
%     h = get(gca,'xtick');
%     set(gca,'xticklabel',(h-153)*0.29776);

    
    x=-12:0.01:500;
    for kk = 1:1:imageSize(2)
        gg{kk} = fit(tim.',Intensity(:,kk),'gauss2');
        f(kk,:) = gg{kk}.a1*exp(-((x-gg{kk}.b1)/gg{kk}.c1).^2) + gg{kk}.a2*exp(-((x-gg{kk}.b2)/gg{kk}.c2).^2);
        [m(kk),Tc(kk)] = max(f(kk,:));
%         [m(kk),Tc(kk)] = max(Intensity(:,kk));
        Tc(kk) = x(Tc(kk));
    end
    
    
    
    figure
    plot(tim.',Intensity(:,201))
    hold on
    plot(gg{201})
    ylabel('Amplitude at Waveguide Position 11.62um')
    xlabel('Delayline Position(fs)')
    
    figure
    imagesc(Intensity)
    colorbar
    colormap jet
    xlabel('pulse position in the waveguide(um)')
    ylabel('delayline position(fs)')
    title('Field Amplititude Distribution')
    
    % Ran = (unidrnd(51,[1 2])+39);
    clear KP
    k = 1;
for kk = 1:1:50000
    Ran = (unidrnd(50,[1 2])+202);
    a = min(Ran);
    b = max(Ran);
    if b-a>40
        KP(k,1) = a; 
        KP(k,2) = b;
        P = polyfit(a:b,Zc(a:b),1);
        KP(k,3) = P(1);
        KP(k,4) = P(2);
        k = k+1;
    end
end
inc = mean(1./(KP(:,3)))
ub = mean(1./(KP(:,3)))-min(1./(KP(:,3)))
lb = mean(1./(KP(:,3)))-max(1./(KP(:,3)))

figure
P = polyfit(185:228,Zc(185:228),1)
plot(Zc(185:228),'.')
hold on
plot(P(1)*(185:228)+P(2))
h = get(gca,'xtick');
set(gca,'xticklabel',(h+50-8)*3.33);
xlabel('delayline position(fs)')
ylabel('pulse center(um)')

figure
plot(Zc,'.')
h = get(gca,'xtick');
set(gca,'xticklabel',(h-8)*3.33);
h = get(gca,'ytick');
set(gca,'yticklabel',(h-109*0.2992));
xlabel('delayline position(fs)')
ylabel('pulse position in the waveguide(um)')






















function computeFFT    % simple version ... may make a new variable Dat.fft.image in future ...

global Dat

Dat.proc.image = fftshift(fft2(fftshift(Dat.proc.image))); 

dx				= diff(Dat.orig.x);
dx				= dx(1);
dy				= diff(Dat.orig.y);
dy				= dy(1);
end
% Dat.proc.x		= linspace(-(1./dx)/2,(1./dx)/2,length(Dat.orig.x));
% Dat.proc.y		= linspace(-(1./dy)/2,(1./dy)/2,length(Dat.orig.y));

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
end



function shift(scCOM)

global Dat Props

imageCenter      = [(Dat.proc.y(end) - Dat.proc.y(1)) /2,...    % [Y, X]
                    (Dat.proc.x(end) - Dat.proc.x(1)) /2   ] + ...
                   [Dat.proc.y(1) Dat.proc.x(1) ] ;             % offset by roi       


		Dx           = round( imageCenter - scCOM );
    
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
end

function  [h, r] = DesignFilter(ImageSize)

global Props

xyratio	  = Props.filt.xyratio;
order     = Props.filt.order;
cutoff    = Props.filt.cutoff;
% windStr   = get(Props.filt.tag.window,'string');
wind      = 'Blackman';
% methodStr = get(Props.filt.tag.method,'string');
method      = 'fwind1';
Highpass  = Props.filt.highpass;

%Create desired frequency responce
[f1,f2]     =       freqspace(order, 'meshgrid');
d           =       find(f1.^2 + f2.^2 < cutoff^2);
Hd          =       zeros(order);
Hd(d)       =       1;
if Highpass
    Hd          =       1-Hd;   %Highpass filter
end
%Design 

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
end


function I = ColorRangeAdujust(J,n)
T = size(J);
I = J;
for k = 1:1:T(1)
    for kk = 1:1:T(2)
        if I(k,kk) < n
            I(k,kk) = 0;
        end
    end
end
end
    

function [Rxx]=autom(x)
% [Rxx]=autom(x)
% This function Estimates the autocorrelation of the sequence
of
% random variables given in x as: Rxx(1), Rxx(2),…,Rxx(N),
where N is
% Number of samples in x.
N=length(x);
Rxx=zeros(1,N);
for m=1: N+1
for n=1: N-m+1
Rxx(m)=Rxx(m)+x(n)*x(n+m-1);
end
end
end