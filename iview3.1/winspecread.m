%WINSPECREAD reads binary data file to vector or matrix
% Read WinSpec32 (Princeton Instruments/Roper Scientific)
% Binary Data file
%
% Syntax: dat = winspecread(isPlot,filename);
%
% Inputs: isNative: 1- same format as .spe file, 0- double
%         isPlot: plot (image) the file in question (1- yes, 0- no)
%         filename: optional input file name
%
% Outputs: DAT is the output matrix (or vector), with double precision
%
% User .M files: none
% User .MAT files: none
% Subfunctions: none
%
% See also: none
%
% dat = winspecread;

% Author: Kevin Tetz
% Last revision: 07-Nov-2003
% Revision History:
% To add:
% Known Bugs:

function varargout = winspecread(varargin);

isPlot = 0; %% doesn't plot
isNative = 1; %% reads in in whatever format the data is saved in
nFrames = 1; %% default to single frame

if mod(length(varargin),2) ~= 0
    % input args have not com in pairs, woe is me
    error('Arguments to WINSPECREAD must come param/value in pairs.')
end

for ii=1:2:length(varargin)
    switch lower(varargin{ii})
     case 'isplot'
         isPlot = deal(varargin{ii+1});
     case 'isnative'
         isNative = deal(varargin{ii+1});
     case 'filename'
         filename = deal(varargin{ii+1});
     case 'nframes'
         nFrames = deal(varargin{ii+1});
     otherwise
      error(['Unknown parameter name passed to WINSPECREAD.  Name was ' varargin{ii}])
     end
end
 
wd = cd;

if ~exist(filename)
    [fname, pname] = uigetfile({'*.spe', 'Winspec32 Data Files (*.spe)';'*.*', 'All files (*.*)'}, ...
        'Open Winspec32 data file');
    if isequal(fname,0)|isequal(pname,0)
        error('File not found')
    else
        filename = fullfile(pname,fname);	 
    end
end

fid = fopen(filename, 'r', 'n');

% % check the size of the array (byte 0)
% % This value is funny - is the 'x' value for image, and 
% % data must be read in this way 
fseek(fid,656,'bof');
yCCD = fread(fid,1,'short');

%% Why can't I read type 'long' from 660-668 ?? %%

fseek(fid,42,'bof');
xCCD = fread(fid,1,'short');

% % check the data type (byte 108 in header)
fseek(fid,108,'bof');
dType = fread(fid,1,'short');

if dType == 0
   dTypeStr = 'float32';
elseif dType == 1
   dTypeStr = 'int32'; 
elseif dType == 3
   dTypeStr = 'uint16';
elseif dType == 6
   dTypeStr = 'uint8';
else 
   error('Data Type Not Yet Accounted For... Please add (check manual for type)');
end

if isNative, dTypeStr = ['*',dTypeStr]; end %% changes how fread gets data

% % go to beginning of data, and read appropriate number of bytes
fseek(fid,4100,'bof'); %% 4100 beginning of data
datRaw = fread(fid, [xCCD inf], dTypeStr);

fclose(fid); %% close file

% % format column vector into correct 
dat = datRaw'; clear datRaw;
[mRaw,nRaw] = size(dat);

if mRaw ~= yCCD %% making 3-D array if multiple images saved
	dat = reshape(dat,yCCD,xCCD,mRaw/yCCD);
end

[m,n,o] = size(dat);

if isPlot
   f1 = figure; 
   mess = 'Press any key to continue ... <ctrl+c to end> ';
   disp(mess)
   if ( m == 1 | n == 1 ) & o == 1
      plot(dat), title(mess), pause, else
      for ii = 1:o
      imagesc(dat(:,:,ii)); colorbar; title(mess), pause; end
   end

   delete(f1)
end

if nargout == 1
   varargout = {dat};
end
   