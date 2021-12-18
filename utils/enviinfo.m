function hdrInfo = enviinfo(file)
%ENVIINFO Read metadata from ENVI header file
%   INFO = enviinfo(FILE) reads ENVI header file and returns a structure
%   whose fields contain information about ENVI data file and other
%   information related to data in ENVI data file. FILE must be a string
%   scalar or character vector that specifies the name of the ENVI
%   header(.hdr) file. If the file extension is not specified, this
%   function will look for a '.hdr' file.
%
%   INFO must contain these following fields information:
%   
%     Height       - An integer indicating the height of the image or 
%                    number of rows present in ENVI data file.                        
%
%     Width        - An integer indicating the width of the image or number 
%                    of columns present in ENVI data file.                       
%
%     Bands        - Numeric value specifying the number of bands present 
%                    in ENVI data file.
%
%     DataType     - String specifying the type of data representation in
%                    image file. Supported data types are:
%                        1  - "uint8" (8-bit unsigned integer)
%                        2  - "int16" (16-bit signed integer)
%                        3  - "int32" (32-bit signed integer)
%                        4  - "single" (32-bit single-precision)
%                        5  - "double" (64-bit double-precision floating-point)
%                        12 - "uint16" (16-bit unsigned integer)
%                        13 - "uint32" (32-bit unsigned integer)
%                        14 - "int64" (64-bit signed integer)
%                        15 - "uint64" (64-bit unsigned integer)
%                     FILE must contain any one of the above integer as
%                     data type.
%
%     Interleave   - String indicating data interleave. It must be "bsq" or
%                    "bil" or "bip".
%
%     HeaderOffset - Numeric value indicating the zero-based location of 
%                    the first data element in the image file. This value
%                    represents number of bytes from the beginning of the
%                    image file to where the data begins.
%
%     ByteOrder    - String specifying the byte ordering (machine format)
%                    in which the data is stored in image file. Supported
%                    byte order values are:
%                        0 for Little-endian("ieee-le")
%                        1 for Big-endian("ieee-be")
%
%   Example 1
%   ---------
%   % Read ENVI header file
%   info = enviinfo('paviaU.hdr');
%   hcube = hypercube(info.Filename);
%
%   Example 2
%   ---------
%   % Read ENVI image and header files
%   info = enviinfo('indian_pines');
%   data = multibandread('indian_pines.dat',[info.Height, info.Width, info.Bands],...
%          info.DataType, info.HeaderOffset, info.Interleave, info.ByteOrder);
%
%   See also hypercube, multibandread.

%   Copyright 2020 The MathWorks, Inc.


% Read file
arguments
    file {mustBeStringOrChar}
end

% For compiler tests
if isdeployed
    rootDir = ctfroot;
else
    rootDir = matlabroot;
end
% Register resources/hyperspectral.
matlab.internal.msgcat.setAdditionalResourceLocation(rootDir);

file = matlab.images.internal.stringToChar(file);
file = findActualFileName(file);

% Make sure that input file is ENVI header file. 
if ~hyper.internal.isenvi(file)
    error(message('hyperspectral:enviinfo:invalidHeaderFile'));
end

fid = fopen(file, 'r');
str = fread(fid,'uint8=>char');
fclose(fid);

% Find indices of curly braces '{','}' and Remove End of Line characters in
% between curly braces
startIdx = find(str == 123);
endIdx = find(str == 125);
% Number of open and close curly braces should be equal.
tf = ~isequal(numel(startIdx),numel(endIdx)) || any(~(endIdx>startIdx));
if(~tf)
    numBraces = numel(startIdx);
    for braceNo = 1:numBraces
        % [13,10] are ASCII values of End of Line character(\r\n)
        % Find indices of End of Line characters in between curly braces
        eolIdx = find(ismember(str(startIdx(braceNo):endIdx(braceNo)),[10,13]))+startIdx(braceNo)-1;
        % Replace with space
        str(eolIdx) = 32;  % 32 is ASCII value for space
    end
end
% Replace curly braces with spaces
str(startIdx) = 32;
str(endIdx) = 32;

% Split string into lines
lines = splitlines(str');
% Remove inline comments (Key-value pairs may have comments appended to
% them by preceding the comment with a semicolon <;>. Third-party tools can
% ignore all characters including and following a semicolon <;> to the end
% of line code. Where no key is stated, for example when an ASCII line
% starts with a semicolon, the whole line is treated as comment.)
hasComment = strfind(lines,';'); 
commentIdx = find(~cellfun(@isempty,hasComment));
for idx = 1:numel(commentIdx)
    % Remove all characters until end of line
    lines{commentIdx(idx)}(hasComment{commentIdx(idx)}:end) = [];
end

if ~any(~cellfun(@isempty,lines(2:end)))
    % Throw error if file is empty
    error(message('hyperspectral:enviinfo:emptyHeader'));
end

% Split lines into list of key-value pairs
[key, value] = getKeyValuePair(lines);
% Remove empty key-value pairs
hasNonEmptyKeys = ~cellfun(@isempty,key);

key = key(hasNonEmptyKeys);
value = value(hasNonEmptyKeys);
% Convert cell array to structure
tempHeader = cell2struct(value',key',2);
% Get all required field names
hdrInfo = hyper.internal.getReqFieldNames();
hdrInfo.DataType = [];
hdrInfo.ByteOrder = [];

% File details
hdrInfo = getFileDetails(file, hdrInfo);
% Combine structures
for fieldNo = 1:numel(fieldnames(tempHeader))
    hdrInfo.(key{fieldNo}) = tempHeader.(key{fieldNo});
end

% Get geographic map information
hasMap = strcmpi(key,'MapInfo');
if(any(hasMap))
    hdrInfo.MapInfo = getmap(value{hasMap});
end

% Validate required fields of metadata
hdrInfo = validateMetadata(hdrInfo);
hdrInfo.DataType = getDataType(hdrInfo.DataType);
hdrInfo.ByteOrder = getByteOrder(hdrInfo.ByteOrder);
hdrInfo.RasterFormat = "ENVI";
end


function file = findActualFileName(origFile)

[~,~,fileExt] = fileparts(origFile);
baseFileName = erase(origFile,fileExt);

if isequal(exist([baseFileName, '.hdr'],'file'),2)
    file = [baseFileName, '.hdr'];
elseif isequal(exist([origFile, '.hdr'],'file'),2)
    file = [origFile, '.hdr'];
else
    error(message('hyperspectral:hypercube:invalidFile', [baseFileName, '.hdr']));
end
end


function [key, value] = getKeyValuePair(lines)

% Initialize key-value pair
key = cell(numel(lines)-1,1);
value = cell(numel(lines)-1,1);

% lineNo starts with 2 as first line contains 'ENVI' key word
for lineNo = 2:numel(lines)
    % Split lines at first equal
    eqIdx = find(ismember(lines{lineNo},'='),1);
    if isempty(eqIdx) && any(isletter(lines{lineNo}))
        error(message('hyperspectral:enviinfo:noSeparator', num2str(lineNo)));
    end
    
    tempKey = lines{lineNo}(1:eqIdx-1);
    tempVal = strtrim(lines{lineNo}(eqIdx+1:end));
    % Convert tempKey character array to valid structure field name and
    % format the value datatype
    [tempKey, tempVal] = modifyKeyValue(tempKey,tempVal);
    
    key{lineNo-1} = tempKey;
    value{lineNo-1} = tempVal;
end
end


function [key,value] = modifyKeyValue(key, value)

% Convert ENVI key string to structure field
key = regexprep(key,'(\<[a-z])','${upper($1)}');
key(isspace(key)) = [];

value = string(value);
switch key
    case {'AcquisitionTime','BandNames','ClassNames','ReadProcedures','SpectraNames','ZPlotTitles'}
        value = string(strtrim(split(value,',')));              % String array 
    case {'Bands','ByteOrder','Classes','CloudCover','ComplexFunction',...
            'DataIgnoreValue','DataType','DefaultBands','DemBand','HeaderOffset',...
            'SolarIrradiance','SunAzimuth','SunElevation','Wavelength','XStart',...
            'YStart'}
        value = str2double(strtrim(split(value,',')));          % Numeric array
    case 'Bbl'
        key = 'BadBands';                                       % Logical multiplier for each band
        value = logical(str2double(strtrim(split(value,','))));
    case 'ClassLookup'                                          % List of class colors (R,G,B)
        value = reshape(str2double(strtrim(split(value,','))),3,[])';
    case 'DataGainValues'
        key = 'Gain';                                           % Gain values for each band
        value = str2double(strtrim(split(value,',')));
    case 'DataOffsetValues'
        key = 'Offset';                                         % Offset values for each band
        value = str2double(strtrim(split(value,',')));
    case 'DataReflectanceGainValues'
        key = 'ReflectanceGain';                                % Reflectance gain values for each band
        value = str2double(strtrim(split(value,',')));
    case 'DataReflectanceOffsetValues'
        key = 'ReflectanceOffset';                              % Reflectance offset values for each band
        value = str2double(strtrim(split(value,',')));
    case 'Fwhm'
        key = 'FWHM';                                           % Numeric value defining FWHM for each band
        value = str2double(strtrim(split(value,',')));
    case 'GeoPoints'                                            % Location (x, y, latitude, longitude)
        value = reshape(str2double(strtrim(split(value,','))),4,[])';
    case 'Lines'
        key = 'Height';                                         % Scalar numeric value
        value = str2double(strtrim(split(value,',')));          
    case {'PixelSize','ZPlotAverage','ZPlotRange'}              % [x,y] values
        value = str2double(strtrim(split(value,',')))';
    case 'RpcInfo'
        key = 'RPCInfo';                                             
        value = string(strtrim(split(value,',')));
    case 'Samples'
        key = 'Width';                                          % Scalar numeric value
        value = str2double(strtrim(split(value,','))); 
    otherwise
        [tempVal, isNum] = str2num(value);                      % All other fields
        if isNum
            value = reshape(tempVal,[],1);;
        end
end
end


function hdrInfo = getFileDetails(filename, hdrInfo)

% Obtain information about the file.
fid = fopen(filename, 'r');
hdrInfo.Filename = string(fopen(fid));
fclose("all");

d = dir(hdrInfo.Filename);

hdrInfo.FileModDate   = string(d.date);
hdrInfo.FileSize      = d.bytes;
hdrInfo.Format        = "HDR";
hdrInfo.FormatVersion = '';
end


function data = getmap(value)

data = [];
% Split at ','
value = string(strtrim(split(value,',')));
try
    % First value represent projection
    data.ProjType = value(1);
    % Later two values contain pixel tie points
    data.PixelTiePoints = [str2double(value(2)),str2double(value(3))];
    % Later two values contain map tie points
    data.MapTiePoints = [str2double(value(4)),str2double(value(5))];
    % Later two values contain pixel size
    data.PixelSize = [str2double(value(6)),str2double(value(7))];
    
    count = 8;
    % If the projection is UTM, later points contain projection zone and
    % direction, otherwise it will be datum and units
    if(strcmpi(data.ProjType,'UTM'))
        data.ProjZone = str2double(value(8));
        data.NorthOrSouth = value(9);
        count = 10;
    end
    % Datum
    data.Datum = value(count);
    
    % Units (units value can be 'units = meter' or 'meter')
    if contains(value(count+1),'=')
        % If the value is in the form of 'units = meter'
        units = strtrim(split(value(count+1),'='));
        if(strcmpi(units(1),'units'))
            data.Units = units(2);
        end
    else
        % If the value just contains units (ex: 'meter')
        data.Units = value(count+1);
    end
catch ME
    % Ignore error regarding exceeding index.
    if ~strcmpi(ME.identifier,'MATLAB:badsubscript')
        % Throw error as it is, If error is not related to exceeding index.
        rethrow(ME);
    end
end
end


function datatype = getDataType(datatype)

switch datatype
    case 1
        datatype = "uint8";
    case 2
        datatype = "int16";
    case 3
        datatype = "int32";
    case 4
        datatype = "single";
    case 5
        datatype = "double";
    case 12
        datatype = "uint16";
    case 13
        datatype = "uint32";
    case 14
        datatype = "int64";
    case 15
        datatype = "uint64";
    otherwise
        error(message('hyperspectral:enviinfo:unknownDataType'));
end
end


function byteorder = getByteOrder(byteOrder)

if(byteOrder == 1)
    byteorder = "ieee-be";
elseif (byteOrder == 0)
    byteorder = "ieee-le";
else
    error(message('hyperspectral:enviinfo:unknownByteOrder'));
end
end


function hdrInfo = validateMetadata(hdrInfo)
% Validate required fields

% Required fields
requiredFields = {'Height', 'Width', 'Bands', 'HeaderOffset',...
    'DataType', 'ByteOrder', 'Interleave'};

for fieldNo = 1 : numel(requiredFields)    
    if ~(fieldNo == 7)    % Other than interleave which is string
        validateattributes(hdrInfo.(requiredFields{fieldNo}), {'numeric'}, ...
            {'nonempty','finite','real','nonnegative','integer','scalar'},'enviinfo',[requiredFields{fieldNo} ' in input header file']);
    end
end

% Interleave
validateattributes(char(hdrInfo.Interleave), {'string','char'}, ...
    {'nonempty','scalartext'},'enviinfo','Interleave in input header file');

hdrInfo.Interleave = validatestring(hdrInfo.Interleave, ["bsq","bip","bil"],'enviinfo','Interleave in input header file');
end


% Validator
function mustBeStringOrChar(file)
validateattributes(file, {'char','string'}, {'nonempty','scalartext'},'','FILE');
end
