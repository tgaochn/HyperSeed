function info = readmtlfile(file)
%READMTLFILE Read MTL text file
%   INFO = readmtlfile(FILE) reads Landsat or EO satellite MTL text FILE
%   and returns a structure whose fields contain information about TIFF
%   file names exist in MTL file and their corresponding band wavelength,
%   gain and offset values and also the other information describing data
%   in TIFF files.
%
%   NOTE
%   ----
%   The input FILE must be a text file and contains MTL string in it's
%   filename.
%
%   See also hypercube, enviinfo.

%   Copyright 2020 The MathWorks, Inc.


% Read file
arguments
    file {mustBeStringOrChar}
end
file = matlab.images.internal.stringToChar(file);

[filepath,~,fileExt] = fileparts(file);

% Make sure that input file is MTL text file, MTL text file has .txt
% extension and filename includes the 'MTL' character array.
if ~(strcmpi(fileExt,'.txt') && contains(file,'MTL','IgnoreCase',true))
    error(message('hyperspectral:hypercube:invalidMTL'));
end

% Read file
[fid, msg] = fopen(file,'r');
if (fid < 0)
    error(msg);
end
str = fread(fid,'uint8=>char');
fclose(fid);

% Split into lines
lines = splitlines(str');

% Split lines at '=' and make a list of key-value pairs.
hasEq = contains(lines,'=');
keyValue = strtrim(split(lines(hasEq),'='));
% Remove '"' characters
keyValue = erase(keyValue,'"');

% Error out if the file don't contain any key value pairs
if(size(keyValue,2) ~= 2)
    error(message('hyperspectral:hypercube:invalidText'));
end
key = keyValue(:,1);
value = keyValue(:,2);

info = hyper.internal.getReqFieldNames();

% Get MTL file details
info = getFileDetails(file,info);
% Get spectral band details
info = getSpectralInfo(info,key, value, filepath);
% Get datacube details
info = getRasterInfo(info, key, value, numel(info.TIFFile));
% Get image parameters
info = getImgParam(info, key, value);

% Acquisition time
hasDate = contains(key, ["ACQUISITION_DATE","DATE_ACQUIRED"],'IgnoreCase',true);
if any(hasDate)
    info.AcquisitionTime = string(value(hasDate));
end
end


function info = getFileDetails(filename, info)

% Obtain information about the file.
fid = fopen(filename, 'r');
info.Filename = string(fopen(fid));
fclose("all");

d = dir(info.Filename);

info.FileModDate   = string(d.date);
info.FileSize      = d.bytes;
info.Format        = "txt";
info.FormatVersion = '';
end


function info = getSpectralInfo(info, key, value, filepath)

satID = string([]);
sensorID = string([]);

% Look for spacecraft and sensor IDs
satIdx = strcmpi(key,'SPACECRAFT_ID');
senIdx = strcmpi(key,'SENSOR_ID');

if any(satIdx) && any(senIdx)
    satID = value{satIdx};
    sensorID = value{senIdx};
    info.SensorType = strcat(satID," ",sensorID);
end

% Band 6 of Landsat ETM collects high and low gain for all scenes. So it
% has two bands names Band 6_VCID_1, Band 6_VCID_2. Erasing '_VCID_'
% character to get valid band number for these files.
key = erase(key,'_VCID_');

% Get TIFF files
% bandNum contains existed band numbers in MTL file
[bandNum, info] = getTifFiles(info, key, value, filepath);
% Assign radiance gain and offset values
info = getGainOffset(info, key, value, bandNum);
% Read solar irradiance values
solarIrradiance =  hyper.internal.getSolarIrradiance(info.SensorType);
% Intialise pan and thermal band numbers 
panBandNo = [];
thermalbandNo = [];

% Wavelength values are hard coded from USGS website. (https://www.usgs.gov/)
if(strncmpi(satID,'Landsat',7))
    if(strcmpi(sensorID,'OLI_TIRS'))
        % Landsat 8 (OLI_TIRS) satellite images contain 11 bands
        %  BandNum     BandType         Spatial Resolution
        %  ----------------------------------------------
        %  1-7 and 9   Multispectral         30
        %  8           Panchromatic          15
        %  10-11       Thermal               100
        
        % Define band names and center wavelengths for each band
        origBandNum = 1:11;
        bandNames = ["Coastal aerosol";"Blue";"Green";"Red";"NIR";...
            "SWIR 1";"SWIR 2";"Panchromatic";"Cirrus";"Thermal Infrared 1";"Thermal Infrared 2"];
        wavelength = [440;480;560;655;865;1610;2200;590;1370;1089.5;1200.5];
        
        % Assign wavelengths for existed bands
        info.BandNames = bandNames(bandNum);
        info.Wavelength = wavelength(bandNum);
        info.WavelengthUnits = "Nanometers";

        % Differentiate panchromatic, thermal bands and Multi spectral
        % bands
        panBandNo = 8;
        thermalbandNo = [10,11];
        info = extractNonMS(info, bandNum, panBandNo, thermalbandNo);        
        
    elseif(strcmpi(sensorID,'ETM'))
        % Landsat 7 (ETM+) satellite images contain 8 bands
        %  BandNum     BandType         Spatial Resolution
        %  ----------------------------------------------
        %  1-5 and 7   Multispectral         30
        %  6           Thermal(Low and High) 60
        %  8           Panchromatic          15
           
        % Define band names and center wavelengths for each band
        origBandNum = [1,2,3,4,5,61,62,7,8];
        bandNames = ["Band 1";"Band 2";"Band 3";"Band 4";"Band 5";"Band 6 VCID 1";"Band 6 VCID 2";"Band 7";"Band 8"];
        wavelength = [485;560;660;835;1650;11450;11450;2220;710];
        
        % Assign wavelengths
        [~, existedNumIdx] = ismember(bandNum, origBandNum);
        info.BandNames = bandNames(existedNumIdx);
        info.Wavelength = wavelength(existedNumIdx);
        info.WavelengthUnits = "Nanometers";
        
        % Differentiate panchromatic, thermal bands and Multi spectral
        % bands
        panBandNo = 8;
        thermalbandNo = [61,62];
        info = extractNonMS(info, bandNum, panBandNo, thermalbandNo);
        
    elseif(strcmpi(sensorID,'TM'))
        % Landsat 4-5 (TM) satellite images contain 7 bands
        %  BandNum     BandType         Spatial Resolution
        %  ----------------------------------------------
        %  1-5 and 7   Multispectral         30
        %  6           Thermal               120
        
        % Define band names and center wavelengths for each band
        origBandNum = 1:7;
        bandNames = ["Band 1";"Band 2";"Band 3";"Band 4";"Band 5";"Band 6";"Band 7"];
        wavelength = [485;560;660;830;1650;11450;2215];
        
        % Assign wavelengths
        info.BandNames = bandNames(bandNum);
        info.Wavelength = wavelength(bandNum);
        info.WavelengthUnits = "Nanometers";
        
        % Differentiate panchromatic, thermal bands and Multi spectral
        % bands
        panBandNo = [];
        thermalbandNo = 6;
        info = extractNonMS(info, bandNum, panBandNo, thermalbandNo);
        
    elseif(strcmpi(sensorID,'MSS'))
        % Landsat 1-5 (MSS) satellite images contain 4 multispectral bands
        % with 60m spatial resolution.
        
        % Define band names and center wavelengths for each band
        origBandNum = 1:4;
        landsatNo = str2double(regexp(satID,'[\d.]+','match'));
        % Band names are different for Landsat 1-3 and Landsat 4-5
        if(landsatNo == 4 || landsatNo == 5)
            bandNames = ["Band 1";"Band 2";"Band 3";"Band 4"];
        else
            bandNames = ["Band 4";"Band 5";"Band 6";"Band 7"];
            bandNum = bandNum-3;
        end
        wavelength = [550;650;750;950];
        info.WavelengthUnits = "Nanometers";
        
        % Assign wavelengths
        info.BandNames = bandNames(bandNum);
        info.Wavelength = wavelength(bandNum);
    end
elseif(strcmpi(satID,'EO1'))
    if(strcmpi(sensorID,'HYPERION'))
        % Hyperion level-1 data contains 242 bands (70 VNIR, 172 SWIR)
        
        origBandNum = 1:242;
        bandNames = strings(242,1);
        for i = 1:242
            bandNames(i) = sprintf('Band %d',i);
        end
        
        % Spectral resolution for VNIR and SWIR bands
        VNIR_BW = 10.09;
        SWIR_BW = 10.11;
        
        % Total number of VNIR and SWIR bands
        VNIR_Bands = 70;
        SWIR_Bands = 172;
        
        % Hyperion sensor wavelength range is 357nm to 2576nm
        VNIR_wl = 357+(0:VNIR_Bands-1)*VNIR_BW;
        SWIR_wl = fliplr(2576-(0:SWIR_Bands-1)*SWIR_BW);
        
        % Concatenate VNIR and SWIR wavelengths
        wavelength = [VNIR_wl, SWIR_wl]';
        
        % Define bad band multipliers for Hyperion sensor
        % 1-7, 225-242 are unused bands. 58-76 fall in overlap region,
        % 77,78 contains high noise.
        BadBands = true(242,1);
        noisyBands = [1:7,225:242,58:76,77,78]';
        BadBands(noisyBands) = false;
        
        info.BandNames = bandNames(bandNum);
        info.Wavelength = wavelength(bandNum);
        info.WavelengthUnits = "Nanometers";
        info.BadBands = BadBands(bandNum);
        
        gain = 1./[repelem(info.Gain(1),70,1);repelem(info.Gain(2),172,1)];
        info.Gain = gain(bandNum);
        
        % Simulated transmittance spectra of atmospheric Water vapor,
        % Ozone, Oxygen.
        [waterVapour, oxygen, ozone] = hyper.internal.getGaseousAbsorption(info.SensorType);
        info.WaterVapourAbsorption = waterVapour(bandNum);
        info.OxygenAbsorption = oxygen(bandNum);
        info.OzoneAbsorption = ozone(bandNum);
    elseif(strcmpi(sensorID, 'ALI'))
        % EO-1 (ALI) satellite images contain 10 bands
        %  BandNum     BandType         Spatial Resolution
        %  ----------------------------------------------
        %  1           Panchromatic          10
        %  2-10        Multispectral         30
        
        % Define band names and center wavelengths for each band
        origBandNum = 1:10;
        bandNames = ["Panchromatic";"MS-1p";"MS-1";"MS-2";"MS-3";"MS-4";"MS-4p";"MS-5p";"MS-5";"MS-7"];
        wavelength = [585;443;482.5;565;660;790;867.5;1250;1650;2215];
        
        % Assign wavelengths
        info.BandNames = bandNames(bandNum);
        info.Wavelength = wavelength(bandNum);
        info.WavelengthUnits = "Nanometers";
        
        % Differentiate panchromatic, thermal bands and Multi spectral
        % bands
        panBandNo = 1;
        thermalbandNo = [];
        info = extractNonMS(info, bandNum, panBandNo, thermalbandNo);
    end
end

% Extract existed multispectral indices and assign solar irradiance values
% for those bands only.
if ~isempty(solarIrradiance)
    allMutispectral = setdiff(origBandNum,thermalbandNo,'stable');
    existedMs = setdiff(bandNum,[panBandNo, thermalbandNo],'stable');
    [~, existedBandIdx] = ismember(existedMs,allMutispectral);
    info.SolarIrradiance = solarIrradiance(existedBandIdx);
    
    % Panchromatic 
    if any(ismember(bandNum,panBandNo))
        info.Panchromatic.SolarIrradiance = solarIrradiance(panBandNo);
    end        
end
end


function [bandNum, info] = getTifFiles(info, key, value, filepath)
% Get image file names from MTL file

% Initialize output parameters
bandNum = [];
info.TIFFile = "";

% Output format must be GeoTIFF
hasFormat = contains(key,'OUTPUT_FORMAT');
if ~any(hasFormat) || (any(hasFormat) && ~strcmpi(value{hasFormat},'GEOTIFF'))
    error(message('hyperspectral:hypercube:invalidFormat'));
end

% Image file contains 'file_name' and 'Band' strings
hasBand = contains(key,'BAND','IgnoreCase',true);
hasFile = contains(key,'FILE_NAME','IgnoreCase',true);

tifFileIdx = find(hasBand & hasFile);

if any(tifFileIdx)
    % Read Quality band, Quality band file contains 'BQA' character array in
    % it's file name.
    hasQualityBand = contains(value(tifFileIdx),'BQA','IgnoreCase',true);
    if any(hasQualityBand)
        QualityBand = value{tifFileIdx(hasQualityBand)};
        % Attach entire file path
        info.QualityBand = string(fullfile(filepath,QualityBand));
        % Remove Quality band from file indices as it is not a multi spectral file
        tifFileIdx(hasQualityBand) = [];
    end
    
    bandNum = str2double(string(regexp(key(tifFileIdx),'\d+','match')));
    info.TIFFile = fullfile(filepath,string(value(tifFileIdx)));
end
end


function info = getGainOffset(info, key, value, bandNum)

% Gain and Offset values present in either RADIANCE_SCALING or
% RADIOMETRIC_RESCALING group
hasScalingGroup = find(contains(value,["RADIANCE_SCALING","RADIOMETRIC_RESCALING"],'IgnoreCase',true));
if(hasScalingGroup)
    % Assign scaling group into temporary variables
    tempKey = key(hasScalingGroup(1)+1:hasScalingGroup(2)-1);
    tempVal = value(hasScalingGroup(1)+1:hasScalingGroup(2)-1);
    
    % Get indices of radiance gain values
    hasGain = contains(tempKey,["RADIANCE_MULT","SCALING_FACTOR"],'IgnoreCase',true);
    
    % If the sensor is Hyperion, MTL text file has two scaling factors for
    % VNIR filter and SWIR filter, otherwise scaling factor values defined
    % for each and every band.
    hasVNIR = contains(tempKey,"VNIR");
    hasSWIR = contains(tempKey,"SWIR");
    % VNIR scaling factors key word present in SCALING_FACTOR group and
    % contains VNIR string and it is same for SWIR.
    hasVNIRGain = find(hasGain & hasVNIR);
    hasSWIRGain = find(hasGain & hasSWIR);
    if any(hasVNIRGain) && any(hasSWIRGain)
        info.Gain = [str2double(tempVal(hasVNIRGain));str2double(tempVal(hasSWIRGain))];
        return;
    end
    
    % If the sensor is multispectral, scaling factor will be defined for
    % each and every band
    if any(hasGain)
        % Initialize gain
        info.Gain = nan(numel(bandNum),1);
        % Get band numbers from gain keys
        gainBandNum = str2double(string(regexp(tempKey(hasGain),'\d+','match')));
        % Get gain values
        gain = str2double(tempVal(hasGain));
        % Take out only gain values of existed bands
        [hasNum, order] = ismember(bandNum, gainBandNum);
        % Remove non-matched bands
        order(order == 0) = [];
        info.Gain(hasNum) = gain(order);
    end
    
    % Get indices of radiance offset values
    hasOffset =  contains(tempKey,["RADIANCE_ADD","OFFSET"],'IgnoreCase',true);
    
    if any(hasOffset)
        % Initialize offset
        info.Offset = nan(numel(bandNum),1);
        % Get band numbers from offset keys
        offsetBandNum = str2double(string(regexp(tempKey(hasOffset),'\d+','match')));
        % Get offset values
        offset = str2double(tempVal(hasOffset));
        % Take out only offset values of existed bands
        [hasNum, order] = ismember(bandNum, offsetBandNum);
        % Remove non-matched bands
        order(order == 0) = [];
        info.Offset(hasNum) = offset(order);
    end
    
    % Get indices of reflectance gain values
    hasRefGain = contains(tempKey,'REFLECTANCE_MULT','IgnoreCase',true);
    
    if(any(hasRefGain))
        % Initialize reflectance gain
        info.ReflectanceGain = nan(numel(bandNum),1);
        % Get band numbers from reflectance gain keys
        refGainBandNum = str2double(string(regexp(tempKey(hasRefGain),'\d+','match')));
        % Get reflectance gain values
        refGain = str2double(tempVal(hasRefGain));
        % Take out only reflectance gain values of existed bands
        [hasNum, order] = ismember(bandNum, refGainBandNum);
        % Remove non-matched bands
        order(order == 0) = [];
        info.ReflectanceGain(hasNum) = refGain(order);
    end
    
    % Get indices of reflectance offset values
    hasRefOffset = contains(tempKey,'REFLECTANCE_ADD','IgnoreCase',true);
    
    if(any(hasRefOffset))
        % Initialize reflectance offset
        info.ReflectanceOffset = nan(numel(bandNum),1);
        % Get band numbers from reflectance offset keys
        refOffsetBandNum = str2double(string(regexp(tempKey(hasRefOffset),'\d+','match')));
        % Get reflectance offset values
        refOffset = str2double(tempVal(hasRefOffset));
        % Take out only reflectance offset values of existed bands
        [hasNum, order] = ismember(bandNum, refOffsetBandNum);
        % Remove non-matched bands
        order(order == 0) = [];
        info.ReflectanceOffset(hasNum) = refOffset(order);
    end
end
end


function info = extractNonMS(info, bandNum, panBandNo, thermalbandNo)

% Separate thermal band from multispectral bands
fields = ["TIFFile";"BandNames";"Wavelength";"Gain";"Offset";"ReflectanceGain";"ReflectanceOffset"];
if ~isempty(thermalbandNo)
    hasThermal = ismember(bandNum, thermalbandNo);
    if any(hasThermal)
        for fieldNo = 1:numel(fields)-2
            if ~isempty(info.(fields(fieldNo)))
                info.ThermalBand.(fields(fieldNo)) = info.(fields(fieldNo))(hasThermal);
                info.(fields(fieldNo))(hasThermal) = [];
            end
        end
        % Thermal band don't contain reflectance gain and offsets
        if ~isempty(info.ReflectanceGain)
            info.ReflectanceGain(hasThermal) = [];
        end
        if ~isempty(info.ReflectanceOffset)
            info.ReflectanceOffset(hasThermal) = [];
        end
        
        bandNum(hasThermal) = [];
    end
end

% Separate panchromatic band from multispectral bands
if ~isempty(panBandNo)
    hasPan = ismember(bandNum, panBandNo);
    if any(hasPan)
        for fieldNo = 1:numel(fields)
            if ~isempty(info.(fields(fieldNo)))
                info.Panchromatic.(fields(fieldNo)) = info.(fields(fieldNo))(hasPan);
                info.(fields(fieldNo))(hasPan) = [];
            end
        end
    end
end
end


function info = getRasterInfo(info, key, value, numOfBands)

info.RasterFormat = "tiff";

hasHeight = contains(key,["PRODUCT_LINES_REF","REFLECTIVE_LINES"]);
if ~any(hasHeight)
    hasHeight = contains(key,"PRODUCT_LINES");
end
info.Height = str2double(value(hasHeight));

hasWidth =  contains(key,["REFLECTIVE_SAMPLES","PRODUCT_SAMPLES_REF"]);
if ~any(hasWidth)
    hasWidth = contains(key,"PRODUCT_SAMPLES");
end
info.Width = str2double(value(hasWidth));

info.Bands = numOfBands;
info.Interleave = "bsq";
end


function info = getImgParam(info, key, value)

% Image parameters exist in either PRODUCT_PARAMETER or
% IMAGE_ATTRIBUTES group.
hasImgParam = find(~cellfun(@isempty,(regexpi(value,'PRODUCT_PARAMETER|IMAGE_ATTRIBUTES'))));

if(~isempty(hasImgParam))
    % Assign image parameters group into temporary variables
    tempKey = key(hasImgParam(1)+1:hasImgParam(2)-1);
    tempVal = value(hasImgParam(1)+1:hasImgParam(2)-1);
    
    fields = {'CLOUD_COVER', 'SUN_AZIMUTH', 'SUN_ELEVATION', 'EARTH_SUN_DISTANCE', ...
        'CLOUD_COVER_LAND', 'SENSOR_LOOK_ANGLE'};
    actualNames = {'CloudCover','SunAzimuth','SunElevation','EarthSunDistance',...
        'CloudCoverLand','SensorLookAngle'};
    
    for fieldNo = 1:numel(fields)
        hasField = strcmpi(tempKey,fields{fieldNo});
        
        if(any(hasField))
            info.(actualNames{fieldNo}) = str2double(tempVal(hasField));
        end
    end
end
end


% Validator
function mustBeStringOrChar(file)
validateattributes(file, {'char','string'}, {'nonempty','scalartext'},'','FILE');
end
