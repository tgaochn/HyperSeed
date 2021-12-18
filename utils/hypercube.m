classdef hypercube
    %HYPERCUBE Read hyperspectral imagery and metadata
    %   The hypercube class reads the hyperspectral image file and
    %   constructs a hypercube object with properties describing
    %   hyperspectral image and methods to manipulate the data it contains.
    %
    %   HCUBE = hypercube(FILE) reads input FILE and constructs a hypercube
    %   object HCUBE. FILE is a string scalar or character vector that
    %   specifies the name of the NITF file or Landsat/EO MTL text file or
    %   ENVI image file or ENVI header file or Hyperion L1R dataset file.
    %   If the input FILE is either of ENVI image or header file, both
    %   files must have same base file name and exist in specified file
    %   path. The input ENVI filename can include either the .hdr or .img
    %   or .dat extensions or just the base filename. If the input FILE is
    %   L1R image or header file, both files must have same base file name
    %   and present in specified path. The Hyperion L1R dataset has
    %   extension either .L1R or .hdr. The input FILE must contain
    %   wavelength information.
    %
    %   HCUBE = hypercube(IMGFILE, HDRFILE) reads band interleaved data
    %   from ENVI binary image file IMGFILE or Hyperion L1R dataset file
    %   using metadata from ENVI header file HDRFILE and creates hypercube
    %   object HCUBE. HDRFILE must contain dimension of DataCube,
    %   interleave, byte order, data type, header offset, and wavelength
    %   information.
    %
    %   HCUBE = hypercube(___, WAVELENGTH) reads input file and updates
    %   wavelength with given WAVELENGTH vector. WAVELENGTH is C-by-1
    %   numeric vector which represents the wavelength of each band in
    %   'DataCube' respectively. C is the number of spectral bands in input
    %   image file. The wavelengths must be in nanometer.
    %
    %   HCUBE = hypercube(FILE, WAVELENGTH) reads image from TIFF format
    %   FILE and updates the WAVELENGTH for each band.
    %
    %   HCUBE = hypercube(IMGARRAY, WAVELENGTH) creates a hypercube object
    %   with a hyperspectral datacube IMGARRAY with their respective
    %   wavelengths WAVELENGTH in nanometer. IMGARRAY is M-by-N-by-C
    %   numeric matrix where C is number of spectral bands. WAVELENGTH is a
    %   C-by-1 numeric vector.
    %
    %   HCUBE = hypercube(IMGARRAY, WAVELENGTH, METADATA) creates hypercube
    %   object with specified three properties information. METADATA is a
    %   structure containing IMGARRAY information.
    %
    %   SUPPORTED FILE FORMATS
    %   ----------------------
    %   hypercube supports following file formats:
    %
    %   ENVI format (binary image file & header file(.hdr))
    %   TIFF format
    %   NITF: Version 2.0, 2.1 and NSIF 1.0
    %
    %   hypercube supports following satellite data MTL text files
    %
    %   Earth Observing(EO)-1 (Hyperion & ALI)
    %   Landsat datasets(MSS, TM, ETM+, OLI/TIRS)
    %
    %   hypercube properties:
    %         DataCube      -  Hyperspectral image array, M-by-N-by-C numeric Array.
    %         Wavelength    -  Center wavelengths, C-by-1 double array.
    %         Metadata      -  Structure containing information about DataCube
    %
    %   hypercube methods:
    %         removeBands        - Remove specified spectral bands
    %         assignData         - Assign hypercube with new values
    %         cropData           - Crop DataCube
    %         enviwrite          - Write hypercube object as ENVI files
    %         selectBands        - Select most informative bands from DataCube
    %         colorize           - Estimate colored image from DataCube
    %
    %   NOTE
    %   ----
    %   When reading from ENVI files, both the image and header files must
    %   be available. If only one file is given as an input, the other file
    %   must be in the same directory and have the same filename except for
    %   the extension.
    %
    %   Example 1
    %   ---------
    %   % Load Indian Pines data collected from AVIRIS sensor
    %   load('indian_pines.mat');
    %
    %   % Construct hypercube object with Indian Pines data and wavelengths
    %   hcube = hypercube(indian_pines, wavelengths);
    %
    %   % Visualize RGB image
    %   [RGBImg, bandNum] = colorize(hcube,'method','rgb');
    %   figure, imshow(RGBImg)
    %
    %   % Remove bands affected by water vapor and write corrected data to
    %   % new file
    %   newhcube = removeBands(hcube, 'BandNumber', [104:108,150:163,220]);
    %   % Apply normalization to corrected data
    %   data = rescale(newhcube.DataCube);
    %   % Assign with normalized data
    %   normalizedCube = assignData(newhcube, ':',':',':',data);
    %   % Write normalized corrected data to new ENVI file
    %   enviwrite(normalizedCube, 'normalized_data');
    %
    %   % Find most informative bands and write them to new ENVI image and
    %   % header format file
    %   sig = fippi(hcube, 5);
    %   [newhcube, bandNum] = selectBands(hcube, sig);
    %   enviwrite(newhcube, 'informativeBands');
    %
    %   % Visualize false colored image
    %   [falseColoredImg, bandNum] = colorize(newhcube, [1,2,3]);
    %   figure, imshow(falseColoredImg)
    %
    %   See also multibandread, multibandwrite, nitfread, hdfread, ppi,
    %   enviinfo, fippi, nfindr, estimateAbundanceLS, ndvi, hyperpca,
    %   hypermnf, inverseProjection, countEndmembersHFC.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    
    properties (SetAccess = private, GetAccess = public)
        %DataCube - Hyperspectral image array
        %   DataCube is a M-by-N-by-C numeric array where M and N represent
        %   two spatial dimensions and C represent number of bands
        %   (spectral dimension).
        %   This is a read-only property.
        DataCube (:,:,:) {mustBeNumeric, mustBeNonsparse, mustBeReal}
        
        %Wavelength - Center wavelengths
        %   C-by-1 double vector containing center wavelength values of each
        %   spectral band.
        %   This is a read-only property.
        Wavelength (:,1) double {mustBeNumeric, mustBeNonsparse, mustBeFinite, mustBeReal, mustBePositive}

        %Metadata - Structure containing information about DataCube
        %   Metadata is a structure with fields containing hyperspectral
        %   image properties and spectral information. The set of fields in
        %   Metadata structure depends on the individual file and its
        %   format. The common fields are:
        %     Filename           - String containing the name of the input
        %                          file.
        %
        %     FileModDate        - String containing the modification date
        %                          of the input file.
        %
        %     FileSize           - Integer indicating the size of the file in bytes.
        %
        %     Format             - String containing the file format(HDR/TXT/NTF/TIFF).
        %
        %     FormatVersion      - Character vector or number specifying the file
        %                          format version.
        %
        %     SensorType         - Instrument type, such as Landsat TM, EO Hyperion, and so on.
        %
        %     Description        - String describing the image or the
        %                          processing performed.
        %
        %     AcquisitionTime    - String containing data acquisition time.
        %
        %     RasterFormat       - String containing raster file format:
        %                          ENVI files    - "ENVI"
        %                          NITF file     - "NTF"
        %                          MTL text file - "TIFF"
        %
        %     Height             - An integer indicating the height of the
        %                          image in pixels.
        %
        %     Width              - An integer indicating the width of the
        %                          image in pixels.
        %
        %     Bands              - Numeric value specifying the number of
        %                          bands in hyperspectral data.
        %
        %     DataType           - String specifying the type of data representation.
        %                          For example, "uint8","double","single","int8".
        %
        %     Interleave         - String indicating data interleave.
        %                          It must be "bsq" or "bil" or "bip".
        %
        %     HeaderOffset       - Numeric value indicating the zero-based location
        %                          of the first data element in the file. This value
        %                          represents number of bytes from the beginning of the
        %                          image file to where the data begins.
        %
        %     ByteOrder          - String specifying the byte ordering (machine format) in
        %                          which the data is stored in file: "ieee-le" or "ieee-be".
        %
        %     BandNames          - C-by-1 string array containing 'C' image band names.
        %
        %     FWHM               - C-by-1 double array containing full-width-half-maximum(FWHM)
        %                          values of each band in an image.
        %
        %     Gain               - C-by-1 double array of gain values for each band.
        %
        %     Offset             - C-by-1 double array of offset values for each band.
        %
        %     ReflectanceGain    - C-by-1 double array of reflectance gain values.
        %
        %     ReflectanceOffset  - C-by-1 double array of reflectance offset values.
        %
        %     BadBands           - Bad Band List, C-by-1 logical array containing
        %                          bad band multiplier values of each band in image.
        %                          Logical 0 for bad bands and logical 1 for good bands.
        %
        %     CloudCover         - Percentage of cloud cover within the image.
        %
        %     SunAzimuth         - Angle of the sun (in degrees) from due north
        %                          in a clockwise direction.
        %
        %     SunElevation       - Angle of the sun (in degrees) above the horizon.
        % 
        %     SolarIrradiance    - Top of the atmosphere solar irradiance per band.
        %                          Units are W/(m2 * um)
        %
        %     EarthSunDistance   - Distance between earth and sun in astronomical units
        % 
        %     WavelengthUnits    - Units of wavelength: Nanometers, Micrometers 
        %   This is a read-only property.
        Metadata (1,1) struct
    end
    
    methods
        %% Constructor
        function obj = hypercube(img, hdrfile, wavelength)
            
            arguments
                img {mustBeCharOrNum, mustBeNonempty}
                hdrfile {mustBeCharOrNum} = []
                wavelength = []
            end
            
            % For compiler tests
            if isdeployed
                rootDir = ctfroot;
            else
                rootDir = matlabroot;
            end
            % Register resources/hyperspectral.
            matlab.internal.msgcat.setAdditionalResourceLocation(rootDir);
            
            [imgfile, hdrfile, imgArray, wavelength, imgInfo] = parseInputs(img, hdrfile, wavelength);
            [~, ~, imgext] = fileparts(imgfile);
            
            % ENVI file contains one of these extensions:'.img','.hdr','.dat'
            % L1R dataset contains either .L1R or .hdr extension
            isEnviOrL1RFile = ~isempty(imgext) && any(strcmpi(imgext, {'.img','.dat','.hdr','.L1R'}));
            isBaseFile = isempty(imgext) && exist([imgfile,'.hdr'], 'file');
            
            if ~isempty(imgArray)
                % Numeric array
                obj.DataCube = imgArray;
                
                % Required fields for numeric array
                obj.Metadata = struct('Height',[],...
                    'Width',[],...
                    'Bands',[],...
                    'DataType',string([]),...
                    'Interleave',string([]),...
                    'HeaderOffset',[]);
                
                % Fill with known information
                [obj.Metadata.Height, obj.Metadata.Width, obj.Metadata.Bands] = size(imgArray);
                % 'DataCube' is in band sequential order, first two
                % dimensions represent spatial and third dimension
                % represent spectral.
                obj.Metadata.Interleave = "bsq";
                % Default is class of DataCube
                obj.Metadata.DataType = string(class(imgArray));
                obj.Metadata.HeaderOffset = 0;      % Default is 0
                
                % Assign user metadata if given
                if ~isempty(imgInfo)
                    % Remove file info, raster format, as hypercube is not
                    % reading file.
                    % Remove DataCube related fields from input metadata as
                    % they are already estimated.
                    % Second input  argument represents wavelength, so
                    % remove wavelength field from input metadata
                    allFields = fieldnames(imgInfo);
                    existedFields = ["Filename","FileModDate","FileSize","Format","FormatVersion",...
                        "Height","Width","Bands","DataType","Interleave","HeaderOffset",....
                        "Wavelength","WavelengthUnits","RasterFormat"];
                    allFields = allFields(~ismember(allFields, existedFields));
                    % Merge structures
                    for fieldNo = 1:numel(allFields)
                        obj.Metadata.(allFields{fieldNo}) = imgInfo.(allFields{fieldNo});
                    end
                    validateMetadata(obj.Metadata);
                end
                
            elseif isEnviOrL1RFile || isBaseFile
                % ENVI format or Hyperion L1R dataset
                
                if ~isempty(hdrfile)
                    if isequal(imgext, '.hdr')
                        % Swap image and header files, if in case first input
                        % file is header and second input file is image
                        tempFile = hdrfile;
                        hdrfile = imgfile;
                        imgfile = tempFile;
                    end
                else
                    hdrfile = imgfile;
                end
                
                % Read header
                obj.Metadata = getMetaforImg(hdrfile, 'envi');
                validateMetadata(obj.Metadata);
                
                imgfile = findActualImgFileName(imgfile);
                [~, ~, imgext] = fileparts(imgfile);
                
                if isequal(imgext, '.L1R')
                    validateImageFileFormat(imgfile);
                    % Extra bytes containing L1R header information
                    l1rOffset = 4840;
                    validateImageFile(imgfile, obj.Metadata,l1rOffset);
                    obj = readL1RFile(obj, imgfile);
                else
                    l1rOffset = 0;
                    validateImageFile(imgfile, obj.Metadata,l1rOffset);
                    obj = enviread(obj, imgfile);
                end
                
            elseif isnitf(imgfile)
                % Read first image from the NITF 2.0 or 2.1 file
                
                % Read Metadata
                obj.Metadata = getMetaforImg(imgfile, 'ntf');
                % Read DataCube
                obj.DataCube = nitfread(imgfile);
                % Update number of Bands and DataType
                obj.Metadata.Bands = size(obj.DataCube,3);
                obj.Metadata.DataType = string(class(obj.DataCube));
                
            elseif isTIFF(imgfile)
                % TIFF format
                
                % Read Metadata
                obj.Metadata = getMetaforImg(imgfile, 'tif');
                % Read DataCube
                obj.DataCube = imread(imgfile);
                % Update number of Bands and DataType
                obj.Metadata.Bands = size(obj.DataCube,3);
                obj.Metadata.DataType = string(class(obj.DataCube));
                
            elseif strcmpi(imgext,'.txt')
                % Satellite dataset (MTL file)
                
                % Read Metadata
                obj.Metadata = getMetaforImg(imgfile, 'mtl');
                
                obj.Metadata.SensorType = validateSensorType(obj.Metadata.SensorType);
                
                if isempty(obj.Metadata.Bands)
                    obj.DataCube = [];
                else
                    % Read DataCube
                    if obj.Metadata.Bands >= 1
                        try
                            cubeInCell = cell(obj.Metadata.Bands,1);
                            
                            for bandNum = 1 : obj.Metadata.Bands
                                cubeInCell{bandNum} = imread(obj.Metadata.TIFFile{bandNum},'tif');
                                if ~ismatrix(cubeInCell{bandNum})
                                    error(message('hyperspectral:hypercube:mustBeGrayImage', obj.Metadata.TIFFile(bandNum)));
                                end
                                if isempty(cubeInCell{bandNum})
                                    error(message('hyperspectral:hypercube:emptyFile', obj.Metadata.TIFFile(bandNum)));
                                end
                            end
                            obj.DataCube = cat(3,cubeInCell{:});
                        catch ME
                            switch ME.identifier
                                case 'MATLAB:catenate:dimensionMismatch'
                                    throwAsCaller(ME);
                                otherwise
                                    if ~any(strfind(ME.message, obj.Metadata.TIFFile(bandNum)))
                                        msg = getString(message('hyperspectral:hypercube:getFileName',obj.Metadata.TIFFile(bandNum)));
                                        causeException = MException(ME.identifier,msg);
                                        ME = addCause(ME,causeException);
                                    end
                                    throwAsCaller(ME);
                            end
                        end
                        
                        % Update image information
                        obj.Metadata.Height = size(obj.DataCube,1);
                        obj.Metadata.Width = size(obj.DataCube,2);
                        
                        % Update DataType, ByteOrder
                        obj.Metadata.DataType = string(class(obj.DataCube));
                        tempHeader = imfinfo(obj.Metadata.TIFFile{1});
                        
                        if(strcmp(tempHeader.ByteOrder,'big-endian'))
                            obj.Metadata.ByteOrder = "ieee-be";
                        else
                            obj.Metadata.ByteOrder = "ieee-le";
                        end
                    end
                end
            else
                error(message('hyperspectral:hypercube:unknownImageFormat'));
            end
            
            % Wavelength
            if isfield(obj.Metadata, 'Wavelength')
                % If the wavelength is in Metadata property, Store it in
                % Wavelength property
                obj.Wavelength = obj.Metadata.Wavelength;
                obj.Metadata = rmfield(obj.Metadata, 'Wavelength');
            end
            
            if ~isempty(wavelength)
                % If wavelength is provided as input, update existing
                % wavelengths with given wavelengths
                obj.Wavelength = wavelength;
            end
            
            if isempty(obj.Wavelength)
                % If the wavelength is not present in file and not provided
                % as input
                error(message('hyperspectral:hypercube:mustHaveWavelength'));
            end
            
            if ~(numel(unique(obj.Wavelength)) == size(obj.DataCube, 3))
                error(message('hyperspectral:hypercube:mustBeOfSize',size(obj.DataCube, 3),1));
            end
            
            % Assign wavelength units (Visible to thermalIR (200nm -
            % 100,000nm) or (0.2um to 100um)
            if ~any(~(obj.Wavelength > 200 & obj.Wavelength < 10^5))
                obj.Metadata.WavelengthUnits = "Nanometers";
            elseif ~any(~(obj.Wavelength > 0.2 & obj.Wavelength < 100))
                obj.Metadata.WavelengthUnits = "Micrometers";
            else
                obj.Metadata.WavelengthUnits = "";
            end
        end
        
        %% Array manipulation
        function obj = assignData(obj, row, col, band, newValues)
            %assignData Assign hypercube with new values
            %   NEWHCUBE = assignData(HCUBE, ROW, COLUMN, BAND, NEWVALUES)
            %   assigns NEWVALUES into the specified elements of 'DataCube'
            %   in hypercube object HCUBE and returns a NEWHCUBE having
            %   assigned data.
            %
            %   NEWVALUES must be a scalar value or a 3-D array of size
            %   numel(ROW)-by-numel(COLUMN)-by-numel(BAND). If it is
            %   scalar, all specified elements of 'DataCube' will be
            %   assigned with specified scalar value. ROW can be a colon
            %   character or vector of positive integers less than or equal
            %   to the number of rows in 'DataCube', colon means all rows
            %   in 'DataCube'. COLUMN can be a colon character or vector of
            %   positive integers less than or equal to the number of
            %   columns in the 'DataCube', colon means all columns in
            %   'DataCube'. BAND can be a colon character or vector of
            %   positive integers less than or equal to the number of
            %   spectral bands in the 'DataCube', colon means all bands in
            %   'DataCube'.
            %
            %   Example
            %   -------
            %   % Read pavia university hyperspectral data
            %   hcube = hypercube('paviaU');
            %
            %   % Assign normalized data
            %   data = rescale(hcube.DataCube);
            %   normalizedCube = assignData(hcube, ':',':',':',data);
            %
            %   See also hypercube
            
            
            % Validate hypercube object
            validateattributes(obj,{'hypercube'},{'scalar'},'obj','hypercube');
            if isempty(obj.DataCube)
                error(message('hyperspectral:hypercube:emptyDataCube'));
            end
            
            % Make a copy of input arguments to able to assign empty matrix
            tempRow = row;
            tempCol = col;
            tempBand = band;
            
            % ':' means all
            if strcmpi(row,':')
                row = 1:size(obj.DataCube,1);
            end
            if strcmpi(col,':')
                col = 1:size(obj.DataCube,2);
            end
            if strcmpi(band,':')
                band = 1:size(obj.DataCube,3);
            end
            
            dim = [size(obj.DataCube,1),size(obj.DataCube,2),size(obj.DataCube,3)];
            validateDimensions(row, col, band, dim);
            
            if ~isempty(newValues)
                if isscalar(newValues)
                    validateattributes(newValues, {'numeric'},...
                        {'nonsparse','real'},...
                        '','New values');
                else
                    validateattributes(newValues, {'numeric'},...
                        {'nonsparse','real','size',[numel(row),numel(col),numel(band)]},...
                        '','New values');
                end
            end
            
            if isempty(newValues)
                % Assign with empty matrix
                
                obj.DataCube(tempRow, tempCol, tempBand) = [];
                bandsToRemove = [];
                if (size(obj.DataCube,3) ~= obj.Metadata.Bands) || isempty(obj.DataCube)
                    % Update wavelength
                    obj.Wavelength(band) = [];
                    bandsToRemove = band;
                end
                
                if isempty(obj.DataCube)
                    dim = [0,0,0];
                else
                    dim = [size(obj.DataCube,1), size(obj.DataCube,2), size(obj.DataCube,3)];
                end
                % Update Metadata fields
                obj.Metadata = modifyMetadata(obj.Metadata, dim, bandsToRemove);
                
            else
                obj.DataCube(row, col, band) = newValues;
            end
        end
        
        function obj = cropData(obj, row, col, varargin)
            %cropData Crop DataCube
            %   NEWHCUBE = cropData(HCUBE, ROW, COLUMN) crops 'DataCube' in
            %   spatial dimension and returns a new hypercube object
            %   NEWHCUBE containing specified region.
            %
            %   NEWHCUBE = cropData(HCUBE, ROW, COLUMN, BAND) crops the
            %   'DataCube' in spatial and spectral dimension and returns a
            %   new hypercube object NEWHCUBE containing only specified
            %   region.
            %
            %   ROW can be a colon character or vector of positive integers
            %   less than or equal to the number of rows in 'DataCube',
            %   colon means all rows in 'DataCube'. COLUMN can be a colon
            %   character or vector of positive integers less than or equal
            %   to the number of columns in the 'DataCube', colon means all
            %   columns in 'DataCube'. BAND can be a colon character or
            %   vector of positive integers less than or equal to the
            %   number of spectral bands in the 'DataCube', colon means all
            %   bands in 'DataCube'.
            %
            %   Example
            %   -------
            %   % Read pavia university hyperspectral data
            %   hcube = hypercube('paviaU');
            %
            %   % Crop input volume to desired target size
            %   croppedData = cropData(hcube, [1:305],[1:170]);
            %
            %   See also hypercube
            
            narginchk(3,4);
            matlab.images.internal.errorIfgpuArray(obj, row, col, varargin{:});
            
            % Validate hypercube object
            validateattributes(obj,{'hypercube'},{'scalar'},'obj','hypercube');
            if isempty(obj.DataCube)
                error(message('hyperspectral:hypercube:emptyDataCube'));
            end
            
            % ':' means all
            if strcmpi(row,':')
                row = 1:size(obj.DataCube,1);
            end
            if strcmpi(col,':')
                col = 1:size(obj.DataCube,2);
            end
            if isempty(varargin) || (~isempty(varargin) && strcmpi(varargin{1},':'))
                band = 1:size(obj.DataCube,3);
            else
                band = varargin{1};
            end
            
            % Validate rows, columns, bands
            dim = [size(obj.DataCube,1),size(obj.DataCube,2),size(obj.DataCube,3)];
            validateDimensions(row, col, band, dim);
            
            obj.DataCube = obj.DataCube(unique(row, 'stable'), unique(col, 'stable'), unique(band, 'stable'));
            
            bandsToRemove = setdiff(1:obj.Metadata.Bands, band);
            
            % Update wavelength
            obj.Wavelength(bandsToRemove) = [];
            
            dim = [size(obj.DataCube,1), size(obj.DataCube,2), size(obj.DataCube,3)];
            % Update Metadata
            obj.Metadata = modifyMetadata(obj.Metadata,dim,bandsToRemove);
        end
        
        %% Remove bands
        function obj = removeBands(obj, method, band)
            %removeBands Remove specified spectral bands
            %    newhcube = removeBands(hcube,'Wavelength',range) removes
            %    hyperspectral data in the given wavelength range from
            %    hypercube object and returns a new hypercube newhcube with
            %    modified information. range is N-by-2 vector containing N
            %    number of wavelength ranges, each row represents minimum
            %    and maximum wavelength range values to be removed. These
            %    wavelength values must have same units as 'Wavelength'
            %    value in hypercube object.
            %
            %    newhcube = removeBands(hcube, 'BandNumber', band) removes
            %    the specified spectral band images and their corresponding
            %    metadata information from hypercube object and returns a
            %    new object with modified information. band is a positive
            %    integer vector with values less than or equal to number of
            %    bands in hypercube object hcube.
            %
            %   Example
            %   -------
            %   % Read Indian Pines hyperspectral data
            %   hcube = hypercube('indian_pines');
            %
            %   % Find band numbers of most informative bands
            %   sig = fippi(hcube, 5);
            %   [~, bandNum] = selectBands(hcube, sig,  'NumberOfBands', 200);
            %   % Get band numbers of least informative bands and remove
            %   % them from hypercube
            %   unUsableBands = setdiff(1:hcube.Metadata.Bands, bandNum);
            %   newHcube = removeBands(hcube, "BandNumber",unUsableBands);
            %
            %   See also hypercube, fippi.
            
            arguments
                obj (1,1) hypercube
                method {validatestring(method,{'BandNumber','Wavelength'})}
                band {mustBeNonempty, mustBeNonsparse, mustBeFinite, mustBeNumeric, mustBeReal, mustBePositive, mustBeOfSize(band, method)}
            end
            
            if isempty(obj.DataCube)
                error(message('hyperspectral:hypercube:emptyDataCube'));
            end
            
            method = validatestring(method,{'BandNumber','Wavelength'});
            
            if(strcmpi(method, 'wavelength'))
                if any(band(:,1) > band(:,2))
                    error(message('hyperspectral:hypercube:invalidRange'));
                end
                
                bandsToRemove = cell(size(band,1),1);
                for i = 1:size(band,1)
                    wlAboveMin = obj.Wavelength >= band(i,1);
                    wlBelowMax = obj.Wavelength <= band(i,2);
                    
                    bandsToRemove{i} = (find(wlAboveMin & wlBelowMax))';
                end
                % All bands which lies in given wavelength range
                band = [bandsToRemove{:}];
                band = unique(band);
                % Throw error if all provided wavelength regions are
                % non-overlapped with 'Wavelength' property.
                if(isempty(band))
                    error(message('hyperspectral:hypercube:invalidRegion',floor(obj.Wavelength(1)), ceil(obj.Wavelength(end))));
                end
                % Throw warning if there is at least one non-overlapped region
                if any(cellfun(@isempty, bandsToRemove))
                    warning(message('hyperspectral:hypercube:noOverlapRegion'));
                end
            end
            
            if any(band > size(obj.DataCube,3))
                error(message('hyperspectral:hypercube:mustBeLessThanOrEqual',size(obj.DataCube,3)))
            end
            
            obj.DataCube(:,:,band) = [];
            dim = [size(obj.DataCube,1), size(obj.DataCube,2), size(obj.DataCube,3)];
            if isempty(obj.DataCube)
                dim = [0,0,0];
            end
            obj.Wavelength(band) = [];
            obj.Metadata = modifyMetadata(obj.Metadata,dim,band);
        end
        
        %% Write in ENVI format
        function enviwrite(obj,filename, namedargs)
            %ENVIWRITE Write hypercube object as ENVI files
            %   enviwrite(HCUBE, FILENAME) writes hypercube object HCUBE to
            %   the file specified by FILENAME in ENVI file format.
            %   FILENAME is a character vector or string scalar that
            %   specifies the name of the file. ENVI file format consist of
            %   two files: binary image file, header file. The method
            %   enviwrite writes DataCube property to binary image file and
            %   Wavelength and Metadata properties to header file.
            %
            %   enviwrite(_____, PARAM, VALUE) writes hypercube object using
            %   specified parameter and value pairs.
            %
            %   Parameter - Value Pairs:
            %
            %    Interleave   - Data interleave, must be 'bil' or 'bip' or
            %                   'bsq'. (Default is 'bsq')
            %
            %    DataType     - Character vector or string defining the type
            %                   of data representation. See the help for
            %                   FWRITE for a list of valid datatypes.
            %                   Default is class of DataCube property.
            %
            %    ByteOrder    - Character vector or string defining byte
            %                   ordering for the file: 'ieee-le' or 'ieee-be'.
            %                   (Default is local machine format)
            %
            %    HeaderOffset - Numeric number, specifying number of bytes
            %                   to skip before the first data element. If
            %                   the file does not exist, ASCII null values
            %                   will be written to fill the space by default.
            %                   Default is 0.
            %
            %   Notes
            %   -----
            %   The method enviwrite writes two files with same filename
            %   but with different extensions having .dat for binary image
            %   file and .hdr for header file.
            %
            %   Example
            %   -------
            %   % Read Pavia University dataset
            %   hcube = hypercube('paviaU.hdr');
            %
            %   % Find first 20 most informative bands
            %   sig = fippi(hcube, 5);
            %   [newhcube, bandNum] = selectBands(hcube, sig, 'NumberOfBands', 20);
            %
            %   % Write informative bands to new ENVI file
            %   enviwrite(newhcube, 'newData');
            %
            %   See also hypercube, multibandwrite, fippi.
            
            
            arguments
                obj (1,1) hypercube
                filename (1,1) string {mustBeNonempty}
                namedargs.DataType (1,:) char {validatestring(namedargs.DataType,{'uint8','uint16','uint32','uint64','int16','int32','int64','single','double'})} = class(obj.DataCube)
                namedargs.ByteOrder (1,:) char {validatestring(namedargs.ByteOrder, {'ieee-le','ieee-be'})} = defaultByteOrder
                namedargs.HeaderOffset (1,1) double {mustBeNonsparse,mustBeReal, mustBeFinite,mustBeInteger,mustBeNonnegative} = 0
                namedargs.Interleave (1,:) char {validatestring(namedargs.Interleave, {'bsq','bip','bil'})} = 'bsq'
            end
            
            if isempty(obj.DataCube)
                error(message('hyperspectral:hypercube:emptyDataCube'));
            end
            
            interleave = validatestring(namedargs.Interleave,{'bsq','bip','bil'});
            dataType = validatestring(namedargs.DataType,{'uint8','uint16','uint32','uint64','int16','int32','int64','single','double'});
            byteOrder = validatestring(namedargs.ByteOrder, {'ieee-le','ieee-be'});
            offset = namedargs.HeaderOffset;
            
            % Base file name
            [~, ~, ext] = fileparts(filename);
            baseFileName = erase(filename, ext);
            
            % Write DataCube
            imgFile = strcat(baseFileName,".dat");
            multibandwrite(obj.DataCube,imgFile,interleave,'PRECISION',dataType,'OFFSET',offset,'MACHFMT',byteOrder);
            
            % Write header
            obj.Metadata.Height = size(obj.DataCube,1);
            obj.Metadata.Width = size(obj.DataCube,2);
            obj.Metadata.Bands = size(obj.DataCube,3);
            obj.Metadata.Interleave = interleave;
            % Convert data type and byte order to numeric
            obj.Metadata.DataType = setDataType(dataType);
            obj.Metadata.ByteOrder = setByteOrder(byteOrder);
            
            obj.Metadata.HeaderOffset = offset;
            obj.Metadata.Wavelength = obj.Wavelength;
            
            hdrFile = strcat(baseFileName,".hdr");
            enviheaderwrite(obj.Metadata, hdrFile);
        end
        
        %% Band select
        function [obj, bandNum] = selectBands(obj, signature, option)
            %selectBands Select most informative bands from DataCube.
            %
            % NEWHCUBE = SELECTBANDS(HCUBE, SIGNATURES) returns new object
            % NEWHCUBE with most informative bands by minimizing the linear
            % prediction error between original and estimated band. Input
            % HCUBE contains the hyperspectral data of size (M-by-N-by-C)
            % and SIGNATURES is the estimated endmembers with the size of
            % C-by-K matrix. Where K describes number of endmember
            % signatures present in the HCUBE. Output NEWHCUBE provides the
            % most informative cube of size (M-by-N-by-P) where P describes
            % reduced dimension.
            %
            % [NEWHCUBE, BANDNUM] = SELECTBANDS(____, 'PARAM', VAL) returns
            % NEWHCUBE with specific number of most informative bands.
            % BANDNUM provides band numbers of the selected bands. Parameter is:
            %
            %    'NumberOfBands' - Number of informative bands required. It
            %                      is a scalar value. Default value is
            %                      number of endmembers present in the
            %                      SIGNATURES.
            %
            % Notes
            % -----
            % 1) Band select algorithm takes 10 percentage of the original
            % 'DataCube' to avoid the computational complexity. Random
            % behaviour will not change the output bands in most of the
            % cases because each band is spatially highly correlated.
            %
            % 2) To select the most distinctive but informative bands,
            % water absorption and low SNR bands need to be preremoved
            % because they can be distinctive but not informative.
            %
            %   Class Support
            %   -------------
            %   The input array SIGNATURES must be one of the following
            %   classes: uint8, uint16, uint32, int8, int16, int32, single,
            %   double, uint64, int64. It must be real, finite and
            %   non-sparse. NUMOFBANDS is a positive integer-valued numeric
            %   scalar.
            %
            %   Example
            %   -------
            %   % Create hypercube object.
            %   hcube = hypercube('paviaU.dat');
            %
            %   % Estimate signatures using fippi algorithm.
            %   E = fippi(hcube, 9);
            %
            %   % Get the new hcube with prominent bands.
            %   newHcube = selectBands(hcube, E,'NumberOfBands',10);
            %
            %   References:
            %   -----------
            %   Q. Du, H. Yang, "Similarity-based unsupervised band
            %   selection for hyperspectral image analysis", IEEE Geoscience
            %   and Remote Sensing Letters, Vol- 5, Oct. 2008.
            %
            %   See also hypercube, fippi
            
            arguments
                obj {mustBeA(obj,'hypercube')}
                signature (:,:) {mustBeNumeric, mustBeNonNan, mustBeFinite, mustBeNonsparse, ...
                    mustBeNonempty, mustBeReal, mustBeLessSize(signature, obj)}
                
                option.NumberOfBands (1,1) {mustBeNumeric, mustBeNonNan, mustBeFinite, mustBeNonsparse, ...
                    mustBeNonempty, mustBeReal, mustBeInteger, mustBePositive, ...
                    mustBeLessSizeOfNumBands(option.NumberOfBands, obj)} = size(signature,2);
            end
            
            % Rescale the input arguments
            cube = obj.DataCube;
            numOfBands = option.NumberOfBands;
            if isinteger(cube)
                cube = single(cube);
                signature = single(signature);
            end
            
            % Get the dimension of the input data.
            [h, w, bands] = size(cube);
            
            if (h < 500 && w < 500)
                % Since hyperspectral data is very huge due to large number
                % of pixels in an image. Fortunately, band selection
                % process requires only small subset of pixels. Paper
                % considers only 10% of the data for the band selection. To
                % obtain the consistency random seeds are fixed.
                rng('default');
                randomPixelsForHeight = randi([1,h],1, round(h*0.1));
                rng('default');
                randomPixelsForWidth = randi([1,w],1,round(w*0.1));
            else
                % Since hyperspectral data has large number of pixels,
                % selecting only 20x20 region. Fortunately, band selection
                % process requires only small subset of pixels.
                rng('default');
                randomPixelsForHeight = randi([1,h],1, 20);
                rng('default');
                randomPixelsForWidth = randi([1,w],1, 20);
            end
            
            % Reduced cube based on randomly selected pixels.
            croppedCube = cube(randomPixelsForHeight', randomPixelsForWidth', :);
            [h, w, ~] = size(croppedCube);
            V = reshape(croppedCube, [h*w bands]);
            
            % Randomness behaviour might miss important signatures since
            % those signatures are really important for band selection
            % criteria. Check randomly selected pixels contain class
            % signature information and add it to reduced cube.
            for idx = 1:size(signature,2)
                B = ismember(V,signature(:,idx)','rows');
                unionCriteria = all(B, 'all');
                if ~unionCriteria
                    V = union(signature(:,idx)',V,'rows');
                end
            end
            
            % Size of the new cube.
            [hw, ~] = size(V);
            
            % Estimate first two bands.
            initialTwoBands = initialBandEstimation(signature);
            
            % OSP based method
            bandNum = initialTwoBands;
            E = V(:,initialTwoBands);
            
            for nn=1:numOfBands-2
                % Orthogonal Subspace Projection
                P = eye(hw) - (E * pinv(E));
                
                % Maximum absolute projection
                maxProj = sum((P*V).^2);
                
                [~, idx] = max(maxProj);
                bandNum = [bandNum; idx];  %#ok<AGROW>
                E = [E V(:,idx)]; %#ok<AGROW>
            end
            
            bandNum(numOfBands+1:end) = [];
            
            totalBands = 1:bands;
            
            % Remove other bands.
            bandToRemove = totalBands(~(ismember(totalBands,bandNum)));
            
            % If bandToRemove not empty
            if ~ isempty(bandToRemove)
                % Modify the obj.
                obj = removeBands(obj, 'BandNumber', bandToRemove);
            end
            
        end
        
        %% Colorize
        function [coloredImg, indices] = colorize(obj, bandNum, option)
            %COLORIZE Estimate colored image from DataCube.
            %
            % COLOREDIMG = COLORIZE(HCUBE) gives false colored image by
            % estimating most informative bands using selectBands method.
            % Input HCUBE contains the hyperspectral cube of size
            % (M-by-N-by-C). COLOREDIMG is false colored output of size
            % (M-by-N-by-3).
            %
            % COLOREDIMG = COLORIZE(____, BANDNUM) gives false colored
            % image based on given band numbers BANDNUM of the prominent
            % bands. Input HCUBE contains the hyperspectral cube of size
            % (M-by-N-by-C). COLOREDIMG is false colored output of size
            % (M-by-N-by-3). BANDNUM is 3-by-1 or 1-by-3 array.
            %
            % [COLOREDIMG, INDICES] = COLORIZE(___, 'PARAM', VAL) gives RGB
            % or CIR image based on the selected method. COLOREDIMG is
            % M-by-N-by-3 matrix, INDICES is 3-by-1 array, provides band
            % indices of the selected bands.
            %
            % Parameter includes:
            %
            % 'Method'   - Specifies the technique used to visualize the
            %              bands. Available options are:
            %
            %   'falsecolored'  Based on the selectBands method, finds
            %                   three most informative bands. Default is
            %                   'falsecolored'.
            %
            %   'rgb'           Divide spectral range into R, G and B
            %                   bands. Find representative bands based on
            %                   the correlation coefficient metric.
            %
            %   'cir'           Divide spectral range into NIR, R and G
            %                   bands. Find the representative bands based
            %                   on correlation coefficient metric.
            %
            % 'ContrastStretching' - Indicates contrast stretching. Choices
            %                        are:
            %
            %       false  - Returns original colored image. Default is
            %                false.
            %
            %       true   - Applies contrast limited adaptive histogram
            %                equalization.
            %
            % Notes
            % -----
            % 1) Sometimes the quality of the color display may be poor in
            % CIR or RGB bands because, those selected bands falling in the
            % specific spectral ranges. In such scenarios false colored
            % image is very much useful since that technique is independent
            % of the wavelength information. False colored technique tries
            % to find most three distinguished bands.
            %
            % 2) False colored technique is suitable if the objective is to
            % generate an image with unique colors. 
            %
            % 3) To create RGB or CIR images from hyperspectral cube, input
            % object must have representative RGB or CIR wavelengths in the
            % metadata.
            %
            % 4) INDICES return band numbers of the representative bands
            % based on user selected method. The elements contain in the
            % INDICES are always three.
            %
            % 5) Red band ranges from 600nm to 700nm, blue band ranges
            % from 400nm to 500nm, green band ranges from 500nm to 600nm
            % and NIR ranges from 760nm to 960nm.
            %
            % 6) Colorize method gives contrast stretched output when
            % ContrastStretching is true.
            %
            %   Class Support
            %   -------------
            %   The input BANDNUM must 1-by-3 or 3-by-1 vector with one of
            %   the following classes: uint8, uint16, uint32, int8, int16,
            %   int32, single, double, uint64, int64. It must be real,
            %   finite and non-sparse. METHOD is a character vector or
            %   string. CONTRASTSTRETCHING must be a logical value.
            %   COLOREDIMG output will be in single or double class.
            %   INDICES will be in double datatype.
            %
            %   Example - 1
            %   -----------
            %   % Create hypercube object.
            %   hcube = hypercube('paviaU.dat');
            %
            %   % Obtain false colored image.
            %   coloredImg = colorize(hcube);
            %
            %   % Display false colored image
            %   imshow(coloredImg)
            %
            %   Example - 2
            %   -----------
            %   % Create hypercube object.
            %   hcube = hypercube('paviaU.dat');
            %
            %   % Obtain contrast stretched RGB image.
            %   coloredImg = colorize(hcube, 'Method','rgb', 'ContrastStretching', true);
            %
            %   % Display RGB image
            %   imshow(coloredImg)
            %
            %   See also hypercube
            
            
            arguments
                obj {mustBeA(obj,'hypercube')}
                bandNum {mustBeNumeric, mustBeNonNan, mustBeFinite, mustBeNonsparse, ...
                    mustBeReal, mustBeInteger, mustBePositive, mustBeInRange(bandNum, obj)}= [];
                option.Method (1,:) char {validatestring(option.Method, {'rgb', 'cir', 'falseColored'})} = 'falseColored'
                option.ContrastStretching (1,1) {mustBeNumericOrLogical, mustBeReal, ...
                    mustBeNonsparse, mustBeNonNan, mustBeFinite, mustBeNonnegative} = false
            end
            
            cube = obj.DataCube;
            wavelength = obj.Wavelength;
            % Convert wavelength units to nanometer
            if ~isempty(char(obj.Metadata.WavelengthUnits))
                wavelength = hyper.internal.wavelengthUnitConversion(wavelength,obj.Metadata.WavelengthUnits);
            end      
            
            method = validatestring(option.Method, {'rgb', 'cir', 'falseColored'});
            contrastStretching = option.ContrastStretching;
            sz = size(cube);
            
            if sz(3) < 3
                error(message('hyperspectral:hypercube:insufficientBands', sz(3)));
            end
            
            if strcmpi(method, 'falseColored')
                wavelength = [];
                if sz(3) == 3
                    bandNum = [1;2;3];
                    [coloredImg, indices] = bandVisualization(cube, wavelength, bandNum, method, contrastStretching);
                else
                    if isempty(bandNum)
                        % Find the number of endmembers. Assuming there are
                        % 3 endmembers present in the data cube.
                        numSig = 3;
                        
                        % Obtain the signatures
                        E = fippi(cube, numSig);
                        
                        % Find the prominent bands using selectBands method
                        bestThreeBands = 3;
                        [~, bandNum] = selectBands(obj, E, 'numberOfBands', bestThreeBands);
                        
                        % Find the false colored image
                        [coloredImg, indices] = bandVisualization(cube, wavelength, bandNum, method, contrastStretching);
                    else
                        [coloredImg, indices] = bandVisualization(cube, wavelength, bandNum, method, contrastStretching);
                    end
                end
                
                
            else
                [coloredImg, indices] = bandVisualization(cube, wavelength, bandNum, method, contrastStretching);
            end
            if isrow(indices)
                indices = indices';
            end
        end
    end
end


function [imgfile, hdrfile, imgArray, wavelength, imgInfo] = parseInputs(img, hdrfileIn, wlIn)
% Initialize output parameters
imgfile = '';
hdrfile = '';
imgArray = [];
wavelength = [];
imgInfo = [];

% First argument can be numeric array or image file name
if isnumeric(img)
    % Numeric array
    imgArray = img;
else
    % File
    imgfile = matlab.images.internal.stringToChar(img);
    if ~isscalar(string(imgfile))
        error(message('hyperspectral:hypercube:mustBeScalarText','First'));
    end
    if ~exist(imgfile, 'file') && ~exist([imgfile,'.dat'], 'file') && ~exist([imgfile,'.img'], 'file')....
            && ~exist([imgfile,'.hdr'], 'file') && ~exist([imgfile,'.L1R'], 'file')
        error(message('hyperspectral:hypercube:invalidFile',imgfile));
    end
end

% Second argument can be Wavelength or header file name
if ~isempty(hdrfileIn)
    if isnumeric(hdrfileIn)
        wavelength = hdrfileIn;
        validateattributes(wavelength,{'numeric'},...
            {'nonsparse','finite','nonnan','real','positive','vector'},'hypercube','wavelength');
    else
        hdrfile = matlab.images.internal.stringToChar(hdrfileIn);
        if ~isscalar(string(hdrfile))
            error(message('hyperspectral:hypercube:mustBeScalarText','Second'));
        end
    end
end

% Third argument can be wavelength or metadata
if ~isempty(wlIn)
    if isempty(wavelength) && isnumeric(wlIn) && isvector(wlIn)
        wavelength = wlIn;
        validateattributes(wavelength,{'numeric'},...
            {'nonsparse','finite','nonnan','real','positive','vector'},'hypercube','wavelength');
    else
        imgInfo = wlIn;
        validateattributes(imgInfo, {'struct'},{},'hypercube','Metadata');
    end
end
end


function tf = isTIFF(file)
persistent istif;
if isempty(istif)
    % Cache the istiff private function
    fmts = imformats('tif');
    istif = fmts.isa;
end
tf = istif(file);
end


function [datatype, bytes] = setDataType(datatype)

switch datatype
    case 'uint8'
        datatype = 1;
        bytes = 1;
    case 'int16'
        datatype = 2;
        bytes = 2;
    case 'int32'
        datatype = 3;
        bytes = 4;
    case 'single'
        datatype = 4;
        bytes = 4;
    case 'double'
        datatype = 5;
        bytes = 8;
    case 'uint16'
        datatype = 12;
        bytes = 2;
    case 'uint32'
        datatype = 13;
        bytes = 4;
    case 'int64'
        datatype = 14;
        bytes = 8;
    case 'uint64'
        datatype = 15;
        bytes = 8;
    otherwise
        error(message('hyperspectral:hypercube:unknownDataType'));
end
end


function byteorder  = setByteOrder(byteOrder)

byteorder = 0;
if strcmpi(byteOrder, 'ieee-be')
    byteorder = 1;
end
end


function header = getMetaforImg(imgFile, fmt)

switch fmt
    case 'tif'
        tempInfo = imfinfo(imgFile);
        tempInfo = tempInfo(1,:);
    case 'ntf'
        tempInfo = nitfinfo(imgFile);
    case 'envi'
        header = enviinfo(imgFile);
        return;
    case 'mtl'
        header = hyper.internal.readmtlfile(imgFile);
        return;
end

header = hyper.internal.getReqFieldNames();

fid = fopen(imgFile, 'r');
header.Filename = string(fopen(fid));
fclose("all");

d = dir(header.Filename);

% File details
header.FileModDate   = string(d.date);
header.FileSize      = d.bytes;
header.Format = string(tempInfo.Format);
header.FormatVersion = tempInfo.FormatVersion;

if strcmpi(header.Format, "tif")
    header.Format = "tiff";
end

if strcmpi(fmt, 'tif')
    
    % DataCube details
    header.RasterFormat = header.Format;
    header.Height = tempInfo.Height;
    header.Width = tempInfo.Width;
    % In TIFF file, bands are stored band by band
    header.Interleave = "bsq";
    if strcmpi(tempInfo.ByteOrder, 'big-endian')
        header.ByteOrder = "ieee-be";
    else
        header.ByteOrder = "ieee-le";
    end
    
elseif strcmpi(fmt, 'ntf')
    
    % DataCube details
    header.RasterFormat = header.Format;
    if isfield(tempInfo, 'ImageSubheaderMetadata') && isfield(tempInfo.ImageSubheaderMetadata, 'ImageSubheader001')
        % Read first image subheader
        header.Height = tempInfo.ImageSubheaderMetadata.ImageSubheader001.NumberOfSignificantRowsInImage;
        header.Width = tempInfo.ImageSubheaderMetadata.ImageSubheader001.NumberOfSignificantColumnsInImage;
        % In NITF file, bands are stored band by band
        header.Interleave = "bsq";
        
        % Wavelength
        if isfield(tempInfo.ImageSubheaderMetadata.ImageSubheader001, 'BandMeta')
            
            bandMetaFields = fieldnames(tempInfo.ImageSubheaderMetadata.ImageSubheader001.BandMeta);
            
            wavelength = zeros(numel(bandMetaFields),1);
            for fieldNo = 1: numel(bandMetaFields)
                subFields = fieldnames(tempInfo.ImageSubheaderMetadata.ImageSubheader001.BandMeta.(bandMetaFields{fieldNo}));
                % BandSubcategory field in BandMeta structure contains wavelength
                hasWl = find(~cellfun(@isempty,(regexpi(subFields,'BandSubcategory\w*'))));
                
                if hasWl
                    tempwl = tempInfo.ImageSubheaderMetadata.ImageSubheader001.BandMeta.(bandMetaFields{fieldNo}).(subFields{hasWl});
                    if ~isempty(tempwl)
                        wavelength(fieldNo) = double(tempwl);
                    end
                end
            end
            if any(wavelength)
                header.Wavelength = wavelength;
            end
        end
    end
end
end


function imgfile = findActualImgFileName(origFile)
[~,~, imgExt] = fileparts(origFile);

if isequal(imgExt, '.hdr')
    origFile = erase(origFile, '.hdr');
end

if isequal(exist(origFile,'file'), 2)
    imgfile = origFile;
elseif isequal(exist([origFile, '.dat'], 'file'),2)
    imgfile = [origFile, '.dat'];
elseif isequal(exist([origFile, '.img'], 'file'),2)
    imgfile = [origFile, '.img'];
elseif isequal(exist([origFile, '.L1R'], 'file'),2)
    imgfile = [origFile, '.L1R'];
else
    error(message('hyperspectral:hypercube:invalidFile', origFile));
end
end


function validateMetadata(info)
% Required fields
fields = {'SensorType','Description','AcquisitionTime','RasterFormat'};
type = {'string','char'};
attrs = {'scalartext'};
fieldValidator(info, fields, type, attrs);

% Extra fields
fields = {'FWHM'};
type = {'single','double'};
attrs = {'nonsparse','finite','real','positive','numel',info.Bands,'vector'};
fieldValidator(info, fields, type, attrs);

fields = {'Gain','DataIgnoreValues','Offset','ReflectanceGain',...
    'ReflectanceOffset','SolarIrrandiance','WaterVapourAbsorption',...
    'OxygenAbsorption', 'OzoneAbsorption'};
type =  {'single','double'};
attrs = {'nonsparse','finite','real','numel',info.Bands,'vector'};
fieldValidator(info, fields, type, attrs);

% BandNames
fields = {'BandNames'};
type = {'string','char'};
attrs = {'numel',info.Bands};
fieldValidator(info, fields, type, attrs);

% Bad Bands
fields = {'BadBands'};
type = {'logical'};
attrs = {'finite','numel',info.Bands,'vector'};
fieldValidator(info, fields, type, attrs);

% Atmospheric parameters
fields = {'CloudCover','SunAzimuth','SunElevation','EarthSunDistance'};
type = {'numeric'};
attrs = {'nonsparse','finite','nonnan','real','scalar'};
fieldValidator(info, fields, type, attrs);
end


function fieldValidator(info, field, type, attr)
for No = 1:numel(field)
    if isfield(info, field{No}) && ~isempty(info.(field{No}))
        fieldName = getString(message('hyperspectral:hypercube:fieldName',field{No}));
        validateattributes(info.(field{No}), type, attr,'hypercube',fieldName);
    end
end
end


function obj = readL1RFile(obj, imgfile)
% Read Hyperion L1R dataset

dataInfo = hdfinfo(imgfile);

% Get interleave
interleave = [];
for structNo = 1: numel(dataInfo.Attributes)
    if strncmpi(dataInfo.Attributes(structNo).Name,'interleave',10)
        interleave = dataInfo.Attributes(structNo).Value;
        break;
    end
end

if isempty(interleave)
    error(message('hyperspectral:hypercube:noInterleave'));
end
interleave = validatestring(deblank(interleave),{'bil','bip','bsq'});
% Assign data order based on interleave
switch interleave
    case 'bsq'
        permOrder = [2,3,1];
    case 'bil'
        permOrder = [1,3,2];
    case 'bip'
        permOrder = [1,2,3];
end
    
obj.DataCube = permute(hdfread(dataInfo.SDS), permOrder);

% Update datacube information
[obj.Metadata.Height, obj.Metadata.Width, obj.Metadata.Bands] = size(obj.DataCube);
obj.Metadata.DataType = string(class(obj.DataCube));

% Find sensor type
if ~isempty(obj.Metadata.Description) && ...
        contains(obj.Metadata.Description, 'Hyperion', 'IgnoreCase', true)
    obj.Metadata.SensorType = "EO1 HYPERION";
    % Find band numbers
    bandNum = str2double(string(regexpi(obj.Metadata.BandNames,'(?<=Band *)[0-9]*','match')));
    validateattributes(bandNum, {'numeric'},...
        {'nonempty','nonsparse','finite','nonnan','real','vector','integer','>',0,'<=',242});
    % Assign solar irradiance, gain and bad band values for existing band numbers
    try         
        origBands = 1:242;
        [~, existedBandIdx] = ismember(bandNum,origBands);
        
        % Define bad band multipliers for Hyperion sensor
        % 1-7, 225-242 are unused bands. 58-76 fall in overlap region,
        % 77,78 contains high noise.
        BadBands = true(242,1);
        noisyBands = [1:7,225:242,58:76,77,78]';
        BadBands(noisyBands) = false;
        obj.Metadata.BadBands = BadBands(existedBandIdx);
        
        % VNIR scaling factor = 40, SWIR scaling factor = 80 for hyperion
        gain = 1./[repelem(40,70,1);repelem(80,172,1)];
        obj.Metadata.Gain = gain(bandNum);
        
        % Solar Irradiance
        solarIrradiance = hyper.internal.getSolarIrradiance(obj.Metadata.SensorType);
        obj.Metadata.SolarIrradiance = solarIrradiance(existedBandIdx); 
        
        % Simulated transmittance spectra of atmospheric Water vapor,
        % Ozone, Oxygen.
        [waterVapour, oxygen, ozone] = hyper.internal.getGaseousAbsorption(obj.Metadata.SensorType);
        obj.Metadata.WaterVapourAbsorption = waterVapour(existedBandIdx);
        obj.Metadata.OxygenAbsorption = oxygen(existedBandIdx);
        obj.Metadata.OzoneAbsorption = ozone(existedBandIdx);
    catch
        % Do nothing 
    end
end
end


function obj = enviread(obj, imgfile)
% Read envi image file

% Get required input arguments for multibandread
dim = [obj.Metadata.Height, obj.Metadata.Width, obj.Metadata.Bands];
% Add '*' to datatype to format the data in same class type (For example, *uint8)
datatype = strcat('*',obj.Metadata.DataType);
interleave = obj.Metadata.Interleave;
offset = obj.Metadata.HeaderOffset;
byteOrder = obj.Metadata.ByteOrder;

% Read hyperspectral datacube
obj.DataCube = multibandread(imgfile, dim, datatype, offset, interleave, byteOrder);
end


function validateImageFileFormat(imgfile)
%Get full path of the file
fid = fopen(imgfile);
imgfile_out = fopen(fid);
fclose('all');

if ~hdfh('ishdf',imgfile_out)
  error(message('MATLAB:imagesci:hdfinfo:invalidFile', imgfile));
end
end


function validateImageFile(imgfile, meta, l1rOffset)
% Image file size must be equal to the size computed using header
% information

fid = fopen(imgfile, 'r');
filename = fopen(fid);
fclose("all");

imgInfo = dir(filename);
imgFileSize = imgInfo.bytes;

% Size computed using header information
[~, bytes] = setDataType(meta.DataType);
calculatedSize = meta.Height*meta.Width*meta.Bands*bytes+meta.HeaderOffset+l1rOffset;

if ~isequal(imgFileSize,calculatedSize)
    error(message('hyperspectral:hypercube:incorrectFileSize',num2str(calculatedSize),num2str(imgFileSize)));
end
end


function header = modifyMetadata(header,dim,bandsToRemove)

% Raster details
header.Height = dim(1);
header.Width = dim(2);
header.Bands = dim(3);

if ~isempty(bandsToRemove)
    % Spectral details
    fields = {'BandNames', 'FWHM', 'Gain', 'DataIgnoreValues', 'Offset',...
        'ReflectanceGain', 'ReflectanceOffset', 'BadBands', 'TIFFile', 'SolarIrradiance',...
        'WaterVapourAbsorption', 'OxygenAbsorption', 'OzoneAbsorption'};
    for fieldNo = 1: numel(fields)
        if(isfield(header,fields{fieldNo}) && ~isempty(header.(fields{fieldNo})))
            header.(fields{fieldNo})(bandsToRemove) = [];
        end
    end
end
end


function enviheaderwrite(header, filename)

fid = fopen(filename,'w');

% Remove extra fields added by hypercube, Write only multispectral image
% data
extraFields = {'Filename','FileModDate','FileSize','Format','FormatVersion',...
    'RasterFormat', 'TIFFile', 'Panchromatic', 'ThermalBand', 'QualityBand'};
for fieldNo = 1: numel(extraFields)
    if isfield(header, extraFields{fieldNo})
        header = rmfield(header,extraFields{fieldNo});
    end
end

% Replace Metadata field names with ENVI supported field names
header = replaceFieldNames(header);

% Get map info stream
mapStream = [];
if isfield(header,'MapInfo')
    mapinfo = header.MapInfo;
    header = rmfield(header,'MapInfo');
    mapStream = getMapStream(mapinfo);
    if ~isempty(mapStream)
        % 32 - ASCII value for space, 61 - ASCII value for '='
        mapStream = [uint8('map info'),32,61,32, mapStream,10];
    end
end

% Get all fields stream (other than map info)
stream = getHeaderStream(header);

% [69,78,86,73] are ASCII values for 'ENVI'
% 10 is ASCII value for new line character
stream = [69,78,86,73,10,stream, mapStream];

fwrite(fid, stream);
fclose(fid);
end


function stream = getHeaderStream(header)

allFields = fieldnames(header);
% Initialize stream with cell arrays
tempStream = cell(numel(allFields),1);
for idx = 1:numel(allFields)
    thisFieldName = allFields{idx};
    thisField = header.(thisFieldName);
    
    if any(strcmpi(thisFieldName,{'FWHM','BBL'}))
        thisFieldName = lower(thisFieldName);
    end
    
    if ~isempty(thisField)
        % Convert Metadata structure field names to ENVI header fields.
        % (i.e. Convert field name to lower case and add space in front of
        % capital letter)
        thisFieldName = regexprep(thisFieldName, '([A-Z])', ' ${lower($1)}');
        thisFieldName = strtrim(thisFieldName);
        if iscell(thisField)
            % Cell array
            if isnumeric(thisField{1}) || islogical(thisField{1})
                % Numeric cell array
                
                % Extract numeric values from cells
                thisField = cell2mat(thisField);
                % Reshape into column array and convert numerical values to
                % string array
                thisField = reshape(thisField',[],1);
                thisField = num2str(thisField);
                % Convert string array to uint8 stream
                data = array2Stream(thisField);
                % (space - 32, equal - 61) [ex: lines = 100]
                tempStream{idx} = [uint8(thisFieldName),32,61,32,data,10];
            else
                % cell array containing character array or string
                if isstring(thisField{1})
                    thisField = matlab.images.internal.stringToChar(thisField);
                end
                % Convert char or string array to uint8 stream
                data = array2Stream(thisField);
                % (space - 32, equal - 61) [ex: interleave = bsq]
                tempStream{idx} = [uint8(thisFieldName),32,61,32,data,10];
            end
        elseif isnumeric(thisField) || islogical(thisField)
            % Numeric array
            % Reshape into column array and convert numerical values to
            % string array
            thisField = reshape(thisField',[],1);
            thisField = num2str(thisField);
            % Convert string array to uint8 stream
            data = array2Stream(thisField);
            tempStream{idx} = [uint8(thisFieldName),32,61,32,data,10];
        else
            % String array or character array
            thisField = matlab.images.internal.stringToChar(thisField);
            % Convert char or string array to uint8 stream
            data = array2Stream(thisField);
            % (space - 32, equal - 61) [ex: interleave = bsq]
            tempStream{idx} = [uint8(thisFieldName),32,61,32,data,10];
        end
    end
end
stream = [tempStream{:}];
end


function stream = array2Stream(array)
% Convert cell array containing string or character vector to uint8 stream
openBrace = 123;
closeBrace = 125;
if size(array,1) > 1
    % If the character array is having more than one row
    if ischar(array)
        array = cellstr(array);
    end
    newLine = 10;
    % Initialize
    tempStream = cell(size(array,1)-1,1);
    for idx = 1:size(array,1)-1
        % Add newline character after every 6 cells
        ind = ~logical(mod(idx,6));
        % [value, comma, space, newline]
        tempStream{idx} = [uint8(array{idx}),44,32, newLine(ind)];
    end
    stream = [openBrace, tempStream{:},32,uint8(array{idx+1}),closeBrace];
else
    % If the character array is having one row
    stream = uint8(array);
    if any(ismember(stream,44))
        % If there is comma(',') in stream, add curly braces[{,}] to stream
        stream = [openBrace, stream, closeBrace];
    end
end
end


function mapStream = getMapStream(mapinfo)
allFields = fieldnames(mapinfo);
data = cell(numel(allFields),1);
for fieldNo = 1:numel(allFields)
    thisField = mapinfo.(allFields{fieldNo});
    if ~isempty(thisField)
        if isnumeric(thisField)
            thisField = reshape(thisField',[],1);
            thisField = num2str(thisField);
            if numel(thisField) > 1
                subData = cell(size(thisField,1),1);
                for i = 1:size(thisField,1)
                    subData{i} = [uint8(thisField(i,:)),44,32];
                end
                data{fieldNo} = [subData{:}];
            else
                data{fieldNo} = [uint8(thisField),44,32]; % (comma - 44, space - 32)
            end
        else
            data{fieldNo} = [uint8(char(thisField)),44,32];
        end
    end
end
mapStream = [123, data{:}, 125]; % Put the map info in '{' and  '}'

if numel(mapStream)>2
    mapStream(end-2:end-1) = [];
else
    mapStream = [];
end
end


function header = replaceFieldNames(header)
fieldNames = {'Height','Width','Gain','Offset','ReflectanceGain','ReflectanceOffset','BadBands'};
supportedNames = {'Lines','Samples','DataGainValues','DataOffsetValues',...
    'DataReflectanceGainValues','DataReflectanceOffsetValues','BBL'};

totFields = numel(fieldNames);
for i = 1:totFields
    if isfield(header,fieldNames{i})
        header.(supportedNames{i}) = header.(fieldNames{i});
        header = rmfield(header,fieldNames{i});
    end
end
end


function initialTwoBands = initialBandEstimation(bandSignatures)
% Instead of exhaustively searching for the best initial or using a random
% initial, this algorithm initialized using two bands whose dissimilarity
% is largest based on maximum linear prediction error.
V = bandSignatures';

% paper suggest random selection, but sometimes it gets different results
% for the multiple call of the function. To avoid that we are finding the
% maximum length in the signature domain. It will help us in
% initialization.
initialSig = sum(V.^2);

[~, A1] = max(initialSig);
newBand = V(:,A1);

% Prepare a data without A1 band.
newdata = [V(:,1:A1-1) V(:,A1+1:end)];

% Create an orthogonal subspace w.r.t selected band.
createSubspace = eye(size(bandSignatures,2)) - newBand * pinv(newBand);

% Orthogonal projections to selected band.
orthogonalProj= createSubspace'*newdata;
l2norm = vecnorm(orthogonalProj);

% Find A2 band which is most dissimilar to A1.
[~, A2] = max(l2norm);

if A1 <= A2
    A2 = A2 + 1;
end

bandNum = [A1;A2];

% Projecting all other bands to orthogonal subspaces and find the A3.
% Continue this process till A1 and A3 are same. A1 = A3 means A1 and A2
% are the most significant dissimilar pair.
while (1)
    
    newBand = V(:,A2);
    newdata = [V(:,1:A2-1) V(:,A2+1:end)];
    createSubspace = eye(size(bandSignatures,2)) - newBand * pinv(newBand);
    orthogonalProj= createSubspace'*newdata;
    projections2 = vecnorm(orthogonalProj);
    [~, A3] = max(projections2);
    
    if A1 == A3
        break;
    end
    if A2 <= A3
        A1 = A2;
        A2 = A3 + 1;
    else
        A1 = A2;
        A2 = A3;
        
    end
    bandNum = [bandNum;A3]; %#ok<AGROW>
    
end

% Initial estimation of two bands
if length(bandNum) == 2
    initialTwoBands = bandNum;
else
    initialTwoBands = [bandNum(end-2);bandNum(end-1)];
end
end


function [I,bandNum] = bandVisualization(V, W, bandNum, method, contrastStretching)
% Channels covering the wavelengths from 400 nm to 700 nm for the RGB.
switch method
    case {'rgb', 'cir'}
        % Check if data is having RGB or CIR wavelengths.
        if isempty(bandNum)
            %  Divide R, G and B bands according to their wavelength range.
            rbands = hyper.internal.hasWavelength(W, {'red'});
            gbands = hyper.internal.hasWavelength(W, {'green'});
            
            % Find the closest wavelength from 650 nm for R band and 550 nm
            % for G band.
            [~, nearestRband] = min(abs(W(rbands)-650));
            [~, nearestGband] = min(abs(W(gbands)-550));
            
            % Get the closest band location of R and G channels.
            nearestRband = nearestRband + find(rbands,1) - 1;
            nearestGband = nearestGband + find(gbands,1) - 1;
            
            if strcmpi(method, 'rgb')
                % Divide B band according to their wavelength range.
                bbands = hyper.internal.hasWavelength(W, {'blue'});
                
                % Find the closest wavelength from 470 nm for B band.
                [~, nearestBband] = min(abs(W(bbands)-470));
                
                % Get the closest band location of B channels.
                nearestBband = nearestBband + find(bbands,1) - 1;
            else
                % Colored Infrared (CIR) color composite.
                % Divide spectral ranges into near infrared(NIR), R and G.
                % This method is more suitable to show vegetation in an
                % image.
                bbands = hyper.internal.hasWavelength(W, {'NIR'});
                
                % Find the closest wavelength from 860 nm for NIR band.
                [~, nearestBband] = min(abs(W(bbands)-860));
                
                % Get the closest band location of NIR channels.
                nearestBband = nearestBband + find(bbands,1) - 1;
            end
            
            R = nearestRband;
            G = nearestGband;
            B = nearestBband;
            
            % User selected RGB or CIR images.
        else
            
            R = hyper.internal.hasWavelength(W(bandNum), {'Red'});
            G = hyper.internal.hasWavelength(W(bandNum), {'Green'});
            
            if strcmpi(method, 'rgb')
                B = hyper.internal.hasWavelength(W(bandNum), {'Blue'});
            else
                B = hyper.internal.hasWavelength(W(bandNum), {'NIR'});
            end
            
            R = bandNum(R);
            G = bandNum(G);
            B = bandNum(B);
        end
        
        if strcmpi(method, 'rgb')
            % Check if B band is empty or not.
            if isempty(B)
                error(message('hyperspectral:hypercube:mustHaveBlueChannelWavelength'));
            end
        else
            % Check if NIR band is empty or not.
            if isempty(B)
                error(message('hyperspectral:hypercube:mustHaveNIRChannelWavelength'));
            end
        end
        
        % Check if G band is empty or not.
        if isempty(G)
            error(message('hyperspectral:hypercube:mustHaveGreenChannelWavelength'));
        end
        
        % Check if R band is empty or not.
        if isempty(R)
            error(message('hyperspectral:hypercube:mustHaveRedChannelWavelength'));
        end
        
        
        if strcmpi(method, 'rgb')
            bandNum = [R G B];
            [image1, image2, image3] = getImagePlanes(V, bandNum);
            
        else
            
            % For CIR method order should be NIR-R-G bands.
            bandNum = [B R G];
            [image1, image2, image3] = getImagePlanes(V, bandNum);
            
        end
        
        
    case 'falseColored'
        % Considering first three bands are most uncorrelated three bands in
        % hyperspectral data. So, extract first three bands from the
        % selectBands method.
        [image1, image2, image3] = getImagePlanes(V, bandNum);
        
end

I = stretching(image1, image2, image3, contrastStretching);
end


function [image1, image2, image3] = getImagePlanes(V, bandNum)
image1 = V(:,:,bandNum(1));
image2 = V(:,:,bandNum(2));
image3 = V(:,:,bandNum(3));
end


function I = stretching(image1, image2, image3, contrastStretching)
% Enhance the colored image using adaptive histogram equalization method.

% Concatenate three images
I = cat(3,image1,image2,image3);

if isinteger(I)
    I=single(I);
end

if contrastStretching
    for n=1:3
        stretchedImg = rescale(I(:,:,n));
        I(:,:,n)=adapthisteq(stretchedImg); % Contrast stretched image.
    end
    
else
    for n=1:3
        stretchedImg = rescale(I(:,:,n));
        I(:,:,n)=stretchedImg; % Original Image
    end
end
end


% Validators
function mustBeCharOrNum(img)
validateattributes(img, {'char','string','numeric'},{'nonsparse'},'hypercube');
end


function mustBeLessSize(arg2, arg1)

if ~ isequal(size(arg2,1),size(arg1.DataCube,3))
    error(message('hyperspectral:hypercube:invalidEndmembers', size(arg1.DataCube,3)));
end
end


function mustBeLessSizeOfNumBands(arg2, arg1)

if  ((0 > arg2) || (arg2 > size(arg1.DataCube,3)))
    error(message('hyperspectral:hypercube:mustBeLessThan', size(arg1.DataCube,3)));
end
end


function mustBeInRange(arg1, arg2)

if ~isempty(arg1)
    
    sz = size(arg2.DataCube);
    validateattributes(arg1, {'numeric'},...
        {'>', 0, '<=',sz(3), 'vector'},mfilename, 'band numbers',2);
    
    if ~(length(unique(arg1)) == 3)
        error(message('hyperspectral:hypercube:notEnoughSize'));
    end
end
end


function mustBeA(hypercube,classIn)

% Input must be a hypercube object.
validateattributes(hypercube,{classIn},{'scalar'},'obj','hypercube');
validateattributes(hypercube.DataCube, ...
    {'numeric'},{'nonempty','nonnan', 'finite', 'ndims', 3}...
    , mfilename, 'hypercube', 1);
end


function validateDimensions(row, col, ch, dim)
validateattributes(row, {'numeric'},...
    {'nonempty','nonsparse','finite','real','integer','vector','positive','<=',dim(1)},...
    '','First dimension');

validateattributes(col, {'numeric'},...
    {'nonempty','nonsparse','finite','real','integer','vector','positive','<=',dim(2)},...
    '','Second dimension');

validateattributes(ch, {'numeric'},...
    {'nonempty','nonsparse','finite','real','integer','vector','positive','<=',dim(3)},...
    '','Third dimension');
end


function mustBeOfSize(band, method)
if strncmpi(method, 'BandNumber', strlength(method))
    validateattributes(band, {'numeric'},{'integer','vector'},'','band numbers');
else
    validateattributes(band, {'numeric'},{'ncols',2},'','wavelength range');
end
end


function sensorType = validateSensorType(sensorType)
% SensorType consist of both satellite and sensor ID.
sensorType = strsplit(sensorType,' ');

if isempty(sensorType{1})
    error(message('hyperspectral:hypercube:mustBeNonempty','SPACECRAFT_ID'));
end

if ~(size(sensorType,2) > 1) || (size(sensorType,2) > 1 && isempty(sensorType{2}))
    error(message('hyperspectral:hypercube:mustBeNonempty','SENSOR_ID'));
end

satID = validatestring(sensorType{1},...
    ["EO1","LANDSAT_8","LANDSAT_7","LANDSAT_6","LANDSAT_5","LANDSAT_4","LANDSAT_3","LANDSAT_2","LANDSAT_1"],...
    mfilename, 'Spacecraft ID');

if strncmpi(satID,'Landsat',7)
    senID = validatestring(sensorType{2}, ["OLI_TIRS","ETM","TM","MSS"], mfilename, 'Sensor ID');
else
    senID = validatestring(sensorType{2}, ["HYPERION","ALI"], mfilename, 'Sensor ID');
end

sensorType = strcat(satID," ",senID);
end


function byteOrder = defaultByteOrder()

[~,~,endian] = computer;

if(endian == 'B')
    byteOrder = 'ieee-be';
else
    byteOrder = 'ieee-le';
end
end
