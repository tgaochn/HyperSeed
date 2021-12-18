function varargout = hasWavelength(wavelengths, bandNames)
%HASWAVELENGTH Check whether a particular wavelength band is present
%   TF = HASWAVELENGTH(WAVELENGTHS, BANDNAMES) returns a logical array with
%   elements set to logical 1 (true), where WAVELENGTHS value belongs to
%   the wavelength range of a particular band specified by BANDNAMES;
%   otherwise, the element is logical 0 (false). WAVELENGTHS is a positive
%   numeric vector. BANDNAMES can either be a character vector or string
%   scalar containing a single band name or can be a string or cell array
%   containing the different band names. BANDNAMES can be specified as one
%   or multiple of these:
%           - "RED"   - Visible red band (600-700 nm)
%           - "GREEN" - Visible green band (500-600 nm)
%           - "BLUE"  - Visible blue band (400-500 nm)
%           - "NIR"   - Near infrared band (760-960 nm)
%           - "SWIR1" - Short wave infrared 1 band (1550-1750 nm)
%           - "SWIR2" - Short wave infrared 2 band (2080-2350 nm)
%
%   Class Support
%   -------------
%   Supported classes for the input WAVELENGTHS are single, double, uint8,
%   uint16, uint32, uint64, int8, int16, int32, int64.
%
%   Notes
%   -----
%   [1] Wavelengths values must be in Nanometers.
%
%   Standard Band Definition(in nanometers)
%   --------------------------------------
%   Band      Min(LB)      Center      Max(UB)
%   ----      ------       ------      ------
%   Blue      400          470         500
%   Green     500          550         600
%   Red       600          650         700
%   NIR       760          860         960
%   SWIR1     1550         1650        1750
%   SWIR2     2080         2220        2350
%
%   Example 1
%   ---------
%   % Read hyperspectral data
%   hCube = hypercube('jasperRidge2_R198.img');
%
%   % Check if (R, G, B, NIR, SWIR1, SWIR2) bands are present
%   [hasRed, hasGreen, hasBlue, hasNIR, hasSWIR1, hasSWIR2] =...
%       hyper.internal.hasWavelength(hCube.Wavelength,...
%       {'red', 'green', 'blue', 'nir', 'swir1', 'swir2'});
%
%   % Visualize the result
%   figure
%   plot(hCube.Wavelength, hasRed, 'r', hCube.Wavelength, hasGreen, 'g',...
%       hCube.Wavelength, hasBlue, 'b', hCube.Wavelength, hasNIR,...
%       hCube.Wavelength, hasSWIR1, hCube.Wavelength, hasSWIR2,...
%       'LineWidth', 1.5)
%   legend({'Red', 'Green', 'Blue', 'NIR', 'SWIR1', 'SWIR2'},...
%       'Location', 'bestoutside')
%   yticks([0 1])
%   yticklabels({'Not present', 'Present'})
%   xlabel('Wavelength (Nanometers)')
%
%   Example 2
%   ---------
%   % Read hyperspectral data
%   hCube = hypercube('paviaU.dat');
%
%   % Check if SWIR-1 band is present
%   hasSWIR1 = hyper.internal.hasWavelength(hCube.Wavelength, 'swir1');

%  Copyright 2020 The MathWorks, Inc.

arguments
    wavelengths (:,1) {mustBeNumeric, mustBeNonNan, mustBeFinite,...
        mustBeNonsparse, mustBeNonempty, mustBeReal, mustBePositive}
    bandNames (:, 1) string {mustBeNonempty, mustBeASubset(bandNames)}
end

% Convert index names to lower case
bandNames = lower(bandNames);

% Calculate number of output variables
numOut = numel(bandNames);
varargout = cell(1, numOut);

for idx = 1:numOut
    switch bandNames(idx)
        case 'red'
            lb = 600;
            ub = 700;
        case 'green'
            lb = 500;
            ub = 600;
        case 'blue'
            lb = 400;
            ub = 500;
        case 'nir'
            lb = 760;
            ub = 960;
        case 'swir1'
            lb = 1550;
            ub = 1750;
        case 'swir2'
            lb = 2080;
            ub = 2350;
    end
    varargout{idx} = and(wavelengths >= lb, wavelengths <= ub);
end
end

function mustBeASubset(input)
for idx = 1:length(input)
    validatestring(input(idx), ["red", "green", "blue", "nir", "swir1",...
        "swir2"]);
end
end