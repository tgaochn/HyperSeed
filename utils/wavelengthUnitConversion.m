function convWavelength = wavelengthUnitConversion(wavelength, currentUnit, targetUnit)
%WAVELENGTHUNITCONVERSION Convert wavelength units
%   CONVWAVELENGTH = WAVELENGTHUNITCONVERSION(__) converts wavelength unit
%   of WAVELENGTH from CURRENTUNIT to TARGETUNIT.
%   CURRENTUNIT and TARGETUNIT can be:
%       - "micrometers", "micrometer", "um"
%       - "nanometers", "nanometer", "nm"
%       - "millimeters", "millimeter", "mm"
%       - "centimeters", "centimeter", "cm"
%       - "meters", "meter", "m"
%   WAVELENGTHUNITCONVERSION returns the WAVELENGTH unchanged if
%   CURRENTUNIT is anything other than the specified list. WAVELENGTH can
%   be any numeric matrix. Default TARGETUNIT is "nanometer".
%
%   Example 1
%   ---------
%   % Read hyperspectral data
%   hCube = hypercube('jasperRidge2_R198.img');
%
%   % Convert wavelength unit to default unit (nanometer)
%   convWavelength = hyper.internal.wavelengthUnitConversion(hCube.Wavelength,...
%               hCube.Metadata.WavelengthUnits);
%
%   Example 2
%   ---------
%   % Convert array unit
%   convWavelength = hyper.internal.wavelengthUnitConversion([10 20 30 40 50],...
%       'm', 'cm');
%

%  Copyright 2020 The MathWorks, Inc.
arguments
    wavelength double {mustBeNumeric, mustBeNonNan, mustBeFinite,...
        mustBeNonsparse, mustBeNonempty, mustBeReal, mustBePositive}
    currentUnit (1, 1) string {mustBeNonempty}
    targetUnit (1, 1) string {validatestring(targetUnit, ["micrometers",...
        "micrometer", "um", "nanometers", "nanometer", "nm",...
        "millimeters", "millimeter", "mm", "centimeters", "centimeter",...
        "cm", "meters", "meter", "m"])} = "nanometer"
end

% Convert currentUnit and targetUnit to lower case
currentUnit = lower(currentUnit);
targetUnit = lower(targetUnit);

% Perform unit conversion
if strcmp(targetUnit, "micrometers")|| strcmp(targetUnit, "micrometer") || strcmp(targetUnit, "um")
    if strcmp(currentUnit, "nanometers") || strcmp(currentUnit, "nanometer") || strcmp(currentUnit, "nm")
        convWavelength = wavelength / 1000;
    elseif strcmp(currentUnit, "millimeters") || strcmp(currentUnit, "millimeter") || strcmp(currentUnit, "mm")
        convWavelength = wavelength * 1000;
    elseif strcmp(currentUnit, "centimeters") || strcmp(currentUnit, "centimeter") || strcmp(currentUnit, "cm")
        convWavelength = wavelength * 10000;
    elseif strcmp(currentUnit, "meters") || strcmp(currentUnit, "meter") || strcmp(currentUnit, "m")
        convWavelength = wavelength * 10^6;
    else
        convWavelength = wavelength;
    end
elseif strcmp(targetUnit, "nanometers")|| strcmp(targetUnit, "nanometer") || strcmp(targetUnit, "nm")
    if strcmp(currentUnit, "micrometers") || strcmp(currentUnit, "micrometer") || strcmp(currentUnit, "um")
        convWavelength = wavelength * 1000;
    elseif strcmp(currentUnit, "millimeters") || strcmp(currentUnit, "millimeter") || strcmp(currentUnit, "mm")
        convWavelength = wavelength * 10^6;
    elseif strcmp(currentUnit, "centimeters") || strcmp(currentUnit, "centimeter") || strcmp(currentUnit, "cm")
        convWavelength = wavelength * 10^7;
    elseif strcmp(currentUnit, "meters") || strcmp(currentUnit, "meter") || strcmp(currentUnit, "m")
        convWavelength = wavelength * 10^9;
    else
        convWavelength = wavelength;
    end
elseif strcmp(targetUnit, "millimeters")|| strcmp(targetUnit, "millimeter") || strcmp(targetUnit, "mm")
    if strcmp(currentUnit, "micrometers") || strcmp(currentUnit, "micrometer") || strcmp(currentUnit, "um")
        convWavelength = wavelength / 1000 ;
    elseif strcmp(currentUnit, "nanometers") || strcmp(currentUnit, "nanometer") || strcmp(currentUnit, "nm")
        convWavelength = wavelength / 10^6;
    elseif strcmp(currentUnit, "centimeters") || strcmp(currentUnit, "centimeter") || strcmp(currentUnit, "cm")
        convWavelength = wavelength * 10;
    elseif strcmp(currentUnit, "meters") || strcmp(currentUnit, "meter") || strcmp(currentUnit, "m")
        convWavelength = wavelength * 1000;
    else
        convWavelength = wavelength;
    end
elseif strcmp(targetUnit, "centimeters")|| strcmp(targetUnit, "centimeter") || strcmp(targetUnit, "cm")
    if strcmp(currentUnit, "micrometers") || strcmp(currentUnit, "micrometer") || strcmp(currentUnit, "um")
        convWavelength = wavelength / 10000 ;
    elseif strcmp(currentUnit, "nanometers") || strcmp(currentUnit, "nanometer") || strcmp(currentUnit, "nm")
        convWavelength = wavelength / 10^7;
    elseif strcmp(currentUnit, "millimeters") || strcmp(currentUnit, "millimeter") || strcmp(currentUnit, "mm")
        convWavelength = wavelength / 10;
    elseif strcmp(currentUnit, "meters") || strcmp(currentUnit, "meter") || strcmp(currentUnit, "m")
        convWavelength = wavelength * 100;
    else
        convWavelength = wavelength;
    end
else
    if strcmp(currentUnit, "micrometers") || strcmp(currentUnit, "micrometer") || strcmp(currentUnit, "um")
        convWavelength = wavelength / 10^6 ;
    elseif strcmp(currentUnit, "nanometers") || strcmp(currentUnit, "nanometer") || strcmp(currentUnit, "nm")
        convWavelength = wavelength / 1e+9;
    elseif strcmp(currentUnit, "millimeters") || strcmp(currentUnit, "millimeter") || strcmp(currentUnit, "mm")
        convWavelength = wavelength / 1000;
    elseif strcmp(currentUnit, "centimeters") || strcmp(currentUnit, "centimeter") || strcmp(currentUnit, "cm")
        convWavelength = wavelength / 100;
    else
        convWavelength = wavelength;
    end
end
end