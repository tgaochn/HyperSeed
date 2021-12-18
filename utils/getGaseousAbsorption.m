function [waterVapour, oxygen, ozone] = getGaseousAbsorption(sensorType)
% Simulated transmittance spectra of atmospheric Water vapor, Ozone,
% Oxygen. The path is specified for a solar zenith angle of 50°, nadir
% viewing, a sea-level surface, and a sensor at the top of the atmosphere.
% Typical amounts of each gas are used. The spectra are simulated at a
% sampling interval of 1 nm and with aspectral resolution of 5 nm.

% Reference
% ---------
% https://www.sciencedirect.com/science/article/abs/pii/S0034425709000741

waterVapour = [];
oxygen = [];
ozone = [];

if contains(sensorType,'HYPERION','IgnoreCase',true)
    waterVapour = [1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 0.96 0.98 0.88 ,...
        1.00 1.00 1.00 0.98 1.00 0.90 0.95 1.00 1.00 0.84 0.85 0.93 0.61 0.71,...
        0.90 1.00 1.00 1.00 1.00 0.91 0.96 0.90 0.72 0.80 0.88 0.96 1.00 1.00,...
        0.98 0.70 0.76 0.60 0.77 0.10 0.40 0.28 0.47 0.76 0.80 0.96 0.94 0.93,...
        0.95 1.00 1.00 0.98 0.92 1.00 1.00 0.98 0.92 0.60 0.60 0.60 0.40 0.20,...
        0.20 0.40 0.60 0.80 0.92 0.97 0.97 0.98 1.00 1.00 1.00 0.98 0.95 0.87,...
        0.79 0.58 0.52 0.46 0.40 0.42 0.53 0.64 0.75 0.88 0.89 0.91 0.93 0.96,...
        0.98 1.00 1.00 1.00 1.00 1.00 1.00 0.96 0.91 0.86 0.80 0.40 0.02 0.02,...
        0.05 0.02 0.05 0.05 0.07 0.15 0.25 0.40 0.50 0.58 0.35 0.55 0.75 0.86,...
        0.90 0.95 0.98 0.99 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 0.99 0.97 0.94 0.88 0.88 0.85 0.82 0.65 0.50,...
        0.38 0.16 0.10 0.02 0.02 0.02 0.02 0.06 0.06 0.02 0.02 0.02 0.08 0.08,...
        0.10 0.18 0.29 0.42 0.50 0.70 0.90 0.92 0.96 0.97 0.98 0.96 0.97 0.97,...
        0.98 0.98 0.97 0.98 0.97 1.00 1.00 0.98 0.98 0.96 0.97 0.96 0.97 0.98,...
        0.99 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 0.99 0.98 0.98 0.94 0.90,...
        0.80 0.83 0.75 0.82 0.85 0.69 0.63 0.90 0.59 0.85 0.42 0.84 0.55 0.30,...
        0.18 0.15 0.12 0.09 0.06 0.03 0.03 0.03 ]';
    
    oxygen = [1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 0.84 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 0.22 0.60 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 0.97 0.93 0.90 0.80 0.90 0.93 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00]';
    
    ozone = [1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 0.99 0.98,...
        0.98 0.99 0.98 0.98 0.97 0.96 0.95 0.94 0.92 0.92 0.92 0.93 0.92 0.92,...
        0.93 0.95 0.96 0.97 0.98 0.98 0.98 0.99 0.99 0.99 0.99 0.99 0.99 0.99,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00,...
        1.00 1.00 1.00 1.00 1.00]';
end
end