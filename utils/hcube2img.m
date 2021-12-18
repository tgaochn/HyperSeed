% hcube2img.m
% Author      : Tian Gao (tgaochn@gmail.com)
% Link        : 
% Date        : Tue, 01/19/2021, 01:09
% Description :
%   generate a false-color image of the hyperspectral data
% 
%   call function:
%   img = hcube2img(hcube)
%%

% 
function coloredImg = hcube2img(hcube)
    %  ! ================================================ Header =====================================================
    %  ! ================================================ Main =======================================================
    coloredImg = colorize(hcube);
    coloredImg = uint8(coloredImg * 256);
    %
end
