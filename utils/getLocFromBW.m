% getLocFromBW.m
% Author      : Tian Gao (tgaochn@gmail.com)
% Link        :
% Date        : Mon, 02/15/2021, 18:38
% Description :
%   get the important location from a BW image, including:
%       minX, maxX, minY, maxY, avgX, avgY
%
%   call function:
%   [minX, maxX, minY, maxY, avgX, avgY] = getLocFromBW(BW)
%%

%
function [minX, maxX, minY, maxY, avgX, avgY] = getLocFromBW(BW)
    %  ! ================================================ Header =====================================================
    %  ! ================================================ Main =======================================================
    [J, I] = find(BW);
    minX = min(I);
    maxX = max(I);
    minY = min(J);
    maxY = max(J);
    avgX = mean(I);
    avgY = mean(J);
end



