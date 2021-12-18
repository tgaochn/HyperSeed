% label2rainbow.m
% Author      : Tian Gao (tgaochn@gmail.com)
% Link        : 
% Date        : Sat, 02/13/2021, 19:06
% Description :
%   turn a labeled BW into rainbow colored image
%
%   call function:
%   coloredLabelImg = label2rainbow(LabeledBW)
%%

%
function coloredLabelImg = label2rainbow(LabeledBW)
    %  ! ================================================ Header =====================================================
    %  ! ================================================ Main =======================================================
    coloredLabelImg = label2rgb(LabeledBW, @jet, [.5 .5 .5]);
end
