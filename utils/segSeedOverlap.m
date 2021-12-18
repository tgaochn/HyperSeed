% Description :
%   given a BW image
%   output all the segmented overlapping seed considering shape. 
%   similar to clustering
%
%   call function:
%   [L, num] = segSeedOverlap(BW)
%%

function [L, num] = segSeedOverlap(BW)
    %  ! ================================================ Header =====================================================
    AICBIC_SELECTION = 1;  %Set AICBIC_SELECTION = 1, to use AIC is selected else BIC is used
    Overlaping = 0;
    %  ! ================================================ Main =======================================================
    [L, num] = bwlabel(BW);
    S = regionprops(L, 'Area');
    for i = 1: length(S)
        area(i) = S(i).Area;
    end
    MArea = median(area(area > 250));
    IClustTotal = zeros(size(BW, 1), size(BW, 2));
    for i = 1: num
        O = (L == i);

        s = regionprops(O, 'BoundingBox');
        apoX = round(s.BoundingBox(2));
        eosX = min(size(O, 1), apoX + round(s.BoundingBox(4)));
        apoY = round(s.BoundingBox(1));
        eosY = min(size(O, 2), apoY + round(s.BoundingBox(3)));
        Ocrop = O(apoX: eosX, apoY: eosY);

        if area(i) < 0.3 * MArea
            [~, EL, ~] = runMergeFitting1(Ocrop, AICBIC_SELECTION);
            IClust = Ocrop;
            totEll(i).EL = EL;
            totEll(i).NUMEllipses = 1;
        else
            [IClust, EL, NUMEllipses] = runMergeFitting(Ocrop, AICBIC_SELECTION, Overlaping);  %DEFA method
            % [IClust,EL,NUMEllipses] = runSlitFitting(Ocrop,AICBIC_SELECTION);
            totEll(i).EL = EL;
            totEll(i).NUMEllipses = NUMEllipses;
        end
        totEll(i).BoundBox = [apoX eosX apoY eosY];
        M = max(IClustTotal(:));
        Bit = IClust > 0;
        IClust = IClust + M * Bit;
        IClustTotal(apoX: eosX, apoY: eosY) = IClustTotal(apoX: eosX, apoY: eosY) + IClust;
    end

    L = IClustTotal;
    num = max(L(:));
end
