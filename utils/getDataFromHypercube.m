% getDataFromHypercube.m
% Author      : Tian Gao (tgaochn@gmail.com)
% Link        :
% Date        : 2021/04/12, 17:44:29
% Description :
%
%%

function dataMatrix = getDataFromHypercube(hcubeObj, isRef, BWlist)
    %  ! ================================================ Header =====================================================
    p = inputParser;
    addRequired(p, 'hcubeObj');
    addRequired(p, 'isRef');
    if ~isRef
        addRequired(p, 'BWlist')
    end

    %  ! ================================================ Main =======================================================
    dataMatrix = [];  % store averaged data for table
    bandCnt = size(hcubeObj.Wavelength, 1);

    if ~isRef
        % ! for non-ref data, generate an averaged value for each band from each seed
        validSeedCnt = length(BWlist);
        for j = 1:validSeedCnt
            % get the filtered hypercube
            curSeedHcubeData = hcubeObj.DataCube;
            curSeedHcubeData(repmat(~BWlist{j}{1}, [1 1 bandCnt])) = 0;

            % get a seed-based averaged value
            pixelCnt = BWlist{j}{2};
            pixelwiseAvg = sum(curSeedHcubeData, [1, 2]) / pixelCnt;
            seedAvgValOnBands = reshape(pixelwiseAvg, [bandCnt, 1]);
            dataMatrix = [dataMatrix; seedAvgValOnBands'];
        end
    else
        % ! for ref data, generate an image-based averaged value for each band
        hcubeData = hcubeObj.DataCube;
        pixelCnt = size(hcubeData, 1) * size(hcubeData, 2);

        pixelwiseAvg = sum(hcubeData, [1, 2]) / pixelCnt;
        refAvgValOnBands = reshape(pixelwiseAvg, [bandCnt, 1]);
        dataMatrix = [dataMatrix; refAvgValOnBands'];
    end
end
