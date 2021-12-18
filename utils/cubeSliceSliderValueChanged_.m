% Author      : Tian Gao (tgaochn@gmail.com)
% Link        :
% Date        : 2021/05/28, 01:23:54
% Description :
%
%%

%

function cubeSliceSliderValueChanged_(app, event)
    curBandId = floor(app.cubeSliceSlider.Value);
    slicingImg = app.hcubeObj.DataCube(:, :, curBandId);
    imagesc(app.slicedImgFig, slicingImg);
    histogram(app.slicedImgHist, slicingImg);
    curWavelength = app.hcubeObj.Wavelength(curBandId);

    % output the info of current band
    logger = app.logger;
    infoStr = sprintf('current band ID: %d, wavelength: %.2f.', curBandId, curWavelength);
    logger.info(infoStr)
end
