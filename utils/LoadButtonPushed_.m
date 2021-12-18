% Author      : Tian Gao (tgaochn@gmail.com)
% Link        :
% Date        : 2021/05/28, 01:23:54
% Description :
%
%%

%

function LoadButtonPushed_(app, event)
    clc;
    initLogger(app)
    logger = app.logger;

    % load app components
    inputHeaderPathField = app.imageheaderpathhdrEditField;
    inputImgPathField = app.imagedatapathEditField;
    slicedImgHistField = app.slicedImgHist;
    slicedImgField = app.slicedImgFig;
    sld = app.cubeSliceSlider;
    lamp = app.Lamp;
    runButton = app.RunButton;

    % variables
    inputHeaderFilenames = inputHeaderPathField.Value;
    inputImgFilenames = inputImgPathField.Value;
    initBandId = 1;

    % load data
    [~, headerLocalFn, headerExt] = fileparts(inputHeaderFilenames);
    [~, imgLocalFn, imgExt] = fileparts(inputImgFilenames);
    hcubeLocalFn = [imgLocalFn, imgExt];
    hcubeHeaderLocalFn = [headerLocalFn, headerExt];
    files = dir(inputHeaderFilenames);

    if length(files) == 0
        logger.info('No files found! Please check the input folder: ');
        logger.info(inputHeaderFilenames);
        return
    end

    hcubeFn = fullfile(files(1).folder, hcubeLocalFn);
    hcubeHeaderFn = fullfile(files(1).folder, hcubeHeaderLocalFn);

    checkFiles = {hcubeFn, hcubeHeaderFn};
    for curFn = checkFiles
        curFn = curFn{:};
        if ~exist(curFn, 'file')
            logger.info(['files do not exist: ', curFn]);
            return
        end
    end

    hcubeObj = hypercube(hcubeFn, hcubeHeaderFn);
    app.hcubeObj = hcubeObj;
    bandCnt = size(hcubeObj.Wavelength, 1);

    % init setting of slider
    sld.Visible = true;
    sld.Limits = [1 bandCnt];
    sld.Value = initBandId;

    % show init img in seedSegFigField
    slicingImg = hcubeObj.DataCube(:, :, initBandId);
    imagesc(slicedImgField, slicingImg);

    % show init histgram in slicedImgHist
    histogram(slicedImgHistField, slicingImg);

    % Lamp turns greed
    lamp.Color = 'green';

    % enable run button
    runButton.Enable = true;

    infoStr = sprintf('data loaded: %s', hcubeFn);
    logger.info(infoStr);
    curBandId = 1;
    curWavelength = app.hcubeObj.Wavelength(curBandId);
    infoStr = sprintf('current band ID: %d, wavelength: %.2f.', curBandId, curWavelength);
    logger.info(infoStr)
end
