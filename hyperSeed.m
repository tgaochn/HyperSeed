classdef hyperSeed < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        Lamp                            matlab.ui.control.Lamp
        LoadButton                      matlab.ui.control.StateButton
        cubeSliceSlider                 matlab.ui.control.Slider
        UseTestingParaButton            matlab.ui.control.Button
        enableoverwritingCheckBox       matlab.ui.control.CheckBox
        selectfileButton_3              matlab.ui.control.Button
        darkreferenceheaderpathEditField  matlab.ui.control.EditField
        darkreferenceheaderpathEditFieldLabel  matlab.ui.control.Label
        selectfileButton_2              matlab.ui.control.Button
        whitereferenceheaderpathEditField  matlab.ui.control.EditField
        whitereferenceheaderpathEditFieldLabel  matlab.ui.control.Label
        removebandsCheckBox             matlab.ui.control.CheckBox
        selectpathButton                matlab.ui.control.Button
        outputdatapathEditField         matlab.ui.control.EditField
        outputpathLabel                 matlab.ui.control.Label
        selectfileButton                matlab.ui.control.Button
        RunButton                       matlab.ui.control.Button
        imageheaderpathhdrEditField     matlab.ui.control.EditField
        imageheaderpathhdrEditFieldLabel  matlab.ui.control.Label
        logTextArea                     matlab.ui.control.TextArea
        logTextAreaLabel                matlab.ui.control.Label
        inputdatatypeDropDown           matlab.ui.control.DropDown
        inputdatatypeDropDownLabel      matlab.ui.control.Label
        inputdatapathEditFieldLabel_2   matlab.ui.control.Label
        bandID                          matlab.ui.control.EditField
        minintensityLabel               matlab.ui.control.Label
        minInten                        matlab.ui.control.EditField
        maxintensityLabel               matlab.ui.control.Label
        maxInten                        matlab.ui.control.EditField
        histogramLabel                  matlab.ui.control.Label
        slicedimageLabel                matlab.ui.control.Label
        seedsegmentresultLabel          matlab.ui.control.Label
        minpixelforclusteringLabel      matlab.ui.control.Label
        minPixelForClustering           matlab.ui.control.EditField
        enableellipsefittingCheckBox    matlab.ui.control.CheckBox
        imagedatapathEditField          matlab.ui.control.EditField
        imagedatapathEditFieldLabel     matlab.ui.control.Label
        whitereferencedatapathEditField  matlab.ui.control.EditField
        whitereferencedatapathLabel     matlab.ui.control.Label
        darkreferencedatapathEditField  matlab.ui.control.EditField
        darkreferencedatapathEditFieldLabel  matlab.ui.control.Label
        segTestFig                      matlab.ui.control.UIAxes
        slicedImgFig                    matlab.ui.control.UIAxes
        slicedImgHist                   matlab.ui.control.UIAxes
    end

    properties (Access = private)
%         path1 = addpath(genpath('lib'));
%         path2 = addpath(genpath('lib/callback'));
        path3 = addpath(genpath('utils'));
    end

    properties (Access = public)
         hcubeObj;
         loggerCom;
         loggerUI;
         logger = nan;
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: RunButton
        function run(app, event)
            runButtonPushed_(app, event)
        end

        % Button pushed function: selectfileButton
        function selectFile(app, event)
            [file, path] = uigetfile('raw.hdr');
            if isequal(file,0)
                return
            end
            [~, baseName, ~] = fileparts(file); 
            app.imageheaderpathhdrEditField.Value = fullfile(path, file);            
            app.imagedatapathEditField.Value = fullfile(path, baseName);            
        end

        % Callback function
        function debugButtonPushed(app, event)
            debugButtonPushed_(app, event);
        end

        % Button pushed function: selectpathButton
        function selectPath(app, event)
            outputPath = uigetdir();
            if isequal(outputPath,0)
                return
            end
            app.outputdatapathEditField.Value = outputPath;
        end

        % Button pushed function: selectfileButton_2
        function selectfileButton_2Pushed(app, event)
            [file, path] = uigetfile('raw.hdr');
            if isequal(file,0)
                return
            end
            [~, baseName, ~] = fileparts(file); 
            app.whitereferenceheaderpathEditField.Value = fullfile(path, file);            
            app.whitereferencedatapathEditField.Value = fullfile(path, baseName);             
        end

        % Button pushed function: selectfileButton_3
        function selectfileButton_3Pushed(app, event)
            [file, path] = uigetfile('raw.hdr');
            if isequal(file,0)
                return
            end
            [~, baseName, ~] = fileparts(file); 
            app.darkreferenceheaderpathEditField.Value = fullfile(path, file);            
            app.darkreferencedatapathEditField.Value = fullfile(path, baseName);             
        end

        % Value changed function: inputdatatypeDropDown
        function inputdatatypeDropDownValueChanged(app, event)
            inputdatatypeDropDownValueChanged_(app, event);
        end

        % Callback function
        function inputdatatypeDropDownOpening(app, event)
            inputdatatypeDropDownOpening_(app, event);
        end

        % Button pushed function: UseTestingParaButton
        function UseTestingParaButtonPushed(app, event)
            UseTestingParaButtonPushed_(app, event);
        end

        % Value changed function: cubeSliceSlider
        function cubeSliceSliderValueChanged(app, event)
            cubeSliceSliderValueChanged_(app, event);
        end

        % Value changed function: LoadButton
        function LoadButtonPushed(app, event)
            LoadButtonPushed_(app, event);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1294 781];
            app.UIFigure.Name = 'MATLAB App';

            % Create Lamp
            app.Lamp = uilamp(app.UIFigure);
            app.Lamp.Position = [779 732 20 20];
            app.Lamp.Color = [1 0 0];

            % Create LoadButton
            app.LoadButton = uibutton(app.UIFigure, 'state');
            app.LoadButton.ValueChangedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Text = 'Load';
            app.LoadButton.FontSize = 14;
            app.LoadButton.FontWeight = 'bold';
            app.LoadButton.Position = [699 728 59 25];

            % Create cubeSliceSlider
            app.cubeSliceSlider = uislider(app.UIFigure);
            app.cubeSliceSlider.ValueChangedFcn = createCallbackFcn(app, @cubeSliceSliderValueChanged, true);
            app.cubeSliceSlider.Visible = 'off';
            app.cubeSliceSlider.Position = [940 539 150 3];

            % Create UseTestingParaButton
            app.UseTestingParaButton = uibutton(app.UIFigure, 'push');
            app.UseTestingParaButton.ButtonPushedFcn = createCallbackFcn(app, @UseTestingParaButtonPushed, true);
            app.UseTestingParaButton.FontSize = 14;
            app.UseTestingParaButton.FontWeight = 'bold';
            app.UseTestingParaButton.Position = [492 26 128 29];
            app.UseTestingParaButton.Text = 'Use Testing Para';

            % Create enableoverwritingCheckBox
            app.enableoverwritingCheckBox = uicheckbox(app.UIFigure);
            app.enableoverwritingCheckBox.Text = 'enable overwriting';
            app.enableoverwritingCheckBox.FontSize = 14;
            app.enableoverwritingCheckBox.FontWeight = 'bold';
            app.enableoverwritingCheckBox.Position = [5 420 145 22];
            app.enableoverwritingCheckBox.Value = true;

            % Create selectfileButton_3
            app.selectfileButton_3 = uibutton(app.UIFigure, 'push');
            app.selectfileButton_3.ButtonPushedFcn = createCallbackFcn(app, @selectfileButton_3Pushed, true);
            app.selectfileButton_3.FontSize = 14;
            app.selectfileButton_3.FontWeight = 'bold';
            app.selectfileButton_3.Enable = 'off';
            app.selectfileButton_3.Position = [692 568 89 25];
            app.selectfileButton_3.Text = 'select file';

            % Create darkreferenceheaderpathEditField
            app.darkreferenceheaderpathEditField = uieditfield(app.UIFigure, 'text');
            app.darkreferenceheaderpathEditField.HorizontalAlignment = 'center';
            app.darkreferenceheaderpathEditField.Enable = 'off';
            app.darkreferenceheaderpathEditField.Position = [209 567 461 30];

            % Create darkreferenceheaderpathEditFieldLabel
            app.darkreferenceheaderpathEditFieldLabel = uilabel(app.UIFigure);
            app.darkreferenceheaderpathEditFieldLabel.HorizontalAlignment = 'center';
            app.darkreferenceheaderpathEditFieldLabel.FontSize = 14;
            app.darkreferenceheaderpathEditFieldLabel.FontWeight = 'bold';
            app.darkreferenceheaderpathEditFieldLabel.Position = [5 572 185 22];
            app.darkreferenceheaderpathEditFieldLabel.Text = 'dark reference header path';

            % Create selectfileButton_2
            app.selectfileButton_2 = uibutton(app.UIFigure, 'push');
            app.selectfileButton_2.ButtonPushedFcn = createCallbackFcn(app, @selectfileButton_2Pushed, true);
            app.selectfileButton_2.FontSize = 14;
            app.selectfileButton_2.FontWeight = 'bold';
            app.selectfileButton_2.Enable = 'off';
            app.selectfileButton_2.Position = [692 644 89 25];
            app.selectfileButton_2.Text = 'select file';

            % Create whitereferenceheaderpathEditField
            app.whitereferenceheaderpathEditField = uieditfield(app.UIFigure, 'text');
            app.whitereferenceheaderpathEditField.HorizontalAlignment = 'center';
            app.whitereferenceheaderpathEditField.Enable = 'off';
            app.whitereferenceheaderpathEditField.Position = [209 640 461 30];

            % Create whitereferenceheaderpathEditFieldLabel
            app.whitereferenceheaderpathEditFieldLabel = uilabel(app.UIFigure);
            app.whitereferenceheaderpathEditFieldLabel.HorizontalAlignment = 'center';
            app.whitereferenceheaderpathEditFieldLabel.FontSize = 14;
            app.whitereferenceheaderpathEditFieldLabel.FontWeight = 'bold';
            app.whitereferenceheaderpathEditFieldLabel.Position = [5 645 193 22];
            app.whitereferenceheaderpathEditFieldLabel.Text = 'white reference header path';

            % Create removebandsCheckBox
            app.removebandsCheckBox = uicheckbox(app.UIFigure);
            app.removebandsCheckBox.Text = 'remove bands in beginning/end';
            app.removebandsCheckBox.FontSize = 14;
            app.removebandsCheckBox.FontWeight = 'bold';
            app.removebandsCheckBox.Position = [5 395 235 22];
            app.removebandsCheckBox.Value = true;

            % Create selectpathButton
            app.selectpathButton = uibutton(app.UIFigure, 'push');
            app.selectpathButton.ButtonPushedFcn = createCallbackFcn(app, @selectPath, true);
            app.selectpathButton.FontSize = 14;
            app.selectpathButton.FontWeight = 'bold';
            app.selectpathButton.Position = [692 493 89 25];
            app.selectpathButton.Text = 'select path';

            % Create outputdatapathEditField
            app.outputdatapathEditField = uieditfield(app.UIFigure, 'text');
            app.outputdatapathEditField.HorizontalAlignment = 'center';
            app.outputdatapathEditField.FontSize = 14;
            app.outputdatapathEditField.FontWeight = 'bold';
            app.outputdatapathEditField.Position = [209 491 461 30];

            % Create outputpathLabel
            app.outputpathLabel = uilabel(app.UIFigure);
            app.outputpathLabel.HorizontalAlignment = 'center';
            app.outputpathLabel.FontSize = 14;
            app.outputpathLabel.FontWeight = 'bold';
            app.outputpathLabel.Position = [5 496 86 22];
            app.outputpathLabel.Text = 'output path';

            % Create selectfileButton
            app.selectfileButton = uibutton(app.UIFigure, 'push');
            app.selectfileButton.ButtonPushedFcn = createCallbackFcn(app, @selectFile, true);
            app.selectfileButton.FontSize = 14;
            app.selectfileButton.FontWeight = 'bold';
            app.selectfileButton.Position = [602 728 80 25];
            app.selectfileButton.Text = 'select file';

            % Create RunButton
            app.RunButton = uibutton(app.UIFigure, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @run, true);
            app.RunButton.FontSize = 14;
            app.RunButton.FontWeight = 'bold';
            app.RunButton.Enable = 'off';
            app.RunButton.Position = [717 26 96 29];
            app.RunButton.Text = {'Run'; ''};

            % Create imageheaderpathhdrEditField
            app.imageheaderpathhdrEditField = uieditfield(app.UIFigure, 'text');
            app.imageheaderpathhdrEditField.HorizontalAlignment = 'center';
            app.imageheaderpathhdrEditField.FontSize = 14;
            app.imageheaderpathhdrEditField.FontWeight = 'bold';
            app.imageheaderpathhdrEditField.Position = [209 723 375 30];

            % Create imageheaderpathhdrEditFieldLabel
            app.imageheaderpathhdrEditFieldLabel = uilabel(app.UIFigure);
            app.imageheaderpathhdrEditFieldLabel.HorizontalAlignment = 'center';
            app.imageheaderpathhdrEditFieldLabel.FontSize = 14;
            app.imageheaderpathhdrEditFieldLabel.FontWeight = 'bold';
            app.imageheaderpathhdrEditFieldLabel.Position = [5 727 169 22];
            app.imageheaderpathhdrEditFieldLabel.Text = 'image header path (.hdr)';

            % Create logTextArea
            app.logTextArea = uitextarea(app.UIFigure);
            app.logTextArea.Position = [27 81 699 245];

            % Create logTextAreaLabel
            app.logTextAreaLabel = uilabel(app.UIFigure);
            app.logTextAreaLabel.HorizontalAlignment = 'right';
            app.logTextAreaLabel.FontSize = 14;
            app.logTextAreaLabel.FontWeight = 'bold';
            app.logTextAreaLabel.Position = [363 333 26 22];
            app.logTextAreaLabel.Text = 'log';

            % Create inputdatatypeDropDown
            app.inputdatatypeDropDown = uidropdown(app.UIFigure);
            app.inputdatatypeDropDown.Items = {'intensity w/o reference', 'reflectance w/t uniform reference'};
            app.inputdatatypeDropDown.ValueChangedFcn = createCallbackFcn(app, @inputdatatypeDropDownValueChanged, true);
            app.inputdatatypeDropDown.Position = [121 453 158 22];
            app.inputdatatypeDropDown.Value = 'intensity w/o reference';

            % Create inputdatatypeDropDownLabel
            app.inputdatatypeDropDownLabel = uilabel(app.UIFigure);
            app.inputdatatypeDropDownLabel.HorizontalAlignment = 'right';
            app.inputdatatypeDropDownLabel.FontSize = 14;
            app.inputdatatypeDropDownLabel.FontWeight = 'bold';
            app.inputdatatypeDropDownLabel.Position = [5 453 105 22];
            app.inputdatatypeDropDownLabel.Text = 'input data type';

            % Create inputdatapathEditFieldLabel_2
            app.inputdatapathEditFieldLabel_2 = uilabel(app.UIFigure);
            app.inputdatapathEditFieldLabel_2.HorizontalAlignment = 'center';
            app.inputdatapathEditFieldLabel_2.FontSize = 14;
            app.inputdatapathEditFieldLabel_2.FontWeight = 'bold';
            app.inputdatapathEditFieldLabel_2.Position = [295 456 55 22];
            app.inputdatapathEditFieldLabel_2.Text = 'band id';

            % Create bandID
            app.bandID = uieditfield(app.UIFigure, 'text');
            app.bandID.HorizontalAlignment = 'center';
            app.bandID.Position = [403 452 45 30];

            % Create minintensityLabel
            app.minintensityLabel = uilabel(app.UIFigure);
            app.minintensityLabel.HorizontalAlignment = 'center';
            app.minintensityLabel.FontSize = 14;
            app.minintensityLabel.FontWeight = 'bold';
            app.minintensityLabel.Position = [295 416 92 22];
            app.minintensityLabel.Text = 'min intensity';

            % Create minInten
            app.minInten = uieditfield(app.UIFigure, 'text');
            app.minInten.HorizontalAlignment = 'center';
            app.minInten.Position = [403 412 138 30];

            % Create maxintensityLabel
            app.maxintensityLabel = uilabel(app.UIFigure);
            app.maxintensityLabel.HorizontalAlignment = 'center';
            app.maxintensityLabel.FontSize = 14;
            app.maxintensityLabel.FontWeight = 'bold';
            app.maxintensityLabel.Position = [295 374 95 22];
            app.maxintensityLabel.Text = 'max intensity';

            % Create maxInten
            app.maxInten = uieditfield(app.UIFigure, 'text');
            app.maxInten.HorizontalAlignment = 'center';
            app.maxInten.Position = [403 370 138 30];

            % Create histogramLabel
            app.histogramLabel = uilabel(app.UIFigure);
            app.histogramLabel.FontSize = 14;
            app.histogramLabel.FontWeight = 'bold';
            app.histogramLabel.Position = [1101 731 73 22];
            app.histogramLabel.Text = 'histogram';

            % Create slicedimageLabel
            app.slicedimageLabel = uilabel(app.UIFigure);
            app.slicedimageLabel.FontSize = 14;
            app.slicedimageLabel.FontWeight = 'bold';
            app.slicedimageLabel.Position = [872 731 89 22];
            app.slicedimageLabel.Text = 'sliced image';

            % Create seedsegmentresultLabel
            app.seedsegmentresultLabel = uilabel(app.UIFigure);
            app.seedsegmentresultLabel.FontSize = 14;
            app.seedsegmentresultLabel.FontWeight = 'bold';
            app.seedsegmentresultLabel.Position = [978 484 141 22];
            app.seedsegmentresultLabel.Text = 'seed segment result';

            % Create minpixelforclusteringLabel
            app.minpixelforclusteringLabel = uilabel(app.UIFigure);
            app.minpixelforclusteringLabel.HorizontalAlignment = 'center';
            app.minpixelforclusteringLabel.FontSize = 14;
            app.minpixelforclusteringLabel.FontWeight = 'bold';
            app.minpixelforclusteringLabel.Position = [506 453 160 22];
            app.minpixelforclusteringLabel.Text = 'min pixel for clustering';

            % Create minPixelForClustering
            app.minPixelForClustering = uieditfield(app.UIFigure, 'text');
            app.minPixelForClustering.HorizontalAlignment = 'center';
            app.minPixelForClustering.Position = [696 449 79 30];

            % Create enableellipsefittingCheckBox
            app.enableellipsefittingCheckBox = uicheckbox(app.UIFigure);
            app.enableellipsefittingCheckBox.Text = 'enable ellipse fitting';
            app.enableellipsefittingCheckBox.FontSize = 14;
            app.enableellipsefittingCheckBox.FontWeight = 'bold';
            app.enableellipsefittingCheckBox.Position = [5 370 156 22];
            app.enableellipsefittingCheckBox.Value = true;

            % Create imagedatapathEditField
            app.imagedatapathEditField = uieditfield(app.UIFigure, 'text');
            app.imagedatapathEditField.HorizontalAlignment = 'center';
            app.imagedatapathEditField.FontSize = 14;
            app.imagedatapathEditField.FontWeight = 'bold';
            app.imagedatapathEditField.Position = [209 681 375 30];

            % Create imagedatapathEditFieldLabel
            app.imagedatapathEditFieldLabel = uilabel(app.UIFigure);
            app.imagedatapathEditFieldLabel.HorizontalAlignment = 'center';
            app.imagedatapathEditFieldLabel.FontSize = 14;
            app.imagedatapathEditFieldLabel.FontWeight = 'bold';
            app.imagedatapathEditFieldLabel.Position = [5 685 112 22];
            app.imagedatapathEditFieldLabel.Text = 'image data path';

            % Create whitereferencedatapathEditField
            app.whitereferencedatapathEditField = uieditfield(app.UIFigure, 'text');
            app.whitereferencedatapathEditField.HorizontalAlignment = 'center';
            app.whitereferencedatapathEditField.Enable = 'off';
            app.whitereferencedatapathEditField.Position = [209 603 461 30];

            % Create whitereferencedatapathLabel
            app.whitereferencedatapathLabel = uilabel(app.UIFigure);
            app.whitereferencedatapathLabel.HorizontalAlignment = 'center';
            app.whitereferencedatapathLabel.FontSize = 14;
            app.whitereferencedatapathLabel.FontWeight = 'bold';
            app.whitereferencedatapathLabel.Position = [5 608 191 22];
            app.whitereferencedatapathLabel.Text = 'white reference data path    ';

            % Create darkreferencedatapathEditField
            app.darkreferencedatapathEditField = uieditfield(app.UIFigure, 'text');
            app.darkreferencedatapathEditField.HorizontalAlignment = 'center';
            app.darkreferencedatapathEditField.Enable = 'off';
            app.darkreferencedatapathEditField.Position = [209 528 461 30];

            % Create darkreferencedatapathEditFieldLabel
            app.darkreferencedatapathEditFieldLabel = uilabel(app.UIFigure);
            app.darkreferencedatapathEditFieldLabel.HorizontalAlignment = 'center';
            app.darkreferencedatapathEditFieldLabel.FontSize = 14;
            app.darkreferencedatapathEditFieldLabel.FontWeight = 'bold';
            app.darkreferencedatapathEditFieldLabel.Position = [5 533 172 22];
            app.darkreferencedatapathEditFieldLabel.Text = 'dark reference data  path';

            % Create segTestFig
            app.segTestFig = uiaxes(app.UIFigure);
            app.segTestFig.XTick = [];
            app.segTestFig.XTickLabel = {'[ ]'};
            app.segTestFig.YTick = [];
            app.segTestFig.Position = [794 81 481 400];

            % Create slicedImgFig
            app.slicedImgFig = uiaxes(app.UIFigure);
            app.slicedImgFig.XTick = [];
            app.slicedImgFig.XTickLabel = {'[ ]'};
            app.slicedImgFig.YTick = [];
            app.slicedImgFig.Position = [812 551 192 174];

            % Create slicedImgHist
            app.slicedImgHist = uiaxes(app.UIFigure);
            app.slicedImgHist.XTick = [];
            app.slicedImgHist.XTickLabel = {'[ ]'};
            app.slicedImgHist.YTick = [];
            app.slicedImgHist.Position = [1041 551 192 174];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = hyperSeed

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end

