% Author      : Tian Gao (tgaochn@gmail.com)
% Link        :
% Date        : 2021/05/28, 01:23:54
% Description :
%
%%

%

function UseTestingParaButtonPushed_(app, event)
    clc;

    % load app components
    inputHeaderPathField = app.imageheaderpathhdrEditField;
    inputImgPathField = app.imagedatapathEditField;
    whiteReferenceHeaderPathField = app.whitereferenceheaderpathEditField;
    whiteReferenceImgPathField = app.whitereferencedatapathEditField;
    darkReferenceHeaderPathField = app.darkreferenceheaderpathEditField;
    darkReferenceImgPathField = app.darkreferencedatapathEditField;
    outputDataPathField = app.outputdatapathEditField;
    inputDataTypeField = app.inputdatatypeDropDown;
    removeBandsCheckBoxField = app.removebandsCheckBox;
    enableOverwritingCheckBoxField = app.enableoverwritingCheckBox;
    bandIDField = app.bandID;
    minIntenField = app.minInten;
    maxIntenField = app.maxInten;
    selectfileButtonField_2 = app.selectfileButton_2;
    selectfileButtonField_3 = app.selectfileButton_3;
    minPixelInClusteringField = app.minPixelForClustering;
    enableellipsefittingCheckBoxField = app.enableellipsefittingCheckBox;

    % change values
    inputHeaderPathField.Value = 'data/input/seeds/seed1/raw.hdr';
    inputImgPathField.Value = 'data/input/seeds/seed1/raw';
    whiteReferenceHeaderPathField.Value = 'data/input/whiteRef/whiteReference.hdr';
    whiteReferenceImgPathField.Value = 'data/input/whiteRef/whiteReference';
    darkReferenceHeaderPathField.Value = 'data/input/darkRef/darkReference.hdr';
    darkReferenceImgPathField.Value = 'data/input/darkRef/darkReference';
    outputDataPathField.Value = 'data/output';
    inputDataTypeField.Value = 'reflectance w/t uniform reference';
    removeBandsCheckBoxField.Value = true;
    enableOverwritingCheckBoxField.Value = true;
    selectfileButtonField_2.Enable = true;
    selectfileButtonField_3.Enable = true;
    whiteReferenceHeaderPathField.Enable = true;
    whiteReferenceImgPathField.Enable = true;
    darkReferenceHeaderPathField.Enable = true;
    darkReferenceImgPathField.Enable = true;
    bandIDField.Value = int2str(20);
    minIntenField.Value = int2str(400);
    maxIntenField.Value = int2str(2000);
    minPixelInClusteringField.Value = int2str(500);
    enableellipsefittingCheckBoxField.Enable = true;

    initLogger(app)
end
