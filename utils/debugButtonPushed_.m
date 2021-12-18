% Author      : Tian Gao (tgaochn@gmail.com)
% Link        :
% Date        : 2021/05/28, 01:23:54
% Description :
%
%%

%

function debugButtonPushed_(app, event)
    app.inputdatapathEditField.Value = 'E:\0_projData\5_hyperData\testForAPP\input\seed1\raw.hdr';
    app.outputdatapathEditField.Value = 'E:\0_projData\5_hyperData\testForAPP\output';
    app.inputdatatypeDropDown.Value = app.inputdatatypeDropDown.Items{1};
    app.whitereferencepathEditField.Value = '';
    app.darkreferencepathEditField.Value = '';
    app.whitereferencepathEditField.Enable = false;
    app.darkreferencepathEditField.Enable = false;
    app.selectfileButton_2.Enable = false;
    app.selectfileButton_3.Enable = false;
end
