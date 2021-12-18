% Author      : Tian Gao (tgaochn@gmail.com)
% Link        :
% Date        : 2021/05/28, 01:23:54
% Description :
%
%%

%

function inputdatatypeDropDownOpening_(app, event)
    value = app.inputdatatypeDropDown.Value;
    if value == 1
        app.whitereferenceheaderpathEditField.Value = '';
        app.whitereferencedatapathEditField.Value = '';
        app.darkreferenceheaderpathEditField.Value = '';
        app.darkreferencedatapathEditField.Value = '';
    elseif value == 2
        app.whitereferenceheaderpathEditField.Value = 'whiteReference.hdr';
        app.whitereferencedatapathEditField.Value = 'whiteReference';
        app.darkreferenceheaderpathEditField.Value = 'darkReference.hdr';
        app.darkreferencedatapathEditField.Value = 'darkReference';
    end
end
