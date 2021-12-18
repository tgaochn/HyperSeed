% Author      : Tian Gao (tgaochn@gmail.com)
% Link        :
% Date        : 2021/05/28, 01:23:54
% Description :
%
%%

%

function inputdatatypeDropDownValueChanged_(app, event)
    totalDataTypeCnt = length(app.inputdatatypeDropDown.Items);
    value = app.inputdatatypeDropDown.Value;
    for i = 1:totalDataTypeCnt
        curItemVlaue = app.inputdatatypeDropDown.Items{i};
        if strcmp(value, curItemVlaue)
            value = i;
            break
        end
    end
    if value == 1
        app.whitereferenceheaderpathEditField.Value = '';
        app.whitereferencedatapathEditField.Value = '';
        app.darkreferenceheaderpathEditField.Value = '';
        app.darkreferencedatapathEditField.Value = '';
        app.whitereferenceheaderpathEditField.Enable = false;
        app.whitereferencedatapathEditField.Enable = false;
        app.darkreferenceheaderpathEditField.Enable = false;
        app.darkreferencedatapathEditField.Enable = false;
        app.selectfileButton_2.Enable = false;
        app.selectfileButton_3.Enable = false;
    elseif value == 2
        app.whitereferenceheaderpathEditField.Enable = true;
        app.whitereferencedatapathEditField.Enable = true;
        app.darkreferenceheaderpathEditField.Enable = true;
        app.darkreferencedatapathEditField.Enable = true;
        app.selectfileButton_2.Enable = true;
        app.selectfileButton_3.Enable = true;
    elseif value == 3
        if isempty(app.whitereferenceheaderpathEditField.Value)
            app.whitereferenceheaderpathEditField.Value = 'whiteReference.hdr';
            app.whitereferencedatapathEditField.Value = 'whiteReference';
        else
            [path, fn, ext] = fileparts(app.whitereferenceheaderpathEditField.Value);
            if ~isempty(path)
                app.whitereferenceheaderpathEditField.Value = [fn, ext];
                app.whitereferencedatapathEditField.Value = fn;
            end
        end

        if isempty(app.darkreferenceheaderpathEditField.Value)
            app.darkreferenceheaderpathEditField.Value = 'darkReference.hdr';
            app.darkreferencedatapathEditField.Value = 'darkReference';
        else
            [path, fn, ext] = fileparts(app.darkreferenceheaderpathEditField.Value);
            if ~isempty(path)
                app.darkreferenceheaderpathEditField.Value = [fn, ext];
                app.darkreferencedatapathEditField.Value = fn;
            end
        end

        app.whitereferenceheaderpathEditField.Enable = true;
        app.whitereferencedatapathEditField.Enable = true;
        app.darkreferenceheaderpathEditField.Enable = true;
        app.darkreferencedatapathEditField.Enable = true;
        app.selectfileButton_2.Enable = true;
        app.selectfileButton_3.Enable = true;
    end
end
