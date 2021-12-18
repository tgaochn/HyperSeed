% Author      : Tian Gao (tgaochn@gmail.com)
% Link        :
% Date        : 2021/05/28, 01:23:54
% Description :
%
%%

%

function initLogger(app)
    if ~(isa(app.logger, 'log4m') || isa(app.logger, 'log4ui'))
        app.loggerCom = log4m.getLogger();
        app.loggerUI = log4ui.getLogger(app.logTextArea);
        % app.logger = app.loggerCom; % output log to command line for testing
        app.logger = app.loggerUI;
    end
end

