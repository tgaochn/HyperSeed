classdef log4ui < handle

    % addLog(app.logTextArea, '123')
    % addLog(app.logTextArea, '321')

    properties
        logger;
    end

    methods (Access = private)

        function self = log4ui()
        end
    end

    methods (Static)

        function obj = getLogger(ele)
            obj = log4ui();
            obj.logger = ele;
            obj.cleanLog();
        end
    end

    methods

        function info(self, logStr)
            formattedLogStr = sprintf('%s %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'), logStr);
            if isempty(self.logger.Value)
                self.logger.Value = formattedLogStr;
            else
                self.logger.Value = [self.logger.Value; formattedLogStr];
            end
        end

        function cleanLog(self)
            self.logger.Value = '';
        end
    end
end
