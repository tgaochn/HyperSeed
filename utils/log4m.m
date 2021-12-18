classdef log4m < handle
    % !! usage:
    % logger = log4m.getLogger();
    % logger = log4m.getLogger('log/run.log');
    % logger = log4m.getLoggerWithLV('log/run.log', lv);
    % logger = log4m.getLoggerWithLV('log/run.log', lv, disableFileOutput);

    % lv == 3: info warn error
    % lv == 4: warn error
    % lv == 5: error

    % logger.info('begin to preprocess images');

    % if no log file is given or para3 is false, only display log in command line.

    % ! output in run.log
    % c_preprocessImg 2020-09-13 02:08:41 INFO    begin to preprocess images

    % !! =====

    %LOG4M This is a simple logger based on the idea of the popular log4j.
    %
    % Description: Log4m is designed to be relatively fast and very easy to
    % use. It has been designed to work well in a matlab environment.
    % Please contact me (info below) with any questions or suggestions!
    %
    %
    % Author:
    %       Luke Winslow <lawinslow@gmail.com>
    % Heavily modified version of 'log4matlab' which can be found here:
    %       http://www.mathworks.com/matlabcentral/fileexchange/33532-log4matlab
    %

    properties (Constant)
        ALL = 0;
        TRACE = 1;
        DEBUG = 2;
        INFO = 3;
        WARN = 4;
        ERROR = 5;
        FATAL = 6;
        OFF = 7;
    end

    properties (Access = protected)
        logger;
        lFile;
        disableFileOutput;
    end

    properties (SetAccess = protected)
        fullpath = 'log4m.log';  %Default file
        commandWindowLevel = log4m.ALL;
        logLevel = log4m.INFO;
    end

    methods (Static)

        function obj = getLoggerOrig(logPath)
            %GETLOGGER Returns instance unique logger object.
            %   PARAMS:
            %       logPath - Relative or absolute path to desired logfile.
            %   OUTPUT:
            %       obj - Reference to signular logger object.
            %

            if (nargin == 0)
                logPath = 'log4m.log';
            elseif (nargin > 1)
                error('getLogger only accepts one parameter input');
            end

            persistent localObj;
            if isempty(localObj) || ~isvalid(localObj)
                localObj = log4m(logPath);
            end
            obj = localObj;
        end

        function obj = getLogger(logPath)
            if (nargin == 0)
                obj = log4m.getLoggerWithLV();
            else
                obj = log4m.getLoggerWithLV(logPath);
            end
        end

        function obj = getLoggerWithLV(logPath, logLevel, disableFileOutput)
            % function obj = getLoggerWithLV(logPath, logLevel)
            if (nargin == 0)
                logPath = 'run.log';
                logLevel = 0;
                disableFileOutput = true;
            elseif (nargin == 1)
                logLevel = 0;
                disableFileOutput = false;
            elseif (nargin == 2)
                disableFileOutput = false;
            end
            [fPath, ~, ~] = fileparts(logPath);
            if ~exist(fPath, 'dir')
                mkdir(fPath)
            end

            if ~exist(logPath, 'file') && ~disableFileOutput
                fclose(fopen(logPath, 'w'));
            end

            obj = log4m.getLoggerOrig(logPath);
            obj.setCommandWindowLevel(logLevel);
            obj.setLogLevel(logLevel);
            obj.disableFileOutput = disableFileOutput;
            obj.fullpath = logPath;
        end

        function name = getCallerInfo()
            [ST, ~] = dbstack();
            offset = min(size(ST, 1), 3);
            name = ST(offset).name;
        end

        function testSpeed(logPath)
            %TESTSPEED Gives a brief idea of the time required to log.
            %
            %   Description: One major concern with logging is the
            %   performance hit an application takes when heavy logging is
            %   introduced. This function does a quick speed test to give
            %   the user an idea of how various types of logging will
            %   perform on their system.
            %

            L = log4m.getLogger(logPath);

            disp('1e5 logs when logging only to command window');

            L.setCommandWindowLevel(L.TRACE);
            L.setLogLevel(L.OFF);
            tic;
            for i = 1: 1e5
                L.trace('log4mTest', 'test');
            end

            disp('1e5 logs when logging only to command window');
            toc;

            disp('1e6 logs when logging is off');

            L.setCommandWindowLevel(L.OFF);
            L.setLogLevel(L.OFF);
            tic;
            for i = 1: 1e6
                L.trace('log4mTest', 'test');
            end
            toc;

            disp('1e4 logs when logging to file');

            L.setCommandWindowLevel(L.OFF);
            L.setLogLevel(L.TRACE);
            tic;
            for i = 1: 1e4
                L.trace('log4mTest', 'test');
            end
            toc;

        end
    end

    %% Public Methods Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        function setFilename(self, logPath)
            %SETFILENAME Change the location of the text log file.
            %
            %   PARAMETERS:
            %       logPath - Name or full path of desired logfile
            %

            if ~self.disableFileOutput
                [fid, message] = fopen(logPath, 'a');

                if (fid < 0)
                    error(['Problem with supplied logfile path: ' message]);
                end
                fclose(fid);
            end

            self.fullpath = logPath;
        end

        function setCommandWindowLevel(self, loggerIdentifier)
            self.commandWindowLevel = loggerIdentifier;
        end

        function setLogLevel(self, logLevel)
            self.logLevel = logLevel;
        end

        %% The public Logging methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function trace(self, message)
            %TRACE Log a message with the TRACE level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            %
            name = self.getCallerInfo();
            self.writeLog(self.TRACE, name, message);
        end

        function debug(self, message)
            %TRACE Log a message with the DEBUG level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            %
            name = self.getCallerInfo();
            self.writeLog(self.DEBUG, name, message);
        end

        function info(self, message)
            %TRACE Log a message with the INFO level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            %
            name = self.getCallerInfo();
            self.writeLog(self.INFO, name, message);
        end

        function warn(self, message)
            %TRACE Log a message with the WARN level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            %
            name = self.getCallerInfo();
            self.writeLog(self.WARN, name, message);
        end

        function error(self, message)
            %TRACE Log a message with the ERROR level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            %
            name = self.getCallerInfo();
            self.writeLog(self.ERROR, name, message);
        end

        function fatal(self, message)
            %TRACE Log a message with the FATAL level
            %
            %   PARAMETERS:
            %       funcName - Name of the function or location from which
            %       message is coming.
            %       message - Text of message to log.
            %
            name = self.getCallerInfo();
            self.writeLog(self.FATAL, name, message);
        end

    end

    %% Private Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Unless you're modifying this, these should be of little concern to you.
    methods (Access = private)

        function self = log4m(fullpath_passed)

            if (nargin > 0)
                path = fullpath_passed;
            end
            self.setFilename(path);
        end

        %% WriteToFile
        function writeLog(self, level, callerNm, message)

            % set up our level string
            switch level
                    case{self.TRACE}
                    levelStr = 'TRACE';
                    case{self.DEBUG}
                    levelStr = 'DEBUG';
                    case{self.INFO}
                    levelStr = 'INFO';
                    case{self.WARN}
                    levelStr = 'WARN';
                    case{self.ERROR}
                    levelStr = 'ERROR';
                    case{self.FATAL}
                    levelStr = 'FATAL';
                otherwise
                    levelStr = 'UNKNOWN';
            end

            % !! If necessary, display log with all the level on screen
            outputStr = sprintf('%s %s %s %s\n'...
                                , callerNm ...
                            , datestr(now, 'yyyy-mm-dd HH:MM:SS') ...
                , levelStr ...
                , message);

            % use 3 colors showing msg type: info, warn, error
            switch level
                    case{self.WARN}
                    cprintf([1, 0.5, 0], '%s', outputStr) % orange for warn 
                    case{self.ERROR}
                    cprintf('err', '%s', outputStr);  % red for error
                otherwise
                    cprintf('text', '%s', outputStr);  % black for info
            end

            %If currently set log level is too high, just skip this log
            if (self.logLevel > level)
                return;
            end

            % Append new log to log file
            if ~self.disableFileOutput

                try
                    fid = fopen(self.fullpath, 'a');
                    fprintf(fid, ' %s %s %s %s\r\n'...
                            , callerNm ...
                        , datestr(now, 'yyyy-mm-dd HH:MM:SS') ...
                        , levelStr ...
                        , message);
                    fclose(fid);
                catch ME_1
                    disp(ME_1);
                end
            end
        end
    end

end
