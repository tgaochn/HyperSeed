function tf = isenvi(filename)
%ISENVI Check if a file is ENVI header.
%   TF = ISENVI(FILENAME) checks if the input FILENAME is a ENVI header
%   file or not. It returns true if it is, and false otherwise.
%
%   Example 1
%   ---------
%   % This example checks whether the file is ENVI header or not.
%   TF = isenvi('paviaU.hdr');
%
%   See also hypercube, enviinfo.

%   Copyright 2020 The MathWorks, Inc.


% Open the file.
[fid, msg] = fopen(filename, 'r');
if (fid < 0)
  error(msg);
end

header = fread(fid, 4, 'uint8=>uint8');
fclose(fid);

% Verify that the format is ENVI
% For ENVI file format, the first four characters must be 'ENVI'
tf = isequal(header, [69;78;86;73]);  % [69;78;86;73] are ASCII values for 'ENVI'
end
