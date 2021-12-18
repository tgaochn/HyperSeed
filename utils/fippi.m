function endMembers = fippi(cube, num, option)
%FIPPI Extract endmember signatures using fast iterative pixel purity index.
%   ENDMEMBERS = FIPPI(CUBE, NUM) estimates the endmember signatures by
%   projecting a reduced data on initially estimated endmembers. NUM
%   describes the number of endmember signatures present in the CUBE. CUBE
%   is a numeric array with dimensions M-by-N-by-C containing hyperspectral
%   data or a hypercube object. M and N are rows and columns of the
%   hyperspectral cube respectively. C describes number of bands in the
%   hyperspectral cube. NUM is a non-zero positive integered numeric value.
%   Output ENDMEMBERS is a numeric matrix and it provides the endmember
%   signatures of dimensions C-by-K. Where K is final set of endmembers
%   estimated based on the stopping criteria.
%
%   ENDMEMBERS = FIPPI(___,Name, Value) estimates the endmember signatures
%   by projecting a reduced data on initially estimated endmembers. Data
%   reduction is done using two methods. Parameter is:
%
%   'ReductionMethod' - Method used to perform dimension reduction. Choices
%   are:
%
%       'PCA' Based on the maximum data variation Principal Component
%             Analysis (PCA) transforms large set of data into smaller one
%             by maintaining most of the information.
%
%       'MNF' Based on the Noise fractions Maximum Noise Fraction (MNF)
%             transforms large set of data into smaller one by maintaining
%             most of the information. MNF is a method which orders the
%             components according to image quality. The default method is
%             MNF.
%
%   Notes
%   -----
%   1) FIPPI is an iterative algorithm, where iterative rule is developed
%   to improve each of the iteration until it reaches a final set of
%   endmembers.
%   2) FIPPI does not provide NUM number of endmembers as output. It
%   provides K number of endmembers based on the FIPPI stopping criteria.
%   FIPPI has replacement rule to iteratively refine each iteration and is
%   terminated when a set of final skewers will not have any new
%   endmember.
%   3) NUM input argument is helpful to estimate the initial set of
%   endmembers using Automatic Target Generation Process (ATGP) algorithm.
%
%   Class Support
%   ------------- 
%   The input array CUBE must be one of the following classes: uint8,
%   uint16, uint32, int8, int16, int32, single, double, uint64, int64 or
%   hypercube object. It must be real and non-sparse.The output ENDMEMBERS
%   is an numeric matrix of size C-by-K.
%
%   References
%   ----------
%   CI Chang, A. Plaza, "A fast iterative algorithm for implementation of
%   pixel purity index", IEEE Geoscience and Remote Sensing Letters,
%   volume-3, pages: 63-67, 2006.
%
%   Example - 1 : Find endmember signatures.
%   -----------
%   % Load indianpines dataset
%   cube = load('indian_pines.mat');
%  
%   % Estimate number of endmembers.
%   num = countEndmembersHFC(cube.indian_pines);
%  
%   % Estimate and visualize the endmember signatures using FIPPI
%   endMembers = fippi(cube.indian_pines, num);
%   plot(endMembers);
%   xlabel('Band index')
%   ylabel('Signatures')
%   title('Endmembers')
%
%   Example - 2 : Find endmember signatures using 'PCA' dimension reduction
%   % method.
%   -----------
%   % Read hyperspectral cube
%   obj = hypercube('indian_pines.hdr');
%   
%   % Estimate number of endmembers.
%   num = countEndmembersHFC(obj);
%
%   % Estimate the endmember signatures using FIPPI method
%   endMembers = fippi(obj, num, 'ReductionMethod', 'pca');
%   plot(endMembers);
%   xlabel('Band index')
%   ylabel('Signatures')
%   title('Endmembers')
%
%   See also ppi, nfindr, hypercube, countEndmembersHFC.

%   Copyright 2020 The MathWorks, Inc.

arguments
    cube {mustBeA(cube,'hypercube')}
    num (1,1) {mustBeNumeric, mustBeNonNan, mustBeFinite, mustBeNonsparse, ...
        mustBeNonempty, mustBeReal, mustBeInteger, mustBePositive, ...
        mustBeLessSize(num, cube)}
    option.ReductionMethod (1,1) string ...
        {validatestring(option.ReductionMethod,{'MNF', 'PCA'})} = "MNF" 
end

% For compiler tests
if isdeployed
    rootDir = ctfroot;
else
    rootDir = matlabroot;
end
% Register resources/hyperspectral.
matlab.internal.msgcat.setAdditionalResourceLocation(rootDir);

% Input parsing
if isobject(cube)
    cube = cube.DataCube;
end

method = validatestring(option.ReductionMethod, {'MNF', 'PCA'});

% Get the size of a cube.
[rows, cols, channels] = size(cube);
numofPixels = rows*cols;

% Reshape the cube for extracting signatures.
originalCube = reshape(cube, numofPixels, channels);

if isinteger(cube)
    cube = single(cube);
    num = single(num);
end

if strcmpi(method, 'PCA')
    % Principal Component Analysis (PCA) on hyperspectral cube.
    cube = hyperpca(cube, num);
    
else 
    % Maximum noise fraction (MNF) transformation is used to determine the
    % intrinsic dimensionality of the image data, to remove noise in the
    % data and to reduce the computational requirement.
    cube = hypermnf(cube, num);
end

% Initialize skewers using ATGP algorithm
initializeSkewers = atgp(cube,numofPixels, num);

% FIPPI Algorithm
stoppingCriteria = false;
cube = reshape(cube, numofPixels, num);

% Stopping criteria to find the final endmembers.
while ~stoppingCriteria
    
    % Project reduced data on initially estimated endmembers.
    projPoints = abs(initializeSkewers*cube');
    
    % Find the extreme positions to form an extrema set.
    [~, maxLoc] = max(projPoints,[], 2);
    maxLoc = unique(maxLoc);
    
    % Check for the set of new endmembers.
    stoppingCriteria = all(ismember(cube(maxLoc,:), initializeSkewers), 'all');
    initializeSkewers = unique([cube(maxLoc, :);initializeSkewers],'sorted','rows');
end

endMembers = originalCube(maxLoc, :)';
end

function U = atgp(dimRed, numVariables, num)
% Automatic Target Generation Process (ATGP) algorithm estimate initial set
% of signatures based on orthogonal subspace projection.

% Reshape the reduced cube.
dimRed = reshape(dimRed, numVariables, num)';

% Find the first signature index.
[~, idx] = max(sum(dimRed.^2));

% Initialize and get the first endmember signature.
U = zeros(num, class(dimRed));
U(:,1) = dimRed(:,idx);

for n=1:num-1
    % Maximum absolute orthogonal subspace projection  
    [~, idx] = max(sum(((eye(num) - (U * pinv(U)))*dimRed).^2)); 
    U(:,n+1) = dimRed(:,idx);
end
end

function mustBeA(hypercube,classIn)

if isobject(hypercube)
    validateattributes(hypercube,{classIn},{'nonempty'},'fippi','hypercube');
else
    validateattributes(hypercube, ...
    {'uint8','int8','uint16','int16','uint32','int32','single','double', 'uint64', 'int64'},...
    {'nonsparse', 'nonempty', 'real','nonnan', 'finite', 'ndims',3}, mfilename, 'hypercube', 1);
end

end

% Validate the numComponents are less than number of bands.
function mustBeLessSize(a, b)
if isobject(b)
    cube = b.DataCube;
else
    cube = b;
end
sz = size(cube);
validateattributes(a, {'numeric'},...
    {'>', 0, '<=',sz(3)},mfilename, 'num',2);

end