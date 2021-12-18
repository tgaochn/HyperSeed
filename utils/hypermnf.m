function [reducedCube, coeff] = hypermnf(cube, numComponents, method)
%HYPERMNF Maximum Noise Fraction (MNF) transformation on raw data.
%
%   REDUCEDCUBE = HYPERMNF(CUBE, NUMCOMPONENTS) returns first P channels of
%   noise adjusted principal components as REDUCEDCUBE, which is the
%   representation of M-by-N-by-C data matrix or hypercube object CUBE in
%   the noise adjusted principal component space. Rows and columns of
%   REDUCEDCUBE correspond to spatial dimensions, channels to components.
%   HYPERMNF, by default, centers the data. For the non-default options,
%   use the name/value pair arguments. The value P specifies NUMCOMPONENTS.
%
%   [REDUCEDCUBE, COEFF] = HYPERMNF(____) returns the noise adjusted
%   principal component transform coefficients for the M-by-N-by-C data
%   matrix or object CUBE of hypercube class. Each column of COEFF contains
%   coefficients for one noise adjusted principal component.
%
%   [_____] = HYPERMNF(___, Name, Value) specifies optional parameter
%   name/value pairs to control the computations. Parameters are:
%
%   'MeanCentered' - Indicator for centering the CUBE. Choices are:
%
%           true   - The default. HYPERMNF centers CUBE by subtracting off
%           channel mean.
%
%           false  - HYPERMNF without mean centered data.
%
%   Notes
%   -----
%   1) MNF transform is equivalent to a two-stage transformation in which
%   the data are first transformed so that the noise covariance matrix is
%   the identity matrix and the second stage is the principal components
%   transform.
%
%   2) MNF transformation gives better SNR than PCA transformation when
%   signal-dependent noise is present in the input data. PCA transformation
%   works bettter SNR than MNF transformation when gaussian white noise is
%   present in the input data.
%
%   Class Support
%   -------------
%   The input array CUBE must be one of the following classes: uint8,
%   uint16, uint32, uint64, int8, int16, int32, int64, single or double. It
%   must be real and non-sparse. The output COEFF is an array of size
%   C-by-P and REDUCEDCUBE is an array of size M-by-N-by-P.  The output
%   type is single unless the DataCube type is double, in which case the
%   output type is double.
%
%   Reference
%   ---------
%   A. Green, M. Berman, P. Switzer and D. Craig, "A transformation for
%   ordering multispectral data in terms of image quality with implications
%   for noise removal", IEEE Transactions on Geoscience and Remote Sensing,
%   vol. 26, pages 65-74.
%
%   Example
%   -------
%   % Load paviaU dataset
%   cube = load('paviaU.mat');
%
%   % Reduce spectral dimensions using mnf
%   reducedCube= hypermnf(cube.paviaU, 10);
%
%   % Display first 10 bands from original data
%   figure
%   montage(cube.paviaU, 'Indices', 1:10, 'DisplayRange',[]);
%   title('First 10 bands from original data');
%
%   % Display reduced cube
%   figure
%   montage(reducedCube,'DisplayRange',[])
%   title('Reduced Data');
%
%   See also hyperpca, inverseProjection.

% Copyright 2020 The MathWorks, Inc.

arguments
    cube {mustBeA(cube,'hypercube')}
    numComponents (1,1) {mustBeNumeric, mustBeNonNan, mustBeFinite, mustBeNonsparse, ...
        mustBeNonempty, mustBeReal, mustBeInteger, mustBeNonnegative, ...
        mustBeLessSize(numComponents, cube)}
    method.MeanCentered (1,1) {mustBeNumericOrLogical, mustBeReal, ...
        mustBeNonsparse, mustBeNonNan, mustBeFinite, mustBeNonnegative} = true
end

% Mean centered Name-value pair.
centered = method.MeanCentered;

% For compiler tests
if isdeployed
    rootDir = ctfroot;
else
    rootDir = matlabroot;
end
% Register resources/hyperspectral.
matlab.internal.msgcat.setAdditionalResourceLocation(rootDir);

% Get the original data type.
if isobject(cube)
    V = cube.DataCube;
else
    V = cube;
end

if isinteger(V)
    V = single(V);
    numComponents = single(numComponents);
end

% Get the dimension of the input data.
[h, w, channels] = size(V);

cube = reshape(V, [h*w channels]);

if centered
    u = sum(cube,1) ./ (h*w);
    cube = cube - u;
end

% Noise will be uncorrelated compared to signal for (1, 0) small spatial
% lag. (1, 0) or (0, 1) for salt-pepper noise, while (0, 1) is usually
% appropriate for horizontal striping noise.
dX = diff(V);

% Find the covariance of the noise.
V = computeCov(dX);

% Rotation for the noise covariance.
[dX, U, ~] = svd(V, 'econ');

% Consider only nonzero eigen values
getNonZeros = sum(~(diag(U)== 0));

U = U(1:getNonZeros,1:getNonZeros);

if isempty(U)
    % Noise whitening matrix when all the eigen values are zero.
    noiseWhitening = dX;
else
    % Noise whitening matrix
    noiseWhitening = dX;
    
    % Consider only non zero eigen values.
    nonzeroRotations = dX(1:getNonZeros,1:getNonZeros)/(sqrt(U));
    noiseWhitening(1:getNonZeros,1:getNonZeros) = nonzeroRotations;
end

% Get the noise adjusted data.
cube = reshape(cube, h*w, channels);
cube = (noiseWhitening'*cube');

% Second rotation by principal components.
[cube, coeff, ~] = hyperpca(reshape(cube',h, w, channels), numComponents, ...
    'MeanCentered', 0);

% Noise adjusted PCA transformation.
coeff = coeff'/(noiseWhitening);

% Cast to original data type.
reducedCube = cast(cube , class(coeff));
coeff = coeff';
end

function V = computeCov(cube)
[h, w, channels] = size(cube);
cube = reshape(cube, [h*w channels]);
V = cube;
V = (V'*V)/(h*w);
end

function mustBeA(hypercube,classIn)

if isobject(hypercube)
    validateattributes(hypercube,{classIn},{'nonempty'},'hypermnf','hypercube');
    datacube = hypercube.DataCube;
else
    datacube = hypercube;
end
validateattributes(datacube, ...
    {'numeric', classIn},...
    {'nonsparse', 'nonempty', 'real','nonnan', 'finite', 'ndims', 3}, mfilename, 'hypercube', 1);

sz = size(datacube);
if ~(all(sz>1))
    error(message('hyperspectral:hyperpca:requireVolumeData')); 
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
    {'>', 0, '<=',sz(3)},mfilename, 'numComponents',2);

end