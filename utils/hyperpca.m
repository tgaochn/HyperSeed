function [reducedData, coeff, varianceRetained] = hyperpca(cube, numComponents, options)
%HYPERPCA Reduce spectral dimension of hyperspectral data using principal component analysis(PCA).
%   REDUCEDCUBE = HYPERPCA(CUBE, P) returns first P principal component
%   channels as REDUCEDCUBE, which is the representation of data cube
%   present in hypercube object or a 3-D (M-by-N-by-C) array CUBE in 
%   the principal component space.
%
%   [REDUCEDCUBE, COEFF] = HYPERPCA(___) returns the transform coefficients
%   for given input. Each column of COEFF contains coefficients for one
%   principal component and the columns are in descending order of principal
%   component variance.
%
%   [REDUCEDCUBE, COEFF, VARIANCERETAINED] = HYPERPCA(___) also returns the
%   percentage of the total variance retained by each principal component.
%
%   REDUCEDCUBE = HYPERPCA(___, Name, Value) specifies one or more
%   name-value parameters to set the options for principal component analysis.
%
%    'Method' - A character vector or string scalar specifying the method
%               to be used for principal component analysis. The value for
%               method can be any of these:
%        'svd'   - Singular Value Decomposition of the data.
%                  This is the default.
%        'eig'   - Eigenvalue decomposition of the covariance matrix. The
%                  eigen value decomposition method is comparatively faster
%                  than SVD when spatial dimensions of hyperspectral data
%                  is greater than the number of spectral bands, but less
%                  accurate because the condition number of the covariance
%                  is the square of the condition number of data.
%
%    'MeanCentered' - A logical value to set the option for mean-centering
%                     input data in order to compute principal components.
%                     Set the values as:
%         true   - HYPERPCA centers the data by subtracting off channel
%                  means before computing SVD or EIG. This is the default.
%         false  - HYPERPCA does not center the data.
%
%   Class Support
%   -------------
%   The input CUBE must be hypercube object or numeric array. It must be
%   real and non-sparse. P is a positive integer scalar specifying the
%   number of components desired in reducedCube. The output REDUCEDCUBE is
%   a 3D array of size M-by-N-by-P, COEFF is a matrix of size C-by-P and
%   VARIANCERETAINED is a vector of size P-by-1. The output type is single
%   unless the DataCube type is double, in which case the output type is
%   double.
%
%   Example
%   -------
%   % Read ENVI image and header file and create hypercube object
%   hcube = hypercube('paviaU.dat','paviaU.hdr');
%
%   % Reduce spectral dimensions using pca
%   reducedCube = hyperpca(hcube, 10);
%
%   % Display first 10 bands from original data
%   figure, montage(hcube.DataCube, 'Indices', 1:10, 'DisplayRange',[]);
%   title('First 10 bands from original data');
%
%   % Display reduced cube
%   figure, montage(reducedCube,'DisplayRange',[])
%   title('Reduced Data');
%
% See also inverseProjection, hypermnf.

% Copyright 2020 The MathWorks, Inc.

arguments
    cube {mustBeA(cube)}
    numComponents (1,1) {mustBeNumeric, mustBeNonNan, mustBeFinite, ...
        mustBeNonsparse, mustBeInteger, mustBePositive, mustBeReal, ...
        mustBeLessSize(numComponents, cube)}
    options.Method (1,1) string ...
        {validatestring(options.Method,{'svd','eig'})} = "svd"    
    options.MeanCentered (1,1) logical {mustBeNumericOrLogical, ...
        mustBeFinite, mustBeReal, mustBeNonsparse} = true
    
end
meanCentered = options.MeanCentered;
method = validatestring(options.Method,{'svd','eig'});

% For compiler tests
if isdeployed
    rootDir = ctfroot;
else
    rootDir = matlabroot;
end
% Register resources/hyperspectral.
matlab.internal.msgcat.setAdditionalResourceLocation(rootDir);

if isobject(cube)
    I = cube.DataCube;
else
    I = cube;
end

originalSize = size(I);
I = reshape(I,originalSize(1)*originalSize(2),originalSize(3));

[n, c] = size(I);
if isinteger(I)
    I = single(I);
end

infIdx = isinf(I);
nanIdx = isnan(I);
% Rows that contain NaN
rowsWithNaN = any(nanIdx,2);
DOF = max(0,n-meanCentered-sum(rowsWithNaN));

if all(nanIdx(:))
    coeff = NaN;
    reducedData = NaN;
    varianceRetained = NaN;
    return;
end

if all(infIdx(:))
    coeff = [];
    reducedData = [];
    varianceRetained = [];
    return;
end

if meanCentered
    %  sample mean
    if any(rowsWithNaN)
        mu = mean(I,'omitnan');
    else
        mu = mean(I);
    end
    I = I-mu;
end

switch method
    case 'svd'
        % Remove all rows of I that contain NaN values.
        if any(rowsWithNaN)
            I = removeNaNs(rowsWithNaN, I);
        end
        [U,sigma,coeff] = svd(I,'econ');
        sigma = diag(sigma);
        reducedData =  U.*sigma';
        % Insert NaNs back
        if any(rowsWithNaN)
            reducedData = insertNaNs(rowsWithNaN, reducedData);
        end
        
        if nargout > 2
            variance = sigma.^2./DOF;
            varianceRetained = 100*variance/sum(variance);
            varianceRetained(numComponents+1:end,:) = [];
        end
        
    case 'eig'
        ICov = nonCentCov(I, rowsWithNaN, meanCentered);
        [coeff, eigValueDiag] = eig(ICov);
        [eigValues, idx] = sort(diag(eigValueDiag), 'descend');
        coeff = coeff(:, idx);
        if (DOF<c)
            coeff(:, DOF+1:c) = [];
            eigValues(DOF+1:c, :) = [];
        end
        % Set small negative eigen values to zero
        eigValues((eigValues<0)&(abs(eigValues)<(eps(eigValues(1))*length(eigValues)))) = 0;  
        reducedData = I/coeff';
        if nargout > 2
            variance = eigValues;
            varianceRetained = 100*variance/sum(variance);
            varianceRetained(numComponents+1:end,:) = [];
        end
end

% Return first numComponents
if (numComponents<=DOF)
    coeff(:, numComponents+1:end) = [];
    reducedData(:, numComponents+1:end) = [];
end

% Enforce a sign convention on the coefficients -- the largest element in
% each column will have a positive sign.
[~,maxind] = max(abs(coeff), [], 1);
[d1, d2] = size(coeff);
colsign = sign(coeff(maxind + (0:d1:(d2-1)*d1)));
coeff = coeff.*colsign;
reducedData = reducedData.*colsign;
reducedData = reshape(reducedData,originalSize(1),originalSize(2),numComponents);
end


function c = nonCentCov(x,rowsWithNaN,centered)
% Non centered covariance
%
%   C = nonCentCov(X,false) returns X'*X/N, where N is the number of
%   observations after removing rows missing values.
%
%   C = nonCentCov(X,true) returns X'*X/N-1. If data X is centered.

[n, ~] = size(x);
if ~any(rowsWithNaN)
    c = x'*x/(n-centered);
else
    xNotNaN = x((~rowsWithNaN),:);
    denom = max(0, (size(xNotNaN,1)-centered) );
    c = xNotNaN'*xNotNaN/denom;
end
end


function x = removeNaNs(rowsWithNaN, x)
% Remove all rows of X that contain NaN values.
t = ~rowsWithNaN;
x= x(t,:);
end


function X = insertNaNs(rowsWithNaN, x)
% Insert missing values in x and return it as X.
t = ~rowsWithNaN;
len = length(rowsWithNaN);
[~,c] = size(x);
X = nan([len,c],class(x));
X(t,:) = x;
end


function mustBeA(hypercube)
if isobject(hypercube)
    validateattributes(hypercube,{'hypercube'},{'scalar'},'hyperpca',...
        'hypercube');
    cube = hypercube.DataCube;
else
    cube = hypercube;
end
validateattributes(cube, ...
    {'numeric','hypercube'}, {'nonsparse', 'nonempty', 'real', ...
    'ndims',3}, mfilename, 'hypercube', 1);

sizeCube = size(cube);
if ~(all(sizeCube > 1))
    error(message('hyperspectral:hyperpca:requireVolumeData'));
end
end


function mustBeLessSize(numComponents, data)
if isobject(data)
    cube = data.DataCube;
else
    cube = data;
end
if numComponents > size(cube,3)
    error(message('hyperspectral:hyperpca:notLessEqual'));
end
end
