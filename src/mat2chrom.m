function chrom = mat2chrom(mat, quant)
%MAT2CHROM encode a matrix into its corresponding chromosome
%
% chrom = MAT2CHROM(mat, quant)
% 
% INPUT:
%   mat:    matrix to be encoded
%   quant:  (optional) is the matrix quantized? [default = true]
% OUTPUT:
%   chrom:  1D chromosome
%
% See also CHROM2MAT GA_2D LOCAL_OPT

% Copyright 2015-2016 Antoine JUCKLER. All rights reserved.

if nargin < 2
    quant = 0;
end;

% Check that matrix is square
if size(mat, 1) ~= size(mat, 2)
    error 'mat is not a square matrix';
end;

% Check that matrix is symmetric
if ~isempty(mat(mat ~= mat(end:-1:1,:))) || ...
   ~isempty(mat(mat ~= mat(:,end:-1:1)))
    error 'mat is not symmetric';
else
    mat = mat(1:end/2, 1:end/2);
end;

% Check whether matrix is quantized
if sum(sum(mat(1:2:end,:) ~= mat(2:2:end,:))) == 0 ...
        && sum(sum(mat(:,1:2:end) ~= mat(:,2:2:end)))  == 0 ...
        && quant == 1
    chrom = reshape(mat(:,1:2:end), 1, numel(mat)/2);
    chrom = chrom(1:2:end);
elseif quant == 1
    error 'mat is not quantized';
else
    chrom = reshape(mat, 1, numel(mat));
end;

end