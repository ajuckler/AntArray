function chrom = mat2chrom(mat, quant)
% Transforms a matrix into a chromosome
%
% INPUT:
%   mat:    matrix to be encoded
%

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
if (sum(sum(mat(1:2:end,:) ~= mat(2:2:end,:))) || ...
    sum(sum(mat(:,1:2:end) ~= mat(:,2:2:end)))) && quant == 1
    chrom = reshape(mat, 1, numel(mat));
elseif quant == 1
    error 'mat is not quantized';
else
    chrom = reshape(mat(:,1:2:end), 1, numel(mat)/2);
    chrom = chrom(1:2:end);
end;

end