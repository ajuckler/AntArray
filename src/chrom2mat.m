function C = chrom2mat(chrom, quant)
%CHROM2MAT transforms a chromosome to its matrix representation
%
% C = CHROM2MAT(chrom, quant)
% 
% INPUT:
%   chrom:  chromosome to be reshaped
%   quant:  (optional) is the end matrix quantized? [default = true]
% OUTPUT:
%   C:      2D matrix representation of the chromosome
%
% See also MAT2CHROM GA_2D LOCAL_OPT

% Copyright 2015-2016 Antoine JUCKLER. All rights reserved.

if nargin < 2
    quant = 1;
end;

% Check chromosome length
chrom_sz = sqrt(length(chrom));
if abs(round(chrom_sz) - chrom_sz) > 10^-6
    error 'Invalid chromosome length';
else
    chrom_sz = round(chrom_sz);
end;

if ~quant
    C = zeros(chrom_sz*2);
    C(1:end/2, 1:end/2) = reshape(chrom, chrom_sz, chrom_sz);
else
    C = zeros(chrom_sz*4);
    chrom = repmat(chrom, 2, 1);
    chrom = reshape(chrom, chrom_sz*2, chrom_sz);

    C(1:end/2, 1:2:end/2) = chrom;
    C(1:end/2, 2:2:end/2) = chrom;
end;

% Symmetries
C(C(end:-1:1,:)==1) = 1;
C(C(:,end:-1:1)==1) = 1;

end