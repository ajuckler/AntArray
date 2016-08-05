%GENMASK generates a mask based on binary matrix with given radius
%
%   This function generates a mask based on a chromosome: the chromosome is
%   first converted to its 2D representation, a 2D mask with given radius
%   is then generated, which is finally converted back to a 1D vector.
%
%   The radius is the number of cells around an element that should be
%   included in the mask
%
%   mask = genMask(chrom, rad)
%
%   INPUT:
%       chrom:  chromosome (1D binary vector)
%       rad:    (optional) radius (number of cells) [default=1]
%   OUTPUT:
%       mask:   mask (1D binary vector)
%
%   See also GA_2D LOCAL_OPT

%   Copyright 2016, Antoine Juckler. All rights reserved.

function mask = genMask(chrom, rad)
    if nargin < 2 || isempty(rad)
        rad = 1;
    end;
    dim = sqrt(length(chrom));
    mat = reshape(chrom, dim, dim);
    mask = mat;

    for i=1:dim
        for j=1:dim
            if mat(i,j)
                lim_l = max(1, i-rad);
                lim_r = min(dim, i+rad);
                lim_u = max(1, j-rad);
                lim_b = min(dim, j+rad);
                mask(lim_l:lim_r, lim_u:lim_b) = 1;
                
%                 if i > 1
%                     mask(i-1,j) = 1;
%                 end;
%                 if i < dim
%                     mask(i+1,j) = 1;
%                 end;
%                 if j > 1
%                     mask(i,j-1) = 1;
%                 end;
%                 if j < dim
%                     mask(i,j+1) = 1;
%                 end;
            end;
        end;
    end;

    mask = reshape(mask, 1, numel(mask));
end