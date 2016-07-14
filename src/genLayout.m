%GENLAYOUT Create random grid matrices
%
%   ptrns = GENLAYOUT(side, num)
%   INPUT:
%       side:   side length of the matrix
%       num:    number of matrices to generate
%   OUTPUT:
%       ptrns:  cell array containing the generated matrices
%
%   See also GA_2D CHROM2MAT MAT2CHROM

%   Copyright 2016, Antoine Juckler. All rights reserved.


function ptrns = genLayout(side, nb)

    if nargin < 2 || isempty(nb)
        nb = 1;
    end;
    if side < 1 || nb < 1
        error('MyERR:Invalid', 'Invalid input argument');
    end;
    
    ptrns = cell(1, nb);
    
    side1 = side/2;
    side2 = side/4;
    
    for i=1:nb
        tmp = zeros(side);
        
        ypos = randi([1 side2]);
        while ypos <= side
            xpos = randi([1 side1]);
            while xpos <= side
                width = min(randi([1 side2]), side-xpos+1);
                height = min(randi([1 side2]), side-ypos+1);

                tmp(ypos:ypos+height-1, xpos:xpos+width-1)=1;
                xpos = xpos+width+randi([1 side1]);
            end;
            ypos = ypos+randi([1 side2]);
        end;
        
        ptrns{i} = tmp;
    end;
end