%GEN_SQ_ARRGT generates the double square arrangement with specified
%dimension
%
%   ARR = GEN_SQ_ARRGT(DIM, REM)
%   INPUT:
%       DIM:    dimension of the array
%       REM:    (optional) number of lines to remove [default = 0]
%   OUTPUT:
%       ARR:    the generated binary array
%

% Copyright 2016, Antoine Juckler. All rights reserved

function arr = gen_sq_arrgt(dim, rem_els)
    if nargin < 2
        rem_els = 0;
    else
        rem_els = abs(round(rem_els));
    end;
    
    if mod(dim,2) ~= 0
        error 'DIM should be an even number';
    end;
    
    arr = zeros(dim);
    for i=round(dim/4):dim/2-1-rem_els
        arr(i, i:end-i+1)=1;
    end;

    arr(arr(end:-1:1,:)==1)=1;
    arr(arr'==1)=1;

    for j=1:round(dim/4)-rem_els
        for i=1:dim/2
            ln = i+2*(j-1);
            rw = dim/2-i+1;
            if ln <= dim/2 && rw <= dim/2
                arr(ln,rw)=1;
            end;
        end;
    end;

    arr(arr'==1)=1;
    arr(arr(end:-1:1,:)==1)=1;
    arr(arr(:,end:-1:1)==1)=1;
end