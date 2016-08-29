function arr = gen_tr_arrgt(dim, mode, rem)
    %GEN_TR_ARRGT generates the two/four triangle arrangement with specified
    %dimension
    %
    % ARR = GEN_TR_ARRGT(DIM, MODE, REM)
    %
    % INPUT:
    %   DIM:    dimension of the array
    %   MODE:   number of triangles (2 or 4)
    %   REM:    (optional) number of lines to remove [default = 0]
    % OUTPUT:
    %   ARR:    the generated binary array
    %

    % Copyright 2015-2016, Antoine JUCKLER. All rights reserved
    
    if nargin < 3
        rem = 0;
    else
        rem = abs(round(rem));
    end;
    if mode ~= 2 && mode ~= 4
        error 'Invalid MODE argument';
    end;
    
    b_fl = dim/(sqrt(3)+1/sqrt(3))*2;
    b = 2*round(b_fl/2);
    h = round(b_fl/2*sqrt(3));
    
    rem_els = max(1, ceil(h/3 - rem));

    arr = zeros(dim);
    mask = arr;
    for i=1:h
        x = dim/2 - round((i-1)*tan(pi/6));
        arr(i, x) = 1;
        mask(i, x:end) = 1;
    end
    mask(h+1:end,:) = 0;
    for k=1:rem_els
        arr(h-k+1, dim/2:-1:dim/2-b/2+1) = 1;
    end;

    k_max = ceil(rem_els*tan(pi/3));
    for k=2:k_max
        for i=1:h
            x = dim/2 - round((i-1)*tan(pi/6));
            if i+k-1 < h
                arr(i+k-1, x) = 1;
            end;
        end;
    end;

    arr(mask == 0) = 0;

    arr(arr(end:-1:1,:)==1)=1;
    arr(arr(:,end:-1:1)==1)=1;
    
    if mode == 4
        arr(arr'==1)=1;
    end;

end