function gen_init_pop(nb, model)
if nargin < 1 || isempty(nb)
    nb = 64;
end
if nb == 64 && (nargin < 2 || isempty(model))
    model = 1;
elseif ~isempty(model)
    model = model > 0;
end;

if nb == 64
    if model == 1
        rem = [9 4 9 22];
    else
        rem = [9 9 9 22];
    end;
elseif nb == 32
    rem = [0 4 5 8];
elseif nb == 16
    rem = [0 0 0 0];
elseif nb == 8
    rem = [0 0 0 0];
else
    error 'Invalid input argument';
end;

% Squares 2
rem_els = rem(1);
    arr = AntArray(zeros(nb), 60500, [], 0.84, [], 0);
    arr = arr.setName(['sq2_' num2str(nb) '_init']);

    % Create elements' pattern
    tmp = gen_sq_arrgt(nb, rem_els);

    arr = arr.adaptArray(tmp, 100000, 0, 0);
    arr.plotAntArray();

% Triangles
rem_els = rem(2);
    arr = AntArray(zeros(nb), 60500, [], 0.84, [], 0);
    if nb == 64 && model == 0
        arr = arr.setName(['tr_' num2str(nb) '_60']);
    else
        arr = arr.setName(['tr_' num2str(nb) '_init']);
    end;

    % Create elements' pattern
    if nb == 64 && model == 0
        tmp = zeros(60);
        len = size(tmp,1)/4;
        d_ratio = abs(1-tan(pi/3))/sqrt((tan(pi/3))^2-1);
        k_max = size(tmp,1)/2 - round(rem_els/d_ratio);
        k_lim = round(k_max*d_ratio);
        for k=1:k_max
            i_max = size(tmp,1)-len+2-k;
            if k <= k_lim
                tmp(len-1+k, 3+round(k/tan(pi/3)):end-2-round(k/tan(pi/3)))=1;
            end;
            for i=min(k,max(k_lim - 1,1)):i_max
               tmp(len-1+i, 3+round((i+k-1)/tan(pi/3)))=1; 
            end
        end;

        tmp(tmp(end:-1:1,:)==1)=1;
        tmp(tmp(:,end:-1:1)==1)=1;
        mat = zeros(64);
        mat(3:end-2,3:end-2) = tmp;
    else
        mat = gen_tr_arrgt(nb, 2, rem_els);
    end;

    arr = arr.adaptArray(mat, 100000, 0, 0);

    arr.plotAntArray();

% Triangles 2
rem_els = rem(3);
    arr = AntArray(zeros(nb), 60500, [], 0.84, [], 0);
    if nb==64 && model == 0
        arr = arr.setName(['tr2_' num2str(nb) '_60']);
    else
        arr = arr.setName(['tr2_' num2str(nb) '_init']);
    end;

    % Create elements' pattern
    if nb == 64 && model == 0
        tmp = zeros(60);
        len = size(tmp,1)/4;
        d_ratio = abs(1-tan(pi/3))/sqrt((tan(pi/3))^2-1);
        k_max = size(tmp,1)/2 - round(rem_els/d_ratio);
        k_lim = round(k_max*d_ratio);
        for k=1:k_max
            i_max = size(tmp,1)-len+2-k;
            if k <= k_lim
                tmp(len-1+k, 3+round(k/tan(pi/3)):end-2-round(k/tan(pi/3)))=1;
            end;
            for i=min(k,max(k_lim - 1,1)):i_max
               tmp(len-1+i, 3+round((i+k-1)/tan(pi/3)))=1; 
            end
        end;

        tmp(tmp'==1)=1;
        tmp(tmp(end:-1:1,:)==1)=1;
        tmp(tmp(:,end:-1:1)==1)=1;
        mat = zeros(64);
        mat(3:end-2,3:end-2) = tmp;
    else
        mat = gen_tr_arrgt(nb, 4, rem_els);
    end;

    arr = arr.adaptArray(mat, 100000, 0, 0);
    arr.plotAntArray();


% Circles
rem_els = rem(4);
    arr = AntArray(zeros(nb), 60500, [], 0.84, [], 0);
    arr = arr.setName(['circ_' num2str(nb) '_init']);

    % Create elements' pattern
    tmp = drawCircle(nb, nb/2, rem_els);
    arr = arr.adaptArray(tmp, 100000, 0, 0);
    arr.plotAntArray();