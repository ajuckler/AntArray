clear

matlabpool open

% Squares
for rem_els=16:15
    % Init
    arr = AntArray(zeros(64), 60500, [], 0.84);
    arr = arr.setName(['sq_' mat2str(rem_els)]);
    arr = arr.setMax('XY', 30);
    arr = arr.setMax('YZ', 30);

    % Create elements' pattern
    tmp = zeros(64);
    for i=16:32-rem_els
        tmp(i, i:end-i+1)=1;
    end;

    tmp(tmp(end:-1:1,:)==1)=1;
    tmp(tmp'==1)=1;

    for j=1:16-rem_els
        for i=1:32
            ln = i+2*(j-1);
            rw = 32-i+1;
            if ln <= 32 && rw <=32
                tmp(ln,rw)=1;
            end;
        end;
    end;

    tmp(tmp'==1)=1;
    tmp(tmp(end:-1:1,:)==1)=1;
    tmp(tmp(:,end:-1:1)==1)=1;

    arr = arr.adaptArray(tmp, 90000, 0, 0);
    
    el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
    
    arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
        'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%\n' ...
        'Adjusted impedance computations']));

    % Plots
    arr.plotAntArray();
    if rem_els ~= 13
        arr = arr.genPattern(11000, 3000, 'XY', 30);
    end;
    arr = arr.genPattern([], [], 'XY-BW');
    arr = arr.E_strength(15000, 0, 0, 500);
    arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
    arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
end;

% Triangles
for rem_els=0:15
    % Init
    arr = AntArray(zeros(60), 60500, [], 0.84);
    arr = arr.setName(['tr_' num2str(rem_els)]);
    arr = arr.setMax('XY', 50);
    arr = arr.setMax('YZ', 50);

    % Create elements' pattern
    tmp = zeros(60);
    len = size(arr.M,1)/4;
    d_ratio = abs(1-tan(pi/3))/sqrt((tan(pi/3))^2-1);
    k_max = 30 - round(rem_els/d_ratio);
    k_lim = round(k_max*d_ratio);
    for k=1:k_max
        i_max = size(arr.M,1)-len+2-k;
        if k <= k_lim
            tmp(len-1+k, 3+round(k/tan(pi/3)):end-2-round(k/tan(pi/3)))=1;
        end;
        for i=min(k,max(k_lim - 1,1)):i_max
           tmp(len-1+i, 3+round((i+k-1)/tan(pi/3)))=1; 
        end
    end;

    tmp(tmp(end:-1:1,:)==1)=1;
    tmp(tmp(:,end:-1:1)==1)=1;

    arr = arr.adaptArray(tmp, 90000, 0, 0);
    
    el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
    
    arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
        'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));

    % Plots
    arr.plotAntArray();
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    arr = arr.E_strength(15000, 0, 0, 500);
    arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
    arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
end;

% Triangles 2
for rem_els=0:15
    % Init
    arr = AntArray(zeros(60), 60500, [], 0.84);
    arr = arr.setName(['tr2_' num2str(rem_els)]);
    arr = arr.setMax('XY', 50);
    arr = arr.setMax('YZ', 50);

    % Create elements' pattern
    tmp = zeros(60);
    len = size(arr.M,1)/4;
    d_ratio = abs(1-tan(pi/3))/sqrt((tan(pi/3))^2-1);
    k_max = 30 - round(rem_els/d_ratio);
    k_lim = round(k_max*d_ratio);
    for k=1:k_max
        i_max = size(arr.M,1)-len+2-k;
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

    arr = arr.adaptArray(tmp, 90000, 0, 0);
    
    el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
    
    arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
        'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));

    % Plots
    arr.plotAntArray();
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    arr = arr.E_strength(15000, 0, 0, 500);
    arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
    arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
end;

matlabpool close