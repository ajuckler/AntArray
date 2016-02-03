clear

rem_els = 9;
sq_dim = 2;

if verLessThan('matlab','8.2')
    matlabpool open
else
    poolobj = parpool;
end;

for quant_lvl=4:-1:1
    % Init
    arr = AntArray(zeros(60), 60500, [], 0.84);
    arr = arr.setName(['tr9_quant' num2str(quant_lvl)]);
    arr = arr.setMax('XY', 30);
    arr = arr.setMax('YZ', 30);
    arr = arr.setMax('E', 25);
    arr = arr.setMin('XY', -60);
    arr = arr.setMin('YZ', -60);
    arr = arr.setMin('E', -15);

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
    
    % Quantization
    if rem(size(tmp,1),sq_dim) > 10^-6
        error('Dimension of quantization matrix should be a divider of the dimensions of the element''s matrix');
    elseif quant_lvl > sq_dim^2
        error('Quantization level must not exceed the number of elements in the quantization matrix');
    end;
    half_dim = round(size(tmp,1)/sq_dim);
    
    for i=1:half_dim
        y = (i-1)*sq_dim + 1;
        for j=1:half_dim
            x = (j-1)*sq_dim + 1;
            sq = tmp(y:y+sq_dim-1, x:x+sq_dim-1);
            if sum(sum(sq)) < quant_lvl
                sq = zeros(sq_dim);
            else
                sq = ones(sq_dim);
            end;
            tmp(y:y+sq_dim-1, x:x+sq_dim-1) = sq(:,:);
        end;
    end;
            

    arr = arr.adaptArray(tmp, 90000, 0, 0);
    
    el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
    
    arr = arr.setComments(sprintf(['Elements spacing: 0.84$\\lambda$\n' ...
        'Quantization level:' num2str(quant_lvl) '\n' ...
        '\\# of elements: ' el_ratio '\\%%']));

    % Plots
    arr.plotAntArray(1);
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    arr = arr.E_strength(15000, 0, 0, 500);
%     arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
%     arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
    arr = arr.genPattern(10*1000, 3000, 'YZ', 30);
    arr = arr.genPattern(10*1000, [], 'YZ-BW');
end;

if verLessThan('matlab','8.2')
    matlabpool close
else
   delete(poolobj);
end;