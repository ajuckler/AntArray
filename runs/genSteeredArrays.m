% ant = AntArray(ones(10), 60500, [], .64, [], 0);
% ant = ant.setName('temp_test');
% ant.genPattern(200, 500, 'YZ', 10);

function genSteeredArrays(name)

    if isempty(name)
        error 'NAME parameter is required';
    elseif iscell(name)
        if length(name) == 1
            name = name{1};
        else
            for el=name
                genSteeredArrays(el);
            end;
            return;
        end;
    end;

    parallel_pool('start');
    ant = AntArray(zeros(10), 60500, [], .84, [], 0);
    ant = ant.setName(name);
    ant = ant.setMin('YZ', -60);
    ant = ant.setMax('YZ', 35);
    
    if strcmpi(name, 'full')
        tmp = ones(10);
        ant = ant.adaptArray(tmp, 200, 0, 0);
    elseif strcmpi(name, 'border_sq_cross')
        tmp = zeros(10);
        tmp([1 end],:) = 1;
        tmp(:,[1 end]) = 1;
        tmp(5:6,5:6) = 1;
        ant = ant.adaptArray(tmp, 200, 0, 0);
        tmp = zeros(10);
        tmp([2:4,7:9],[2:4,7:9])=1;
        ant = ant.adaptArray(tmp, 200, 200, 0);
        tmp = zeros(10);
        tmp(5:6,[2:4,7:9])=1;
        tmp([2:4,7:9],5:6)=1;
        ant = ant.adaptArray(tmp, 200, 200, 200);
    elseif strcmpi(name, 'big_sq_alt')
        tmp = zeros(10);
        tmp(1:5,1:5) = 1;
        tmp(tmp(end:-1:1,end:-1:1)==1) = 1;
        ant = ant.adaptArray(tmp, 200, 200, 200);
        tmp = tmp(:,end:-1:1);
        ant = ant.adaptArray(tmp, 200, -200, 200);
    elseif strcmpi(name, 'alt_v')
        tmp = zeros(10);
        tmp(:,1:2:end) = 1;
        ant = ant.adaptArray(tmp, 200, 200, 0);
        tmp = zeros(10);
        tmp(:,2:2:end) = 1;
        ant = ant.adaptArray(tmp, 200, -200, 0);
    elseif strcmpi(name, 'alt_h')
        tmp = zeros(10);
        tmp(1:2:end,:) = 1;
        ant = ant.adaptArray(tmp, 200, 200, 0);
        tmp = zeros(10);
        tmp(2:2:end,:) = 1;
        ant = ant.adaptArray(tmp, 200, -200, 0);
    elseif strcmpi(name, 'sq_alt')
        tmp = zeros(10);
        for i=1:2:round(length(tmp)/2)
            tmp(i:end-i+1,[i end-i+1]) = 1;
            tmp([i end-i+1],i:end-i+1) = 1;
        end;
        ant = ant.adaptArray(tmp, 200, 200, 0);
        tmp = zeros(10);
        for i=2:2:round(length(tmp)/2)
            tmp(i:end-i+1,[i end-i+1]) = 1;
            tmp([i end-i+1],i:end-i+1) = 1;
        end;
        ant = ant.adaptArray(tmp, 200, -200, 0);
    elseif strcmpi(name, 'alt')
        tmp = zeros(10);
        for i=1:length(tmp)
            if mod(i,2) == 0
                tmp(1:2:end,i) = 1;
            else
                tmp(2:2:end,i)=1;
            end;
        end;
        ant = ant.adaptArray(tmp, 200, 200, 0);
        tmp = zeros(10);
        for i=1:length(tmp)
            if mod(i,2) == 0
                tmp(2:2:end,i) = 1;
            else
                tmp(1:2:end,i)=1;
            end;
        end;
        ant = ant.adaptArray(tmp, 200, -200, 0);
    elseif strcmpi(name, 'alt_diag')
        tmp = zeros(10);
        for i=1:length(tmp)
            if mod(i,2) == 0
                tmp(1:2:end,i) = 1;
            else
                tmp(2:2:end,i)=1;
            end;
        end;
        ant = ant.adaptArray(tmp, 200, 200, -200);
        tmp = zeros(10);
        for i=1:length(tmp)
            if mod(i,2) == 0
                tmp(2:2:end,i) = 1;
            else
                tmp(1:2:end,i)=1;
            end;
        end;
        ant = ant.adaptArray(tmp, 200, -200, 200);
    elseif strcmpi(name, 'X')
        tmp = zeros(10);
        for i=1:length(tmp)
            tmp(i,max(1,i-1):min(length(tmp),i+1)) = 1;
        end;
        tmp(tmp(end:-1:1,:) == 1) = 1;
        ant = ant.adaptArray(tmp, 200, 0, 0);
    elseif strcmpi(name, 'border')
        tmp = zeros(10);
        tmp([1 end], :) = 1;
        tmp(tmp'==1) = 1;
        ant = ant.adaptArray(tmp, 200, 0, 0);
    elseif strcmpi(name, 'border2')
        tmp = zeros(10);
        tmp([1:2 end-1:end], :) = 1;
        tmp(tmp'==1) = 1;
        ant = ant.adaptArray(tmp, 200, 0, 0);
    elseif strcmpi(name, 'border4')
        tmp = zeros(10);
        tmp([1:4 end-3:end], :) = 1;
        tmp(tmp'==1) = 1;
        ant = ant.adaptArray(tmp, 200, 0, 0);
    elseif strcmpi(name, 'border_alt')
        tmp = zeros(10);
        for i=1:2:length(tmp)/2
            j = 2*i-1;
            tmp([1:2 end-1:end], j:j+1) = 1;
        end;
        tmp(tmp'==1) = 1;
        ant = ant.adaptArray(tmp, 200, 0, 0);
        
        tmp = zeros(10);
        for i=1:2:length(tmp)/2-1
            j = 2*i+1;
            tmp([1:2 end-1:end], j:j+1) = 1;
        end;
        tmp(tmp'==1) = 1;
        ant = ant.adaptArray(tmp, 200, 200, 0);
    elseif strcmpi(name, 'corners')
        tmp = zeros(10);
        tmp([1:3 end-2:end], 1) = 1;
        tmp(tmp(:,end:-1:1) == 1) = 1;
        tmp(tmp' == 1) = 1;
        ant = ant.adaptArray(tmp, 200, 0, 0);
    elseif strcmpi(name, 'corner_squares')
        tmp = zeros(10);
        tmp([2:4 7:9], [2:4 7:9]) = 1;
        ant = ant.adaptArray(tmp, 200, 0, 0);
    elseif strcmpi(name, '+')
        tmp = zeros(10);
        tmp(:,5:6) = 1;
        tmp(tmp'==1) = 1;
        ant = ant.adaptArray(tmp, 200, 0, 0);
    elseif strcmpi(name, '+_centre')
        tmp = zeros(10);
        tmp(2:end-1,5:6) = 1;
        tmp(tmp'==1) = 1;
        tmp(5:6,5:6) = 0;
        ant = ant.adaptArray(tmp, 200, 0, 0);
    elseif strcmpi(name, 'intl_sq')
        tmp = zeros(20);
%         tmp([5 16], 5:16) = 1;
        for i=1:10
            tmp(i, 10-i+1) = 1;
        end;
%         for i=3:10
%             tmp(i, 10-i+3) = 1;
%         end;
        tmp(tmp(end:-1:1, :) == 1) = 1;
        tmp(tmp' == 1) = 1;
        tmp(tmp(:, end:-1:1) == 1) = 1;
        ant = ant.rstArray(zeros(20));
        ant = ant.adaptArray(tmp, 200, 0, 0);
    else
        error 'Unknown option';
    end;
    ant.plotAntArray();
    ant.genPattern(200, 1000, 'YZ', 5);

end