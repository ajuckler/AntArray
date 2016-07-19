% ant = AntArray(ones(10), 60500, [], .64, [], 0);
% ant = ant.setName('temp_test');
% ant.genPattern(200, 500, 'YZ', 10);

function genSteeredArrays(name)

    if isempty(name)
        error 'NAME parameter is required';
    end;

    parallel_pool('start');
    ant = AntArray(zeros(10), 60500, [], .84, [], 0);
    ant = ant.setName(name);

    if strcmpi(name, 'border_sq_cross')
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
    else
        error 'Unknown option';
    end;
    ant.plotAntArray();
    ant.genPattern(200, 1000, 'YZ', 10);

end