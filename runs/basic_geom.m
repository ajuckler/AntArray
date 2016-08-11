function basic_geom(name)

    if isempty(name)
        error 'NAME parameter is required';
    elseif iscell(name)
        if length(name) == 1
            name = name{1};
        else
            for el=name
                basic_geom(el);
            end;
            return;
        end;
    end;

    parallel_pool('start');
    ant = AntArray(zeros(64), 60500, [], .84, [], 0);
    ant = ant.setName(name);
    ant = ant.setMin('YZ', -60);
    ant = ant.setMax('YZ', 20);
    
    if strcmpi(name, 'full')
        tmp = ones(64);
        ant = ant.adaptArray(tmp, 5000, 0, 0);
    elseif strcmpi(name, 'sq2')
        tmp = gen_sq_arrgt(64, 5);
        ant = ant.adaptArray(tmp, 5000, 0, 0);
    elseif strcmpi(name, 'sq_wrong_spacing')
        arr = zeros(64);
        for i=16:26
            arr(i, i:end-i+1)=1;
        end;
        for j=1:11
            for i=j:32
                arr(i, 32-i+j) = 1;
            end;
        end;

        arr(arr(end:-1:1,:)==1)=1;
        arr(arr(:,end:-1:1)==1)=1;
        arr(arr'==1)=1;
        
        ant = ant.adaptArray(arr, 5000, 0, 0);
    elseif strcmpi(name, 'tr') || strcmpi(name, 'tr2') || strcmpi(name, 'tr4')
        tmp = zeros(60);
        rem_els = 5;
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
        
        tmp(tmp(:,end:-1:1)==1) = 1;
        tmp = tmp(end:-1:1,:);
        
        if strcmpi(name, 'tr2')
            tmp(tmp(end:-1:1,:)==1) = 1;
        elseif strcmpi(name, 'tr4')
            tmp(tmp(end:-1:1,:)==1) = 1;
            tmp(tmp'==1) = 1;
        end;
        
        ant = ant.rstArray(zeros(60));
        ant = ant.adaptArray(tmp, 5000, 0, 0);   
    elseif strcmpi(name, 'circ')
        tmp = drawCircle(64, 32, 7);
        ant = ant.adaptArray(tmp, 5000, 0, 0);
    else
        error 'Unknown option';
    end;
    ant.plotAntArray();
    ant = ant.genPattern(5000, 3000, 'YZ', 15);
    ant = ant.genPattern(5000, 3000, 'YZ-BW', 15);
    ant = ant.genPattern(10000, 3000, 'XY', 30);
    ant = ant.genPattern(10000, 3000, 'XY-BW', 30);
    ant = ant.genPattern(10000, 3000, 'theta', 30, pi/4);
    ant = ant.genPattern(10000, 3000, 'theta-BW', 30, pi/4);
end