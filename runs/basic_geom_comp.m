function basic_geom_comp(name)

    if isempty(name)
        error 'NAME parameter is required';
    elseif iscell(name)
        if length(name) == 1
            name = name{1};
        else
            for el=name
                basic_geom_comp(el);
            end;
            return;
        end;
    end;

    parallel_pool('start');
    ant = AntArray(zeros(64), 60500, [], .84, [], 0);
    ant = ant.setName(name);
    ant = ant.setMin('YZ', -60);
    ant = ant.setMax('YZ', 20);
    ant = ant.setMin('XY', -60);
    ant = ant.setMax('XY', 20);
    
    if strcmpi(name, 'sq2')
        for rem_els = [4 8 12];
            tmp = gen_sq_arrgt(64, rem_els);
            ant = ant.adaptArray(tmp, 100000, 0, 0);
            ant = ant.setName([name '_' num2str(rem_els)]);

            ant.plotAntArray();
            ant = ant.genPattern(5000, 3000, 'YZ', 15);
            ant = ant.genPattern(5000, 3000, 'YZ-BW', 15);
            ant = ant.genPattern(10000, 3000, 'XY', 15);
            ant = ant.genPattern(10000, 3000, 'XY-BW', 15);
            ant = ant.genPattern(10000, 3000, 'theta', 15, pi/4);
            ant = ant.genPattern(10000, 3000, 'theta-BW', 15, pi/4);
        end;
    elseif strcmpi(name, 'tr') || strcmpi(name, 'tr2')
        for rem_els = [4 8 12]
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

            tmp(tmp(:,end:-1:1)==1) = 1;
            tmp(tmp(end:-1:1,:)==1) = 1;
            if strcmpi(name, 'tr2')
                tmp(tmp'==1) = 1;
            end;
            
            ant = ant.rstArray(zeros(60));
            ant = ant.adaptArray(tmp, 100000, 0, 0);
            ant = ant.setName([name '_' num2str(rem_els)]);

            ant.plotAntArray();
            ant = ant.genPattern(5000, 3000, 'YZ', 15);
            ant = ant.genPattern(5000, 3000, 'YZ-BW', 15);
            ant = ant.genPattern(10000, 3000, 'XY', 15);
            ant = ant.genPattern(10000, 3000, 'XY-BW', 15);
            if strcmpi(name, 'tr')
                ant = ant.genPattern(10000, 3000, 'theta', 15, 0);
                ant = ant.genPattern(10000, 3000, 'theta-BW', 15, 0);
            end;
        end;
    elseif strcmpi(name, 'circ')
        for rem_els = [9 18 27]
            tmp = drawCircle(64, 32, rem_els);
            ant = ant.adaptArray(tmp, 100000, 0, 0);
            ant = ant.setName([name '_' num2str(rem_els)]);

            ant.plotAntArray();
            ant = ant.genPattern(5000, 3000, 'YZ', 15);
            ant = ant.genPattern(5000, 3000, 'YZ-BW', 15);
            ant = ant.genPattern(10000, 3000, 'XY', 15);
            ant = ant.genPattern(10000, 3000, 'XY-BW', 15);
        end;
    elseif strcmpi(name, 'sq2_best')
        rem_els = 9;
        tmp = gen_sq_arrgt(64, rem_els);
        ant = ant.adaptArray(tmp, 100000, 0, 0);
        
        el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100, 3);
        ant = ant.setComments(sprintf(['\\# of elements: ' el_ratio '\\%%']));

        ant.plotAntArray();
        ant = ant.genPattern(5000, 3000, 'YZ', 15);
        ant = ant.genPattern(5000, 3000, 'YZ-BW', 15);
        ant = ant.genPattern(10000, 3000, 'XY', 15);
        ant = ant.genPattern(10000, 3000, 'XY-BW', 15);
        ant = ant.genPattern(10000, 3000, 'theta', 15, pi/4);
        ant.genPattern(10000, 3000, 'theta-BW', 15, pi/4);
    elseif strcmpi(name, 'tr_best') || strcmpi(name, 'tr2_best')
        rem_els = 9;
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

        tmp(tmp(:,end:-1:1)==1) = 1;
        tmp(tmp(end:-1:1,:)==1) = 1;
        if strcmpi(name, 'tr2_best')
            tmp(tmp'==1) = 1;
        end;

        ant = ant.rstArray(zeros(60));
        ant = ant.adaptArray(tmp, 100000, 0, 0);

        el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100, 3);
        ant = ant.setComments(sprintf(['\\# of elements: ' el_ratio '\\%%']));

        ant.plotAntArray();
        ant = ant.genPattern(5000, 3000, 'YZ', 15);
        ant = ant.genPattern(5000, 3000, 'YZ-BW', 15);
        ant = ant.genPattern(10000, 3000, 'XY', 15);
        ant.genPattern(10000, 3000, 'XY-BW', 15);
        if strcmpi(name, 'tr_best')
            ant = ant.genPattern(10000, 3000, 'theta', 15, 0);
            ant.genPattern(10000, 3000, 'theta-BW', 15, 0);
        end;
    elseif strcmpi(name, 'circ_best')
        rem_els = 22;
        tmp = drawCircle(64, 32, rem_els);
        ant = ant.adaptArray(tmp, 100000, 0, 0);
        
        el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100, 3);
        ant = ant.setComments(sprintf(['\\# of elements: ' el_ratio '\\%%']));

        ant.plotAntArray();
        ant = ant.genPattern(5000, 3000, 'YZ', 15);
        ant = ant.genPattern(5000, 3000, 'YZ-BW', 15);
        ant = ant.genPattern(10000, 3000, 'XY', 15);
        ant.genPattern(10000, 3000, 'XY-BW', 15);
    elseif strcmpi(name, 'tr_spacing')
        intv = .5:.1:.7;
        wArr = zeros(length(intv), 2);
        for it = 1:length(intv)
            spacing = intv(it);
            ant = AntArray(zeros(60), 60500, [], spacing, [], 0);
            ant = ant.setName([name '_sp_' num2str(spacing*100)]);
            ant = ant.setMin('YZ', -60);
            ant = ant.setMax('YZ', 20);
            ant = ant.setMin('XY', -60);
            ant = ant.setMax('XY', 20);
            
            rem_els = 9;
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

            tmp(tmp(:,end:-1:1)==1) = 1;
            tmp(tmp(end:-1:1,:)==1) = 1;
            
            ant = ant.adaptArray(tmp, 100000, 0, 0);
            ant = ant.setComments(sprintf(['Spacing: ' num2str(spacing) '$\\lambda$']));

            ant.plotAntArray();
            ant = ant.genPattern(5000, 3000, 'YZ', 15);
            ant = ant.genPattern(5000, 3000, 'YZ-BW', 15);
            ant = ant.genPattern(10000, 3000, 'XY', 15);
            ant = ant.genPattern(10000, 3000, 'XY-BW', 15);
            wArr(it, 1) = ant.weight('XY');
            ant = ant.genPattern(10000, 3000, 'theta', 15, 0);
            ant.genPattern(10000, 3000, 'theta-BW', 15, 0);
            wArr(it, 2) = ant.weight('theta', 0);
        end;
        export_dat(wArr, ['weights_' name '_' datestr(now, 'yyyymmdd')]);
    elseif strcmpi(name, 'tr_freq')
        intv = 59.4:.1:61.6;
        wArr = zeros(length(intv), 2);
        for it = 1:length(intv)
            freq = intv(it);
            ant = AntArray(zeros(60), freq*1000, [], .84, [], 0);
            ant = ant.setNormFreq(60500);
            ant = ant.setName([name '_fr_' num2str(freq*10)]);
            ant = ant.setMin('YZ', -60);
            ant = ant.setMax('YZ', 20);
            ant = ant.setMin('XY', -60);
            ant = ant.setMax('XY', 20);
            
            rem_els = 9;
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

            tmp(tmp(:,end:-1:1)==1) = 1;
            tmp(tmp(end:-1:1,:)==1) = 1;
            
            ant = ant.adaptArray(tmp, 100000, 0, 0);
            ant = ant.setComments(sprintf(['Frequency: ' num2str(freq) '\,GHz']));

            ant.plotAntArray();
            ant = ant.genPattern(5000, 3000, 'YZ', 15);
            ant = ant.genPattern(5000, 3000, 'YZ-BW', 15);
            ant = ant.genPattern(10000, 3000, 'XY', 15);
            ant = ant.genPattern(10000, 3000, 'XY-BW', 15);
            wArr(it, 1) = ant.weight('XY');
            ant = ant.genPattern(10000, 3000, 'theta', 15, 0);
            ant.genPattern(10000, 3000, 'theta-BW', 15, 0);
            wArr(it, 2) = ant.weight('theta', 0);
        end;
        export_dat(wArr, ['weights_' name '_' datestr(now, 'yyyymmdd')]);
    elseif strcmpi(name, 'tr_quant')
        intv = 1:4;
        wArr = zeros(length(intv), 2);
        for it = 1:length(intv)
            lvl = intv(it);
            ant = AntArray(zeros(60), 60500, [], .84, [], 0);
            ant = ant.setName([name '_quant_' num2str(lvl)]);
            ant = ant.setMin('YZ', -60);
            ant = ant.setMax('YZ', 20);
            ant = ant.setMin('XY', -60);
            ant = ant.setMax('XY', 20);
            
            rem_els = 9;
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

            tmp(tmp(:,end:-1:1)==1) = 1;
            tmp(tmp(end:-1:1,:)==1) = 1;
            
            ant = ant.adaptArray(AntArray.quantize(tmp, 2, lvl), 100000, 0, 0);
            ant = ant.setComments(sprintf(['Quantization level: ' num2str(lvl)]));

            ant.plotAntArray();
            ant = ant.genPattern(5000, 3000, 'YZ', 15);
            ant = ant.genPattern(5000, 3000, 'YZ-BW', 15);
            ant = ant.genPattern(10000, 3000, 'XY', 15);
            ant = ant.genPattern(10000, 3000, 'XY-BW', 15);
            wArr(it, 1) = ant.weight('XY');
            ant = ant.genPattern(10000, 3000, 'theta', 15, 0);
            ant.genPattern(10000, 3000, 'theta-BW', 15, 0);
            wArr(it, 2) = ant.weight('theta', 0);
        end;
        export_dat(wArr, ['weights_' name '_' datestr(now, 'yyyymmdd')]);
    else
        error 'Unknown option';
    end;

end