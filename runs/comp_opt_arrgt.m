function comp_opt_arrgt(opt, type)
    % opt is optimal ANTARRAY
    
    if strcmpi(type, 'spacing')
        intv = .5:.1:.7;
        wArr = zeros(length(intv), 2);
        for it = 1:length(intv)
            spacing = intv(it);
            ant = AntArray(opt.M, 60500, [], spacing, [], 0);
            ant = ant.setName(['sp_' num2str(spacing*100)]);
            ant = ant.setMin('YZ', -60);
            ant = ant.setMax('YZ', 20);
            ant = ant.setMin('XY', -60);
            ant = ant.setMax('XY', 20);
            
            ant = ant.setComments(sprintf(['Spacing: ' num2str(spacing) '$\\lambda$']));

            ant.plotAntArray();
            ant = ant.genPattern(2000, 3000, 'YZ', 15);
            ant = ant.genPattern(2000, 3000, 'YZ-BW', 15);
            ant = ant.genPattern(11000, 3000, 'XY', 15);
            ant = ant.genPattern(11000, 3000, 'XY-BW', 15);
            wArr(it, 1) = ant.weight('XY');
            ant = ant.genPattern(11000, 3000, 'theta', 15, 0);
            ant.genPattern(11000, 3000, 'theta-BW', 15, 0);
            wArr(it, 2) = ant.weight('theta', 0);
        end;
        export_dat(wArr, ['weights_' type '_' datestr(now, 'yyyymmdd')]);
    elseif strcmpi(type, 'freq')
        intv = 59.4:.1:61.6;
        wArr = zeros(length(intv), 2);
        for it = 1:length(intv)
            freq = intv(it);
            ant = AntArray(opt.M, freq*1000, [], .84, [], 0);
            ant = ant.setNormFreq(60500);
            ant = ant.setName(['fr_' num2str(freq*10)]);
            ant = ant.setMin('YZ', -60);
            ant = ant.setMax('YZ', 20);
            ant = ant.setMin('XY', -60);
            ant = ant.setMax('XY', 20);
            
            ant = ant.setComments(sprintf(['Frequency: ' num2str(freq) '\,GHz']));

            ant.plotAntArray();
            ant = ant.genPattern(2000, 3000, 'YZ', 15);
            ant = ant.genPattern(2000, 3000, 'YZ-BW', 15);
            ant = ant.genPattern(11000, 3000, 'XY', 15);
            ant = ant.genPattern(11000, 3000, 'XY-BW', 15);
            wArr(it, 1) = ant.weight('XY');
            ant = ant.genPattern(11000, 3000, 'theta', 15, 0);
            ant.genPattern(11000, 3000, 'theta-BW', 15, 0);
            wArr(it, 2) = ant.weight('theta', 0);
        end;
        export_dat(wArr, ['weights_' type '_' datestr(now, 'yyyymmdd')]);
    elseif strcmpi(type, 'quant')
        intv = 1:4;
        wArr = zeros(length(intv), 2);
        for it = 1:length(intv)
            lvl = intv(it);
            mat = AntArray.quantize(opt.M, 2, lvl);
            ant = AntArray(mat, 60500, [], .84, [], 0);
            ant = ant.setName(['quant_' num2str(lvl)]);
            ant = ant.setMin('YZ', -60);
            ant = ant.setMax('YZ', 20);
            ant = ant.setMin('XY', -60);
            ant = ant.setMax('XY', 20);
            
            ant = ant.setComments(sprintf(['Quantization level: ' num2str(lvl)]));

            ant.plotAntArray();
            ant = ant.genPattern(2000, 3000, 'YZ', 15);
            ant = ant.genPattern(2000, 3000, 'YZ-BW', 15);
            ant = ant.genPattern(11000, 3000, 'XY', 15);
            ant = ant.genPattern(11000, 3000, 'XY-BW', 15);
            wArr(it, 1) = ant.weight('XY');
            ant = ant.genPattern(11000, 3000, 'theta', 15, 0);
            ant.genPattern(11000, 3000, 'theta-BW', 15, 0);
            wArr(it, 2) = ant.weight('theta', 0);
        end;
        export_dat(wArr, ['weights_' type '_' datestr(now, 'yyyymmdd')]);
    else
        error 'Unknown option';
    end;
end