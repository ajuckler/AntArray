function run_GA(start_pop, dist, mode, quant)

    if ~exist('dist', 'var') || isempty(dist)
        dist = [2000, 5000, 10000];
    end;
    if ~exist('mode', 'var') || isempty(mode)
        mode = 0:1;
    else
        mode = (mode > 0);
    end;
    if ~exist('quant', 'var') || isempty(quant)
        quant = 0:1;
    else
        quant = (quant > 0);
    end;

    for disti = dist
        for modei = mode
            for quanti = quant
                ant = ga_2D(disti, start_pop, quanti, modei);
                parallel_pool('start');
                ant.plotAntArray();
                ant = ant.genPattern(disti, 3000, 'YZ', 30);
                ant = ant.genPattern(disti, 3000, 'YZ-BW');
                ant = ant.genPattern(disti, 3000, 'YZ-main', 30);
                ant = ant.genPattern(11000, 3000, 'XY', 30);
                ant = ant.genPattern(11000, 3000, 'XY-BW', 30);
                ant = ant.E_strength(11000);
            end;
        end;
    end;
    parallel_pool('close');
end