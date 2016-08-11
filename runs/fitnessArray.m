function fitnessArray(name, sz)
    
    if strcmpi(name, 'freq');
        intv = 594:616;
        opt = 'fr';
    elseif strcmpi(name, 'spacing');
%         intv = [.5:.1:.6 .7:.02:1]*100;
%         intv = [.7:.02:.8]*100;
        intv = (.7:.02:1)*100;
        opt = 'sp';
    elseif strcmpi(name, 'quant');
        intv = 1:4;
        opt = 'quant';
    else
        error 'Invalid name';
    end;

    wArr = zeros(length(intv), 2);
    for i=1:length(intv)
        it = intv(i);
        fname = ['dat/pattern_i' num2str(sz) '_0_' opt '_' num2str(it) '_2.dat'];
        wArr(i, :) = [it compFitness(fname)];
%         genBW(['i' num2str(sz) '_0_' opt '_' num2str(it)], 'theta');
%         genBW(['i' num2str(sz) '_0_' opt '_' num2str(it)], 'XY');
%         genBW(['i' num2str(sz) '_0_' opt '_' num2str(it)], 'YZ');
%         genColor(['i' num2str(sz) '_0_' opt '_' num2str(it)], 'theta');
%         genColor(['i' num2str(sz) '_0_' opt '_' num2str(it)], 'XY');
%         genColor(['i' num2str(sz) '_0_' opt '_' num2str(it)], 'YZ');
    end;
    export_dat(wArr, ['weights_' name '_' datestr(now, 'yyyymmdd')]);
    
    fig = figure();
    if strcmpi(name, 'freq')
        intv = intv/10;
        plot(intv, wArr(:,2), 'LineWidth', 2);
        hold on; box on;
        xlabel('Frequency [GHz]', 'Interpreter', 'latex', 'FontSize', 22);
        ylabel('Fitness [Vm]', 'Interpreter', 'latex', 'FontSize', 22);
        xlim([intv(1) intv(end)]);
        set(get(fig, 'CurrentAxes'), 'FontSize', 16);

        print_plots(fig, ['i' num2str(sz) '_0_' opt '_fit']);
    elseif strcmpi(name, 'spacing')
        intv = intv/100;
        plot(intv, wArr(:,2), 'LineWidth', 2);
        hold on; box on;
        xlabel('Normalized spacing [$\lambda$]', 'Interpreter', 'latex', 'FontSize', 22);
        ylabel('Fitness [Vm]', 'Interpreter', 'latex', 'FontSize', 22);
        xlim([intv(1) intv(end)]);
        set(get(fig, 'CurrentAxes'), 'FontSize', 16);

        print_plots(fig, ['i' num2str(sz) '_0_' opt '_fit']);
    end;
    close(fig);
end