function genColor(name, opt)
    if nargin < 2 || isempty(opt)
        opt = 'XY';
    end;
    if strcmpi(opt, 'XY') && ~exist(['fig/pattern_XY_' name '.fig'], 'file')
        error 'File not found';
    elseif strcmpi(opt, 'YZ') && ~exist(['fig/pattern_' name '_2.fig'], 'file')
        error 'File not found';
    end;
    if strcmpi(opt, 'XY') || strcmpi(opt, 'theta')
        if strcmpi(opt, 'XY')
            fig = openfig(['fig/pattern_XY_' name]);
        else
            fig = openfig(['fig/pattern_theta_0_' name]);
        end;
        ax = get(fig, 'CurrentAxes');
        kids = get(ax, 'Children');
        xdata = get(kids, 'XData');
        ydata = get(kids, 'YData');
        xdim = 2 * xdata(1,end-1);
        ydim = ydata(end-1,1);
        xlims = get(ax, 'XLim');
        ylims = get(ax, 'YLim');
        set(ax, 'YTick', [ylims(1) ylims(end)/2 ylims(end)], ...
            'YTickLabel', [0 round(ylims(end))/2 round(ylims(end))]);
        set(ax, 'XTick', [xlims(1) (xlims(end)+xlims(1))/2 xlims(end)], ...
            'XTickLabel', [xlims(1) 0 -xlims(1)]);
        set(get(ax, 'Title'), 'FontSize', 24);
        set(get(ax, 'YLabel'), 'FontSize', 22);
        set(get(ax, 'XLabel'), 'FontSize', 22);
        set(ax, 'FontSize', 16);
        
        if strcmpi(opt, 'theta')
            xlabel(ax, '$\rm d_{\rm centre}$ [m]', 'FontSize', 22);
        end;

        set(ax, 'PlotBoxAspectRatio', [xdim ydim 1]);
        if strcmpi(opt, 'XY')
            print_plots(fig, ['pattern_XY_' name]);
        else
            print_plots(fig, ['pattern_theta_0_' name]);
        end;
        close(fig);
    elseif strcmpi(opt, 'YZ')
        fig = openfig(['fig/pattern_' name '_2.fig']);
        ax = get(fig, 'CurrentAxes');
        set(get(ax, 'YLabel'), 'FontSize', 22);
        set(get(ax, 'XLabel'), 'FontSize', 22);
        set(ax, 'FontSize', 16);

        print_plots(fig, ['pattern_' name '_2']);
        close(fig);
    else
        error 'Wrong option';
    end;
end