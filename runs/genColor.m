function genColor(name)
    if ~exist(['fig/pattern_XY_' name '.fig'], 'file')
        error 'File not found';
    end;
    fig = openfig(['fig/pattern_XY_' name]);
    ax = get(fig, 'CurrentAxes');
    kids = get(ax, 'Children');
    xdata = get(kids, 'XData');
    ydata = get(kids, 'YData');
    xdim = 2 * xdata(1,end-1);
    ydim = ydata(end-1,1);
    set(ax, 'Title', []);
    
    set(ax, 'PlotBoxAspectRatio', [xdim ydim 1]);
    print_plots(fig, ['pattern_XY_' name]);
    close(fig);
end