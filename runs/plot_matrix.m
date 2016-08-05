function plot_matrix(mat, name)
    if nargin < 2 || isempty(name)
        name = datestr(now, 'HHMMSS');
    end;
    fig = figure('Visible', 'off');
    hold on; box off;
    imagesc(mat);
    colormap(fig, flipud(gray));
    caxis([0 1.75]);
    set(get(fig, 'CurrentAxes'), 'PlotBoxAspectRatio', [1 1 1]);
    xlim([.5 length(mat)+.5]);
    ylim([.5 length(mat)+.5]);
    plotview = axis;
    for i=1:length(mat)+1
        plot([i-.5 i-.5], [plotview(3) plotview(4)], 'k', 'LineWidth', .8);
        plot([plotview(1) plotview(2)], [i-.5 i-.5], 'k', 'LineWidth', .8);
    end;
    set(get(fig, 'CurrentAxes'), 'XTick', [], 'YTick', []);
    hold off;
    print_plots(fig, [name '_mat_BW']);
    close(fig);

    fig = figure('Visible', 'off');
    hold on; box off;
    textStrings = num2str(mat(:),'%u');
    textStrings = strtrim(cellstr(textStrings));
    [x, y] = meshgrid(1:length(mat), -1:-1:-length(mat));
%     [x, y] = meshgrid(1:length(mat), 1:length(mat));
    text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center', ...
        'FontSize', 20);
    xlim([.5 length(mat)+.5]);
    ylim([-length(mat)-.5 -.5]);
    set(get(fig, 'CurrentAxes'), 'PlotBoxAspectRatio', [1 1 1]);
    plotview = axis;
    for i=1:length(mat)+1
        plot([i-.5 i-.5], [plotview(3) plotview(4)], 'k', 'LineWidth', .8);
        plot([plotview(1) plotview(2)], [-i+.5 -i+.5], 'k', 'LineWidth', .8);
    end;
    set(get(fig, 'CurrentAxes'), 'XTick', [], 'YTick', []);
    hold off;
    print_plots(fig, [name '_mat_bin']);
    close(fig);
end