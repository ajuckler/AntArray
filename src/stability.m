function stability(ant, prob, dist, mode)
    %STABILITY Controls the stability of an antenna array arrangement
    %
    % This function permits to evaluate the stability of an antenna array
    % arrangement. The array will be subject to perturbations that will 
    % turn on/off elements with a certain probability. The fitness of the
    % perturbed array is then evaluated and compared to the one of the
    % original array.
    %
    % [] = STABILITY(ANT, PROB, DIST, MODE)
    % INPUT:
    %   ANT:    ANTARRAY object, its stability will be evaluated
    %   PROB:   probability of the occurence of a perturbation in the 
    %           array arrangement
    %   DIST:   distance from the array plane at which the fitness should
    %           be evaluated [mm]
    %   MODE:   (optional) if set to 1, the fitness will be computed from
    %           the surface of the cut plane through the beam; if set to 0
    %           it will be computed from the volume of the beam
    %           [default = 0]
    %
    % See also FITNESS ANTARRAY

    % Copyright 2015-2016, Antoine JUCKLER. All rights reserved.
    
    if prob >= 1
        error 'PROB should be < 1';
    end;
    if nargin < 4
        mode = 0;
    else
        mode = (mode > 0);
    end;

    max_iter = 100;

    % Start parallel pool
    parallel_pool('start');

    plotdata = zeros(1, max_iter+1);
    plotdata(end) = fitness(ant, dist, mode);

    mat = ant.M;

    parfor i = 1:max_iter
        % Generate noise matrix
        perturbation = rand(size(mat,1));
        perturbation = (perturbation <= prob);

        tmp = mat;
        tmp(perturbation == 1) = abs(tmp(perturbation == 1) - 1);

        plotdata(i) = fitness(AntArray(tmp), dist, mode);
    end;

    plotdata = round(plotdata.*10^4)./10^4;

    % Stop parallel pool
    parallel_pool('stop');

    % Plot
    % ----
    avg = mean(plotdata);
    std_dev = std(plotdata);
    med = median(plotdata);
    disp(['Mean: ' num2str(avg)]);
    disp(['Standard deviation: ' num2str(std_dev)]);
    disp(['Median: ' num2str(med)]);

    fig1 = figure(1);
    axes('Parent', fig1, 'Position', [0.13 0.11 0.65 0.8150], 'Box', 'on');
    hold on;

    if verLessThan('matlab','8.2')
        % http://stats.stackexchange.com/questions/798/calculating-optimal-number-of-bins-in-a-histogram-for-n-where-n-ranges-from-30
        q1 = quantile(plotdata, 0.25);
        q3 = quantile(plotdata, 0.75);
        if q1 == q3
            warning('MyWARN:bins', 'First ant third quartile were equal');
            i = 1;
            while q1 == q3 && i < 25
                q1 = quantile(plotdata, 0.25 - 0.1*i);
                q3 = quantile(plotdata, 0.75 + 0.1*i);
                i = i+1;
            end;
        end;
        h = 2*length(plotdata)^(-1/3)*(q3-q1);
        nbins = round((max(plotdata)-min(plotdata))/h);
        [counts, centers] = hist(plotdata, nbins);
        bar(centers, counts, 'BarWidth', 1, 'FaceColor', [0.4,0.67,0.84]);
        sp = centers(end) - centers(end-1);
        xmin = centers(1) - sp/2;
        xmin = xmin-0.05 + mod(xmin,0.05);
        xmax = centers(end) + sp/2;
        xmax = xmax+0.05 - mod(xmax,0.05);
        xlim([xmin xmax]);    
    else
        histogram(plotdata);
    end;

    fig_dim = axis;

    plot([plotdata(end) plotdata(end)], [fig_dim(3) fig_dim(4)], ...
        '--r', 'Linewidth', 1);
    text(1.02, .5, sprintf(['Mean: ' num2str(avg, 4) ...
        '\nStd: ' num2str(std_dev, 4) ...
        '\nMedian: ' num2str(med, 4) ...
        '\nProb: ' num2str(prob)]), 'Units', 'normalized', ...
         'Interpreter', 'latex', 'FontSize', 18);

    title('\textbf{Stability}', 'Interpreter', 'latex', 'FontSize', 24);
    ylabel('Count', 'Interpreter', 'latex', 'FontSize', 22');
    if mode == 0
        fit_unit = 'Vm';
    else
        fit_unit = '$m^2$';
    end;
    xlabel(['Fitness [' fit_unit ']'], 'Interpreter', 'latex', 'FontSize', 22);
    set(gca, 'FontSize', 16);
    hold off;

    savname = ['stability' ant.name '_' num2str(round(prob*100))];
    print_plots(gcf, savname);
    export_dat(plotdata, savname);

    close all;

end
