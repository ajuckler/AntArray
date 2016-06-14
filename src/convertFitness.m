%CONVERTFITNESS convert the fitness metric
%
% Convert the fitness metric from Vm to m² and inversely.
%
% [ ] = CONVERTFITNESS(infile, dist, inmode)
%
% INPUT
%   inname: input fitness file name
%   dist:   distance from the array where to compute the fitness [mm]
%   inmode: input fitness mode (0 for volume, 1 for surface)
%
% See also FITNESS ANTARRAY

% Copyright 2016, Antoine Juckler. All rights reserved.

function convertFitness(inname, dist, inmode)
    inmode = inmode > 0;
    
    % Check that files exist
    if length(inname) > 4 ...
            && (strcmp(inname(end-3:end), '.pdf') ...
            || strcmp(inname(end-3:end), '.fig'))
        inname = inname(1:end-3);
    end;
    i = 0;
    for i=0:10
        str = datestr(addtodate(now, -i, 'day'), 'yyyymmdd');
        if exist([str '/fig/fitness_' inname '.fig'], 'file')
            infile = [str '/fig/fitness_' inname '.fig'];
            break;
        elseif i == 10
            error 'Could not find specified file';
        end;
    end;

    for j=0:10
        str = datestr(addtodate(now, -i-j, 'day'), 'yyyymmdd');
        if exist([str '/' inname], 'dir')
            break;
        elseif j == 10
            error 'Could not find GA files';
        end;
    end;
    
    GA_fold = [str '/' inname];
    
    % Check for config file availability
    k = max(0, i-1);
    for j=k:i
        str = datestr(addtodate(now, -j, 'day'), 'yyyymmdd');
        if exist([str '/cfg/' inname '.cfg'], 'file')
            break;
        elseif j == i
            error 'Could not find config file';
        end;
    end;
    
    cfg_path = [str '/cfg/' inname '.cfg'];
    
    % Open files
    % ----------
    fig = openfig(infile);
    kids = get(get(fig, 'CurrentAxes'), 'Children');
    
    if length(kids) ~= 2
        error 'Wrong format of fitness plot';
    end;
    
    newdata = zeros(1, length(get(kids(1), 'YData')));
    clearvars kids
    
    parallel_pool('start');
    for i=1:length(newdata)
        currfile = [GA_fold '/' num2str(i) '/arrangement_1.dat'];
        if ~exist(currfile, 'file')
            error(['File not found at iteration ' num2str(i)]);
        end;
        ant = AntArray(currfile, [], [], [], cfg_path, 0);
        newdata(i) = fitness(ant, dist, ~inmode);
    end;
    
    fig = figure(1);
    plot(1:length(newdata), newdata, '-b', 'LineWidth', 2);
    hold on

    xlim([1 length(newdata)]);

    xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 22);
    if ~inmode == 0
        label = 'Fitness [Vm]';
    else
        label = 'Fitness [$m^2$]';
    end;
    ylabel(label, 'Interpreter', 'latex', 'FontSize', 22);
    set(get(fig, 'CurrentAxes'), 'FontSize', 16);
    hold off;

    print_plots(gcf, ['fitness_' inname '_' num2str(~inmode)]);
    close all
    
    parallel_pool('stop');
end