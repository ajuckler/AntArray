%CONVERTFITNESS convert the fitness metric
%
% Convert the fitness metric from Vm to m² and inversely.
%
% [ ] = CONVERTFITNESS(infile, dist, inmode, evth)
%
% INPUT
%   inname: input fitness file name
%   dist:   distance from the array where to compute the fitness [mm]
%   inmode: input fitness mode (0 for volume, 1 for surface)
%   evth:   [optional] convert only max or everything (default=0)
%
% See also FITNESS ANTARRAY

% Copyright 2016, Antoine Juckler. All rights reserved.

function convertFitness(inname, dist, inmode, evth)
    if nargin < 4 || isempty(evth)
        evth = 0;
    else
        evth = evth > 0;
    end;
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
    fig = openfig(infile, 'new', 'invisible');
    kids = get(get(fig, 'CurrentAxes'), 'Children');
    
    if length(kids) ~= 2
        error 'Wrong format of fitness plot';
    end;
    
    % Find number of files
    % --------------------
    num_els = 1;
    if evth
        while exist([GA_fold '/1/arrangement_' num2str(num_els) '.dat'], 'file')
            num_els = num_els + 1;
        end;
        num_els = num_els - 1;
        newdata = zeros(3, length(get(kids(1), 'YData')));
    else
        newdata = zeros(1, length(get(kids(1), 'YData')));
    end;    
    clearvars kids
    close(fig);
    
    parallel_pool('start');
    if ~evth
        parfor i=1:length(newdata)
            currfile = [GA_fold '/' num2str(i) '/arrangement_1.dat'];
            if ~exist(currfile, 'file')
                error(['File not found at iteration ' num2str(i)]);
            end;
            ant = AntArray(currfile, [], [], [], cfg_path, 0);
            newdata(i) = fitness(ant, dist, ~inmode);
        end;
    else
        try
            dial = WaitDialog();
            dial.setMainString('Starting...');
            
            % Create save folder
            dirname = [datestr(now, 'yyyymmdd') '/' inname '_conv'];
            if ~exist(dirname, 'dir')
                mkdir(dirname);
            end;
            
            maxi = length(newdata);
            for i=1:maxi
                dial.setMainString(['Working on population ' ...
                        num2str(i) ' of ' num2str(maxi) '...']);
                vals = zeros(1, maxi);
                parfor j=1:num_els
                    currfile = [GA_fold '/' num2str(i) ...
                        '/arrangement_' num2str(j) '.dat'];
                    if ~exist(currfile, 'file')
                        error(['File ' num2str(j) ...
                            ' not found at iteration ' num2str(i)]);
                    end;
                    ant = AntArray(currfile, [], [], [], cfg_path, 0);
                    vals(j) = fitness(ant, dist, ~inmode);
                end;
                newdata(1, i) = vals(1);
                newdata(2, i) = sum(sum(vals))/numel(vals);
                [newdata(3, i), pos] = max(vals);
                
                % Save max
                subdir = [dirname '/' num2str(i) '/'];
                mkdir(subdir);
                arrgt1 = AntArray([GA_fold '/' num2str(i) ...
                    '/arrangement_1.dat'], [], [], [], cfg_path, 0).M;
                arrgt2 = AntArray([GA_fold '/' num2str(i) ...
                    '/arrangement_' num2str(pos) '.dat'], ...
                    [], [], [], cfg_path, 0).M;
                save([subdir 'arrangement_1.dat'], 'arrgt1', '-ASCII');
                save([subdir 'arrangement_' num2str(pos) '.dat'], ...
                    'arrgt2', '-ASCII');
                
                dial.terminate();
            end;
        catch ME
            switch ME.identifier
                case 'MyERR:Terminated'
                    warning 'Operation terminated by user';
                otherwise
                    rethrow(ME);
            end;
        end;
    end;
    
    fig = figure(1);
    plot(1:length(newdata), newdata(1, :), '-b', 'LineWidth', 2, ...
        'DisplayName', 'Converted max');
    hold on
    if evth
        plot(1:length(newdata), newdata(2, :), '-r', 'LineWidth', 2, ...
            'DisplayName', 'Converted mean');
        plot(1:length(newdata), newdata(3, :), '-g', 'LineWidth', 2, ...
            'DisplayName', 'New max');
        L = legend('Location', 'southeast');
        set(L, 'Interpreter', 'latex', 'FontSize', 20);
    end;  

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
    
    if evth
        savname = ['fitness_' inname '_' num2str(~inmode) '_evth'];
    else
        savname = ['fitness_' inname '_' num2str(~inmode)];
    end;
    print_plots(gcf, savname);
    close all
    
    parallel_pool('stop');
end