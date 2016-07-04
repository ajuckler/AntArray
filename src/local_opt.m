%LOCAL_OPT Optimizes the array arrangement using local search
%
%   The behaviour of the algorithm is determined by 2 hard-coded
%   parameters:
%       dist_prob:  distortion probability
%       th_prob:    remaining probability of finding a better solution
%   The inputted arrangement is distorted before the algorithm is applied
%   to it, so that if the arrangement is the optimum we get it back as
%   function output.
%   The other probability is used to discard several neighbours when the
%   probability of finding a better solution in the remaining set of
%   neighbours is smaller than the given one.
%
%   [OPTIM_SOL, OPTIM_VAL] = LOCAL_OPT(START_POP, DIST, QUANT, MODE)
%   INPUT:
%       START_POP:  initial arrangement [AntArray]
%       DIST:       distance from the array where the fitness function will
%                   be evaluated [mm]
%       QUANT:      (optional) is the array arrangement quantized?
%                   [boolean, default=0]
%       MODE:       (optional) mode used for fitness computations [boolean,
%                   default=0]
%   OUTPUT:
%       OPTIM_SOL:  the optimal array arrangement [AntArray]
%       OPTIM_VAL:  fitness of the optimal arrangement [double]
%
%   [...] = LOCAL_OPT(PATH) resumes the algorithm
%   INPUT:
%       PATH:   path to the folder containing the previous results
%
%   See also GA_2D FITNESS CHROM2MAT MAT2CHROM ANTARRAY

%   Copyright 2016, Antoine Juckler. All rights reserved.

function [optim_sol, optim_val] = local_opt(start_pop, dist, quanti, mode)

    if nargin > 1
        if isempty(quanti)
            quanti = 0;
        else
            quanti = quanti > 0;
        end;
        if isempty(mode)
            mode = 0;
        else
            mode = mode > 0;
        end;

        if ~isempty(dist) && (~isa(dist, 'numeric') || dist <= 0)
            error('Invalid DIST argument');
        end;
    end;
    
    if ~isa(start_pop, 'AntArray') && ~isa(start_pop, 'char')
        error('Invalid initial solution');
    end;
    
    dist_prob = 0.05;   % Distortion probability
    th_prob = 0.03;
    
    if ispref('loc_opt', 'save_folder')
        rmpref('loc_opt');
    end;
    
    parallel_pool('start');

    try
        dial = WaitDialog();
        dial.setMainString('Initialisation...');
        if ~isa(start_pop, 'char')
            chrom = mat2chrom(start_pop.M, quanti);
            % Distort arrangement
            pert = rand(1, length(chrom));
            pert = (pert <= dist_prob);
            
            chrom(pert == 1) = abs(chrom(pert == 1) - 1);
            start_pop = start_pop.rstArray(chrom2mat(chrom, quanti));

            % Save initial arrangement
            [opt_name opt_fold] = sav_loc_state(start_pop);
            start_pop = start_pop.setName(opt_name);
            start_pop.saveCfg();
            opt_name = [opt_fold 'cfg/' opt_name];

            cfg = [dist_prob th_prob dist quanti mode];
            sav_loc_state('cfg', cfg);

            best = start_pop;
            [best_val alt_val] = fitness(best, dist, mode);

            tmp_val = 0;
            visited_list = {};
            progress_data = {};
            progress_data{1} = [best_val alt_val];
            dial.setSubString('Done');
            dial.terminate();
            iter = 1;
        else
            dial.setSubString('Loading old data');
            if ~strcmp(start_pop(end), '/')
                start_pop = [start_pop '/'];
            end;
            opt_name = [start_pop(1:8) '/cfg/' start_pop(end-6:end-1)];
            fname = opt_name(end-5:end);
            opt_fold = start_pop(1:9);

            cfgname = [start_pop 'cfg.loc_opt'];
            if ~exist(cfgname, 'file')
                error 'File not found';
            else
                setpref('local_opt', 'save_folder', start_pop);
            end;
            cfg = dlmread(cfgname);
            dist_prob = cfg(1);
            th_prob = cfg(2);
            dist = cfg(end-2);
            quanti = cfg(end-1);
            mode = cfg(end);

            sav_loc_state();

            % Load savdata
            % ------------
            if ~exist([start_pop 'savdata.mat'], 'file')
                error 'Savdata not found';
            end;

            savdata = load([start_pop 'savdata.mat']);
            fl = fieldnames(savdata);
            if length(fl) == 1
                savdata = savdata.(fl{1});
            end;
            subiter = savdata.subiter;
            best = savdata.best;
            pos = savdata.pos;
            best_val = savdata.best_val;
            tmp_val = savdata.tmp_val;
            alt_val = savdata.alt_val;
            tmp_visited = savdata.tmp_visited;

            % Find last iteration and load chrom_ls
            % -------------------------------------
            iter = 1;
            visited_list = {};
            while exist([start_pop num2str(iter)], 'dir')
                if ~exist([start_pop num2str(iter) '/chroms.dat'])
                    error('MyERR:FileNotFound', ['Chromosome list not ' ...
                        'found at iteration ' num2str(iter)]);
                else
                    tmp = dlmread([start_pop num2str(iter) '/chroms.dat']);
                    tmp_cell = cell(1, size(tmp, 1));
                    for i=1:size(tmp, 1)
                        tmp_cell{i} = tmp(i,:);
                    end;
                    bpos = size(visited_list, 1);
                    visited_list(bpos+1:bpos+length(tmp_cell)) = tmp_cell;
                end;
                iter = iter + 1;
            end;

            % Load progress_data
            % ------------------
            progname = [opt_fold 'dat/fitness_loc_' fname '.dat'];
            if ~exist(progname, 'file')
                error('MyERR:FileNotFound', 'Progress file not found');
            else
                progress_data = dlmread(progname);
            end;
            tmp = cell(1, size(progress_data, 1));
            for i=1:length(tmp)
                tmp{i} = progress_data(i, :);
            end;
            progress_data = tmp;

            % Load best chrom
            % ---------------
            if iter == 1
                filename = [start_pop 'arrgt.dat'];
            else
                filename = [start_pop num2str(iter-1) '/arrgt.dat'];
            end;
            if ~exist(filename, 'file')
                error('MyERR:FileNotFound', 'Last arrangement not found');
            else
                best_chrom = dlmread(filename);
                best_chrom = mat2chrom(best_chrom, quanti);
            end;
            
            dial.terminate();
            dial.setMainString('Resuming last iteration');

            % Resume last iteration
            % ---------------------
            best_val = progress_data{end};
            best_val = best_val(1);
            
            maxiter = ceil((1-th_prob)*length(pos));
            nworkers = parallel_pool('status');
            params = struct();
            params.quanti = quanti;
            params.dist = dist;
            params.mode = mode;
            params.opt_name = opt_name;

            pariter = floor((maxiter-subiter)/2/nworkers);
            stopi = pariter;
            for i=1:stopi
                dial.setSubString(['Subiteration ' num2str(i) ' on ' ...
                    num2str(stopi)]);
                dial.terminate();

                if i <= pariter
                    beginit = subiter + (i-1)*2*nworkers+1;
                    endit = subiter + i*2*nworkers;
                elseif i == pariter+1
                    beginit = subiter + pariter*2*nworkers+1;
                    endit = maxiter;
                elseif i == pariter+2
                    beginit = maxiter;
                    endit = length(pos);
                else
                    error('MyERR:Invalid', 'Unexpected iteration value');
                end;

                structout = insideOpt(beginit, endit, ...
                    best_chrom, pos, visited_list, params);

                [mval, mpos] = max(structout.tmp_vect);
                if mval > tmp_val
                    tmp_val = mval;
                    best = structout.tmp_pop{mpos};
                    alt_val = structout.alt_vect(mpos);
                end;
                
                for ii=1:length(structout.tmp_visited_in)
                    val = structout.tmp_visited_in{ii};
                    if ~isempty(val)
                        tmp_visited{length(tmp_visited)+1} = val;
                    end;
                end;

                if i == pariter && pariter*2*nworkers+subiter ~= maxiter
                    stopi = stopi+1;
                elseif i >= pariter && tmp_val == best_val
                    stopi = stopi+2;
                    i = stopi+1;
                end;
                
                dial.terminate();
            end
            % while i <= maxiter
            %     dial.setSubString([num2str(i) ' on ' num2str(maxiter)]);
            %     dial.terminate();
            %     tmp_chrom = best_chrom;
            %     tmp_chrom(pos(i)) = ~tmp_chrom(pos(i));

            %     skip = knownNeighbour(visited_list, tmp_chrom);

            %     if ~skip
            %         tmp_visited{length(tmp_visited)+1} = tmp_chrom;
            %         tmp = AntArray(chrom2mat(tmp_chrom, quanti), ...
            %             [], [], [], opt_name, 0);

            %         [loc_val tmp_alt_val] = fitness(tmp, dist, mode);

            %          if loc_val > tmp_val
            %             tmp_val = loc_val;
            %             best = tmp;
            %             alt_val = tmp_alt_val;
            %         end;
            %     end

            %     i = i + 1;
            %     if i == maxiter && tmp_val == best_val
            %         maxiter = length(pos);
            %     end;
            %     dial.terminate();
            % end;
            
            delete([start_pop 'savdata.mat']);
            
            fname = sav_loc_state(best, tmp_visited, iter);
            progress_data{length(progress_data)+1} = [tmp_val alt_val];

            startpos = length(visited_list);
            visited_list(startpos+1:startpos+length(tmp_visited)) = tmp_visited;
            iter = iter + 1;
            dial.setSubString('Iteration data saved');
            dial.terminate();
        end;
    
        % ------------------------------ %
        %         Start / resume         %
        % ------------------------------ %
        fname = opt_name(end-5:end);
        params = struct();
        params.quanti = quanti;
        params.mode = mode;
        params.dist = dist;
        params.opt_name = opt_name;
        while tmp_val ~= best_val
            dial.setMainString(['Iteration ' num2str(iter) '...']);
            dial.terminate();
            if iter == 1
                tmp_val = best_val;
            end;
            best_val = tmp_val;
            best_chrom = mat2chrom(best.M, quanti);
            pos = genPairs(genMask(best_chrom));

            tmp_visited = {};

            maxiter = ceil((1-th_prob)*length(pos));
            nworkers = parallel_pool('status');
            pariter = floor(maxiter/2/nworkers);
            stopi = pariter;
            for i=1:stopi
                dial.setSubString(['Subiteration ' num2str(i) ' of ' ...
                    num2str(stopi)]);
                dial.terminate();

                if i <= pariter
                    beginit = (i-1)*2*nworkers+1;
                    endit = i*2*nworkers;
                elseif i == pariter+1
                    beginit = pariter*2*nworkers+1;
                    endit = maxiter;
                elseif i == pariter+2
                    beginit = maxiter;
                    endit = length(pos);
                else
                    error('MyERR:Invalid', 'Unexpected iteration value');
                end;
                
                structout = insideOpt(beginit, endit, ...
                    best_chrom, pos, visited_list, params);
                
                [mval, mpos] = max(structout.tmp_vect);
                if mval > tmp_val
                    tmp_val = mval;
                    best = structout.tmp_pop{mpos};
                    alt_val = structout.alt_vect(mpos);
                end;
                
                for ii=1:length(structout.tmp_visited_in)
                    val = structout.tmp_visited_in{ii};
                    if ~isempty(val)
                        tmp_visited{length(tmp_visited)+1} = val;
                    end;
                end;

                if i == pariter && pariter*2*nworkers+subiter ~= maxiter
                    stopi = stopi+1;
                elseif i >= pariter && tmp_val == best_val
                    stopi = stopi+2;
                    i = stopi+1;
                end;
                
                dial.terminate();
            end
            
            % % Run for rest
            % dial.setSubString(['Subiteration ' num2str(pariter+1) ' of ' ...
            %         num2str(pariter+1)]);
            % structout = insideOpt(pariter*2*nworkers+1, maxiter, ...
            %     best_chrom, pos, visited_list, params);
            % [mval, mpos] = max(structout.tmp_vect);
            % if mval > tmp_val
            %     tmp_val = mval;
            %     best = structout.tmp_pop{mpos};
            %     alt_val = structout.alt_vect(mpos);
            % end;
            % for ii=1:length(structout.tmp_visited_in)
            %     val = structout.tmp_visited_in{ii};
            %     if ~isempty(val)
            %         tmp_visited{length(tmp_visited)+1} = val;
            %     end;
            % end;
            % dial.terminate();
            
            % % If no max found, run more
            % if tmp_val == best_val
            %     dial.setSubString(['Subiteration ' num2str(pariter+2) ' of ' ...
            %         num2str(pariter+2)]);
            %     structout = insideOpt(maxiter+1, length(pos), ...
            %         best_chrom, pos, visited_list, params);
            %     [mval, mpos] = max(structout.tmp_vect);
            %     if mval > tmp_val
            %         tmp_val = mval;
            %         best = structout.tmp_pop{mpos};
            %         alt_val = structout.alt_vect(mpos);
            %     end;
            %     for ii=1:length(structout.tmp_visited_in)
            %         val = structout.tmp_visited_in{ii};
            %         if ~isempty(val)
            %             tmp_visited{length(tmp_visited)+1} = val;
            %         end;
            %     end;
            % end;

            fname = sav_loc_state(best, tmp_visited, iter);
            progress_data{length(progress_data)+1} = [tmp_val alt_val];

            startpos = length(visited_list);
            visited_list(startpos+1:startpos+length(tmp_visited)) = tmp_visited;
            iter = iter + 1;
            dial.terminate();
        end;
        delete(dial);
    catch ME
        switch ME.identifier
            case 'MyERR:Terminated'
                warning('MyWARN:Terminated', 'Operation terminated by user');
                savdata = struct();
                if i >= pariter+1
                    i == maxiter;
                else
                    savdata.subiter = i*2*nworkers;
                end;
                savdata.best = best;
                savdata.pos = pos;
                savdata.best_val = best_val;
                savdata.tmp_val = tmp_val;
                savdata.alt_val = alt_val;
                savdata.tmp_visited = tmp_visited;
                save([opt_fold opt_name(end-5:end) '/savdata.mat'], ...
                    'savdata', '-MAT', '-v7.3');

                % Convert to double matrix
                tmp = zeros(length(progress_data), length(progress_data{1}));
                for i=1:size(tmp,1)
                    tmp(i,:) = progress_data{i};
                end;
                progress_data_m = tmp;
                export_dat(progress_data_m, ...
                    ['fitness_loc_' opt_name(end-5:end)]);
            otherwise
                rethrow(ME);
        end;
    end;
    optim_sol = best;
    optim_sol = optim_sol.setName(fname);
    optim_sol = optim_sol.setComments(sprintf([...
        'Dist: ' num2str(dist/1000) '\nQuant: ' num2str(quanti) ...
        '\nMode: ' num2str(mode)]));
    optim_val = best_val;
    parallel_pool('stop');

    % Plot progress
    % -------------
    if ~exist('fname', 'var') || iter == 1
        warning('MyWARN:NoPlot', 'Progress plot no generated');
        return;
    end;
    % Convert to double matrix
    tmp = zeros(length(progress_data), length(progress_data{1}));
    for i=1:size(tmp,1)
        tmp(i,:) = progress_data{i};
    end;
    progress_data = tmp;

    fig = figure();
    hold on; box on;
    plot(1:size(progress_data,1), progress_data(:,1), '-b', ...
        'LineWidth', 2)

    xlim([1 length(progress_data)]);
    xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 22);
    if mode == 0
        label = 'Fitness [Vm]';
    else
        label = 'Fitness [$m^2$]';
    end;
    ylabel(label, 'Interpreter', 'latex', 'FontSize', 22);
    set(get(fig, 'CurrentAxes'), 'FontSize', 16);
    hold off;

    print_plots(fig, ['fitness_loc_' fname]);
    export_dat(progress_data, ['fitness_loc_' fname]);
    close all;
end

function [fname fdir] = sav_loc_state(best, chrom_ls, iter)
    config = 0;
    firstit = 0;

    if nargin < 2
        firstit = 1;
        iter = 0;
    elseif isa(best, 'char')
        if ~strcmp(best, 'cfg')
            error 'Invalid argument';
        else
            config = 1;
            iter = 0;
        end;
    else
        if min(size(chrom_ls)) ~= 1
            error 'Incorrect chromosome list';
        elseif iter < 1
            error 'Invalid iteration';
        end;
    end;

    persistent prefix_folder_name;

    if firstit
        if ispref('local_opt', 'save_folder')
            prefix_folder_name = getpref('local_opt', 'save_folder');
            rmpref('local_opt', 'save_folder');
        else
            prefix_folder_name = datestr(now, 'yyyymmdd/HHMMSS/');
        end;
    end;

    if iter ~= 0
        sav_loc = [prefix_folder_name num2str(iter)];
    else
        sav_loc = prefix_folder_name;
    end;

    if ~exist(sav_loc, 'dir')
        mkdir(sav_loc);
    end;

    if config
        save([sav_loc '/cfg.loc_opt'], 'chrom_ls', '-ASCII');
    elseif nargin > 0
        arrgt = AntArray.quantize(best.M, 1);
        save([sav_loc '/arrgt.dat'], 'arrgt', '-ASCII');

        if iter ~= 0
            chrom_mat = zeros(length(chrom_ls), length(chrom_ls{1}));
            for i=1:size(chrom_mat,1)
                chrom_mat(i,:) = chrom_ls{i};
            end;
            save([sav_loc '/chroms.dat'], 'chrom_mat', '-ASCII');
        end;
    end;

    fname = prefix_folder_name(end-6:end-1);
    fdir = prefix_folder_name(1:end-7);    
end

function skip = knownNeighbour(visited_list, tmp_chrom)
    skip = 0;
    for i=1:length(visited_list)
        if sum(visited_list{i} == tmp_chrom) == length(tmp_chrom)
            skip = 1;
            break;
        end;
    end;
end

function mask = genMask(chrom)
    dim = sqrt(length(chrom));
    mat = reshape(chrom, dim, dim);
    mask = mat;

    for i=1:dim
        for j=1:dim
            if mat(i,j)
                if i > 1
                    mask(i-1,j) = 1;
                end;
                if i < dim
                    mask(i+1,j) = 1;
                end;
                if j > 1
                    mask(i,j-1) = 1;
                end;
                if j < dim
                    mask(i,j+1) = 1;
                end;
            end;
        end;
    end;

    mask = reshape(mask, 1, numel(mask));
end

function pairs = genPairs(chrom)
    pairs = 1:length(chrom);

    pairs = pairs(chrom ~= 0);
    pairs = pairs(randperm(length(pairs)));
end

function structout = insideOpt(beginit, endit, best_chrom, pos, visited, params)
    quanti = params.quanti;
    dist = params.dist;
    mode = params.mode;
    opt_name = params.opt_name;
    
    rangeiter = endit-beginit+1;
    tmp_visited_in = cell(1, rangeiter);
    tmp_pop = cell(1, rangeiter);
    tmp_vect = zeros(1, rangeiter);
    alt_vect = zeros(1, rangeiter);
    tmp_pos = pos(beginit:endit);
    
    parfor ii=1:rangeiter
        currpos = tmp_pos(ii);
        tmp_chrom = best_chrom;
        tmp_chrom(currpos) = ~tmp_chrom(currpos);

        skip = knownNeighbour(visited, tmp_chrom);

        if ~skip
            tmp_visited_in{ii} = tmp_chrom;
            tmp_pop{ii} = AntArray(chrom2mat(tmp_chrom, quanti), ...
                [], [], [], opt_name, 0);

            [loc_val tmp_alt_val] = fitness(tmp_pop{ii}, dist, mode);
            
            tmp_vect(ii) = loc_val;
            alt_vect(ii) = tmp_alt_val;
        end
    end

    structout = struct();
    structout.tmp_visited_in = tmp_visited_in;
    structout.tmp_pop = tmp_pop;
    structout.tmp_vect = tmp_vect;
    structout.alt_vect = alt_vect;
end
