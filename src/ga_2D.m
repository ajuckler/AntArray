%GA_2D Optimize the array arrangement using a genetic algorithm
%
%   The algorithm parameters are hard-coded as follow:
%   chromosome size:    
%       if QUANT(see further):  1/4 of matrix size in START_POP or of
%                               AntArray default
%       if ~QUANT:              1/2 of matrix size in START_POP or of 
%                               AntArray default
%   population size:            50
%   tournament participants:    2
%   chromosomes passed through: 2
%   mutation probability:       0.001
%   maximum iteration:          50
%
%   The pattern for 1-'point' crossover is generated with the
%   GENCROSSOVERPATTERN function.
%   The fitness of the individuals is evaluated with the FITNESS function.
%
%
%   [OPTIM_SOL, OPTIM_VAL] = GA_2D(DIST, START_POP, QUANT, MODE)
%   INPUT:
%       DIST:       distance from the array where the fitness function will
%                   be evaluated [mm]
%       START_POP:  (optional) individuals to be included in the population
%                   [cell array of AntArray elements] OR array size
%       QUANT:      (optional) is the array arrangement quantized?
%                   [boolean]
%       MODE:       (optional) mode used for fitness computations
%
%   OUTPUT:
%       OPTIM_SOL:  the optimal array arrangement [AntArray]
%       OPTIM_VAL:  fitness of the optimal arrangement [double]
%
%   [...] = GA_2D(PATH) resumes the algorithm
%   INPUT:
%       PATH:   path to the folder containing the previous results
%
%   See also FITNESS GENCROSSOVERPATTERN CHROM2MAT MAT2CHROM

%   Copyright 2016, Antoine Juckler. All rights reserved.

function [optim_sol, optim_val] = ga_2D(dist, start_pop, quant, mode)

if nargin < 4
    mode = 0;
else
    mode = (mode > 0);
end;
if nargin < 3
    quant = 1;
else
    quant = (quant > 0);
end;
if nargin < 2
    start_pop = [];
end;

% Parameters
% ----------
chrom_sz = [];      % chromosome size
                    %   empty = 1/4 of start_pop or of AntArray default

pop_sz = 50;        % population size
trn_sz = 2;         % tournament selection size
fit_sz = trn_sz;    % number of elements passed through
                    %   must be a multiple of trn_sz

mut_prob_df = 0.001;
mut_prob = mut_prob_df; % mutation probability

max_iter = 50;      % max number of iterations

if mod(pop_sz, trn_sz) ~= 0
    pop_sz = pop_sz + mod(pop_sz, trn_sz);
end;

cfg = [pop_sz, trn_sz, fit_sz, mut_prob, max_iter];

% Start parallel pool
persistent poolobj

if verLessThan('matlab','8.2')
    if ~matlabpool('size')
        matlabpool open
    end;
else
    if isempty(gcp('nocreate'))
        poolobj = parpool;
    end;
end;

% Clear preferences
if ispref('ga_2D', 'save_folder')
    rmpref('ga_2D');
end;

try
    % Init wait window
    dial = WaitDialog();
    dial.setMainString('Starting...');
    
    if nargin > 1 || isa(dist, 'numeric')
        dial.setSubString('Parsing initial population');
        
        cfg = [cfg mode quant dist];
        save_state(cfg);
        
        % Generate population
        % -------------------
        pop = cell(1, pop_sz);
        if ~isempty(start_pop) % if some elements already given, check
            if iscell(start_pop)
                for i=1:length(start_pop)
                    if ~isa(start_pop{i}, 'AntArray')
                        error('START_POP is not of type AntArray');
                    elseif quant
                        start_pop{i} = AntArray(...
                                       AntArray.quantize(...
                                       round(abs(start_pop{i}.M)), 2, 1));
                    else
                        start_pop{i} = AntArray(...
                                       round(abs(start_pop{i}.M)));
                    end;
                end;
            elseif isa(start_pop, 'numeric')
                if mod(start_pop, 2)
                    error('START_POP should be an even number');
                else
                    chrom_sz = start_pop;
                end;
            else
                if ~isa(start_pop, 'AntArray')
                    error('START_POP is not of type AntArray or numeric');
                elseif quant
                    tmp = AntArray(...
                                   AntArray.quantize(...
                                   round(abs(start_pop.M)), 2, 1));
                    start_pop = cell(1,1);
                    start_pop{1} = tmp;
                else
                    tmp = round(abs(start_pop));
                    start_pop = cell(1,1);
                    start_pop{1} = tmp;
                end;
            end;

            if ~isa(start_pop, 'numeric')
                chrom_sz = size(start_pop{1}.M, 1);
            else
                start_pop = [];
            end;
        end;

        dial.setSubString('Generating population')

        % Continue initial population generation
        if isempty(chrom_sz)
            temp = AntArray();
            chrom_sz = size(temp.M, 1);
            clearvars temp;
        end;

        if mod(chrom_sz, 4) == 0 && quant
            chrom_sz = chrom_sz/4;
        elseif mod(chrom_sz, 2) == 0 && ~quant
            chrom_sz = chrom_sz/2;
        end;

        end_index = min(length(start_pop), pop_sz);
        if end_index ~= 0
            pop(1:end_index) = start_pop(1:end_index);
        end;
        for i=end_index+1:pop_sz
            chrom = randi([0 1], 1, chrom_sz*chrom_sz);
            pop{i} = AntArray(chrom2mat(chrom, quant));
            dial.terminate();
        end;

        % Reshape for future parallel implementation
        v_dim = pop_sz/trn_sz;
        pop = reshape(pop, [trn_sz, v_dim]);
        eva = zeros(trn_sz, v_dim);

        dial.setMainString('Initial population generated');
        dial.terminate();

        % First evaluation
        % ----------------
        dial.setSubString('Computing fitnesses');
        parfor i=1:v_dim
            pop_ln = pop(:, i);
            eva_ln = eva(:, i);
            for j=1:trn_sz
                eva_ln(j) = fitness(pop_ln{j}, dist, mode);
            end;
            eva(:, i) = eva_ln;
        end;

        dial.terminate();
        dial.setMainString('Initial population evaluated');

        iter = 1;
        progress_data = zeros(max_iter+1, 2);
        save_state(pop, eva, iter);
        dial.setSubString('Initial data saved');
        dial.terminate();
        
    elseif nargin == 1 && isa(dist, 'char')
        dial.setSubString('Parsing old data');
        
        % Find folder
        % -----------
        if ~strcmp(dist(end), '/')
            dist = [dist '/'];
        end;
        if ~exist(dist, 'dir') || ~exist([dist '1'], 'dir')
            error('MyERR:InvalidFolder', ...
                ['Specified folder does not exist or does not ' ...
                'contain genetic algorithm data']);
        else
            setpref('ga_2D', 'save_folder', dist);
        end;
        
        % Parsing cfg data
        % ----------------
        if ~exist([dist 'config.ga_2D'], 'file')
            error 'Configuration file not found in the given folder';
        end;
        cfg = dlmread([dist 'config.ga_2D']);
        dir = dist;
        dist = cfg(end);
        
        if mod(cfg(1), cfg(2)) ~= 0
            error 'Erroneous population or turnament size in config file';
        elseif mod(cfg(3), cfg(2)) ~= 0
            error 'Erroneous turnament or fitness size in config file';
        elseif cfg(4) >= 1
            error 'Too high mutation probability in config file';
        elseif cfg(5) <= 1
            error 'Invalid maximal iteration in config file';
        end;
        quant = cfg(end-1);
        mode = cfg(end-2);
        max_iter = cfg(5);
        
        % Parsing progress_data
        % ---------------------
        progress_data = zeros(max_iter+1, 2);
        iter = 1;
        while exist([dir num2str(iter)], 'dir') && iter < max_iter
            if ~exist([dir num2str(iter) '/fitness.dat'], 'file')
                break;
            end;
            
            fit = dlmread([dir num2str(iter) '/fitness.dat']);
            progress_data(iter, 1) = fit(1,2);
            progress_data(iter, 2) = sum(fit(:,2))/size(fit,1);
            pop_sz = size(fit,1);
            
            iter = iter + 1;
        end;
        
        if pop_sz ~= cfg(1)
            error 'Population size does not correspond with config data';
        else
            trn_sz = cfg(2);
            fit_sz = cfg(3);
            mut_prob = cfg(4);
        end;
        
        % Population
        % ----------
        iter = iter - 1;
        if iter < 1
            error('There is not enough simulation data to resume simulations');
        end;
        pop = cell(1, pop_sz);
        fold_name = [dir num2str(iter) '/arrangement_'];
        for i=1:pop_sz
            tmp = dlmread([fold_name num2str(i) '.dat']);
            pop{i} = AntArray(tmp);
            
            quant_cond = sum(sum(tmp(1:2:end,:)~= tmp(2:2:end,:))) == 0;
            quant_cond = quant_cond && sum(sum(tmp(:,1:2:end) ~= tmp(:,2:2:end))) == 0;
            quant = quant_cond && quant;
            
            chrom_sz = size(tmp,1);
        end;
        
        v_dim = pop_sz/trn_sz;
        pop = reshape(pop, [trn_sz, v_dim]);
        eva = reshape(fit(:,2), [trn_sz, v_dim]);
        
        if mod(chrom_sz, 4) == 0 && quant
            chrom_sz = chrom_sz/4;
        elseif mod(chrom_sz, 2) == 0 && ~quant
            chrom_sz = chrom_sz/2;
        end;
        
        dial.setMainString('Old data parsed');
            
    else
        error('MyERR:InvalidInputArg', 'Invalid input argument DIST');
    end;
    
    % Store max & mean for progress tracking
    progress_data(iter, :) = [max(eva(:)) sum(sum(eva))/numel(eva)];

    condition = 0;
    while iter <= max_iter
        dial.setMainString(['Working on population ' num2str(iter) '...']);

        % Pass best individuals through
        % -----------------------------
        for i=1:fit_sz
            [val, pos] = max(eva(i:end));
            temp_ind = pop{i};
            temp_val = eva(i);
            pop{i} = pop{pos+i-1};
            eva(i) = val;
            pop{pos+i-1} = temp_ind;
            eva(pos+i-1) = temp_val;
        end;

        pop_tmp = pop;  % Tmp variable needed for parfor-loop
        eva_tmp = eva;
        
        dial.terminate();

        parfor i=fit_sz/trn_sz+1:v_dim
            % Selection
            % ---------
            inds = pop(:, i);
            vals = eva(:, i);
            for j=1:trn_sz
                rand_index = trn_sz+j;
                while rand_index == (i-1)*trn_sz+j
                    rand_index = randi(pop_sz, 1);
                end;

                if vals(j) < eva(rand_index)
                    inds{j} = pop{rand_index};
                end;

                % Convert to chromosomes
                inds{j} = mat2chrom(inds{j}.M, quant);
            end;

            % Crossover
            % ---------
            cross_patrn = genCrossoverPattern(chrom_sz, quant);
            cross_patrn = reshape(cross_patrn, 1, numel(cross_patrn));
            chroms = cell(1, trn_sz);
            if trn_sz ~= 2
                for j=1:trn_sz
                    temp_chrom = inds{j};   % Parent 1
                    swap_part = randi([0 1],1);
                    k = j;
                    while k == j
                        k = randi(trn_sz, 1);  % Parent 2
                    end;
                    spouse = inds{k};
                    temp_chrom(cross_patrn == swap_part) = ...
                        spouse(cross_patrn == swap_part);   % Child
                    chroms{j} = temp_chrom;
                end;
                clearvars swap_part spouse
            else
                par1 = inds{1};     % Parents
                par2 = inds{2};
                temp_chrom = par1;
                temp_chrom(cross_patrn == 1) = ...
                    par2(cross_patrn == 1);     % Child 1
                chroms{1} = temp_chrom;
                temp_chrom = par2;
                temp_chrom(cross_patrn == 1) = ...
                    par1(cross_patrn == 1);     % Child 2
                chroms{2} = temp_chrom;
            end;
            % clearvars temp_chrom cross_patrn;
            temp_chrom = [];
            cross_patrn = [];

            % Mutation
            % --------
            if rand <= mut_prob
                pos = randi(chrom_sz*chrom_sz, 1);  % Mutation position
                chrom_nb = randi([1 trn_sz], 1);    % Affected chromosome
                mut_chrom = chroms{chrom_nb};
                mut_chrom(pos) = ~mut_chrom(pos);
                chroms{chrom_nb} = mut_chrom;
%                 clearvars mut_chrom chrom_nb pos
                mut_chrom = [];
                chrom_nb = [];
                pos = [];
            end;

            % Save new pop & evaluate
            % -----------------------
            for j=1:trn_sz
                inds{j} = AntArray(chrom2mat(chroms{j}, quant));
                vals(j) = fitness(inds{j}, dist, mode);
            end;
            eva_tmp(:, i) = vals;
            pop_tmp(:, i) = inds;
                
        end;

        eva = eva_tmp;
        pop = pop_tmp;

        iter = iter+1;
        
        % Store max & mean for progress tracking
        progress_data(iter, :) = [max(eva(:)) sum(sum(eva))/numel(eva)];
        
        % Evaluate condition
        % ------------------
        if iter > round(.25*max_iter)
            last_data = progress_data(iter-round(.25*max_iter):iter, :);
            last_max = last_data(:, 1);
            last_mean = last_data(:, 2);
            
            if isempty(last_max(abs(last_max - last_max(end))>10^-6))
                mean_mean = sum(last_mean)/numel(last_mean);
                diff_mean = abs(last_mean - mean_mean);
                if isempty(diff_mean(diff_mean > .05*mean_mean))
                    condition = 1;
                end;
            end;
        end;
        
        % Force mutation if system is stable for too long
        if condition
            mut_prob = 0.25;
            condition = 0;
        else
            mut_prob = mut_prob_df;
        end;

        fname = save_state(pop, eva, iter);
    end;
    dial.terminate();

    [optim_val, pos] = max(eva(:));
    optim_sol = pop{pos};
    
    delete(dial);
    
catch ME
    switch ME.identifier
        case 'MyERR:Terminated'
            warning 'Operation terminated by user';
        otherwise
            rethrow(ME);
    end;
end;

% Stop parallel pool
if verLessThan('matlab','8.2') 
    if matlabpool('size')
        matlabpool close
    end;
else
    if ~isempty(gcp('nocreate'))
        delete(poolobj);
        poolobj = [];
    end;
end;

% Plot progress
% -------------
if ~exist('fname', 'var')
    warning 'Fitness plot not generated';
    return;
else
    optim_sol = optim_sol.setName(fname);
    optim_sol = optim_sol.setComments(sprintf([...
        'Dist: ' num2str(dist/1000) '\nQuant: ' num2str(quant) ...
        '\nMode: ' num2str(mode)]));
end;

figure(1);
plot(1:length(progress_data), progress_data(:,1), '-b', ...
    'LineWidth', 2, 'DisplayName', 'max');
hold on
plot(1:length(progress_data), progress_data(:,2), '-r', ...
    'LineWidth', 2, 'DisplayName', 'mean');

xlim([1 length(progress_data)]);

L = legend('Location', 'southeast');
set(L, 'Interpreter', 'latex', 'FontSize', 20);

xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Fitness', 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'FontSize', 16);
hold off;

print_plots(gcf, ['fitness_' fname]);
close all
    
end