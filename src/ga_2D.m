%GA_2D Optimize the array arrangement using a genetic algorithm
%
%   The algorithm parameters are hard-coded as follow:
%   chromosome side length:    
%       if QUANT(see further):  1/4 of matrix size in START_POP or of
%                               AntArray default
%       if ~QUANT:              1/2 of matrix size in START_POP or of 
%                               AntArray default
%   population size:            70
%   tournament participants:    2
%   chromosomes passed through: 2
%   mutation probability:       0.001 bits
%   maximum iteration:          60 and more
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
%   See also LOCAL_OPT FITNESS CHROM2MAT MAT2CHROM

%   Copyright 2015-2016, Antoine JUCKLER. All rights reserved.

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

optim_sol = [];
optim_val = -1;

% Parameters
% ----------
chrom_sz = [];      % chromosome size
                    %   empty = 1/4 of start_pop or of AntArray default

pop_sz = 70;        % population size
trn_sz = 2;         % tournament selection size
fit_sz = trn_sz;    % number of elements passed through
                    %   must be a multiple of trn_sz

mut_prob_df = 0.001;
mut_prob = mut_prob_df; % mutation probability

max_iter = 60;      % max number of iterations
mut_iter = round(.2*max_iter); % max number of iterations before
                               %   change in mutation probability
off_ratio = .1;     % ratio of max_iter on which convergence will be checked
                                
freq = 60500;
elsp = 0.84;

if mod(pop_sz, trn_sz) ~= 0
    pop_sz = pop_sz + mod(pop_sz, trn_sz);
end;

if trn_sz < 2
    error 'Turnament size must be greater than 1';
end;

cfg = [pop_sz, trn_sz, fit_sz, mut_prob, max_iter];

% Start parallel pool
parallel_pool('start');

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
                                       round(abs(start_pop{i}.M)), 2, 1), ...
                                       freq, [], elsp, [], 0);
                    else
                        start_pop{i} = AntArray(...
                                       round(abs(start_pop{i}.M)), ...
                                       freq, [], elsp, [], 0);
                    end;
                end;
            elseif isa(start_pop, 'numeric')
                if mod(start_pop, 2)
                    error('START_POP should be an even number');
                elseif start_pop > 64
                    error(['Array dimensions too high: ' num2str(start_pop)]);
                else
                    chrom_sz = start_pop;
                end;
            else
                if ~isa(start_pop, 'AntArray')
                    error('START_POP is not of type AntArray or numeric');
                elseif quant
                    tmp = AntArray(...
                                   AntArray.quantize(...
                                   round(abs(start_pop.M)), 2, 1), ...
                                   freq, [], elsp, [], 0);
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
            temp = AntArray([], freq, [], elsp, [], 0);
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
            layout = genLayout(chrom_sz,1);
            chrom = reshape(layout{1}, 1, numel(layout{1}));
            % chrom = randi([0 1], 1, chrom_sz*chrom_sz);
            pop{i} = AntArray(...
                    chrom2mat(chrom, quant), freq, [], elsp, [], 0);
            dial.terminate();
        end;

        % Reshape for future parallel implementation
        v_dim = pop_sz/trn_sz;
        pop = reshape(pop, [trn_sz, v_dim]);
        eva = zeros(trn_sz, v_dim);
        eva_c = eva;

        dial.setMainString('Initial population generated');
        dial.terminate();

        % First evaluation
        % ----------------
        dial.setSubString('Computing fitnesses');
        parfor i=1:v_dim
            pop_ln = pop(:, i);
            eva_ln = eva(:, i);
            eva_c_ln = eva_c(:, i);
            for j=1:trn_sz
                [eva_ln(j), eva_c_ln(j)] = fitness(pop_ln{j}, dist, mode);
            end;
            eva(:, i) = eva_ln;
            eva_c(:, i) = eva_c_ln;
        end;

        dial.terminate();
        dial.setMainString('Initial population evaluated');

        iter = 1;
        off = 0;
        progress_data = zeros(max_iter+1, 2);
        progress_data_conv = progress_data;
        
        cfg = [cfg mode quant dist];
        ant = AntArray(pop{1,1}.M, freq, [], elsp, [], 0);
        save_state(cfg);
        ant.saveCfg();
        save_state(pop, eva, eva_c, 0);
        
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
        conv = 1;
        while exist([dir num2str(iter)], 'dir')
            if ~exist([dir num2str(iter) '/fitness.dat'], 'file')
                break;
            end;
            
            fit = dlmread([dir num2str(iter) '/fitness.dat']);
            progress_data(iter, 1) = fit(1,2);
            progress_data(iter, 2) = sum(fit(:,2))/size(fit,1);
            pop_sz = size(fit,1);
            
            if conv && exist([dir num2str(iter) '/fitness_conv.dat'], 'file')
                fit_c = dlmread([dir num2str(iter) '/fitness_conv.dat']);
            else
                conv = 0;
            end;
            
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
            pop{i} = AntArray(tmp, freq, [], elsp, [], 0);
            
            quant_cond = sum(sum(tmp(1:2:end,:)~= tmp(2:2:end,:))) == 0;
            quant_cond = quant_cond && sum(sum(tmp(:,1:2:end) ~= tmp(:,2:2:end))) == 0;
            quant = quant_cond && quant;
            
            chrom_sz = size(tmp,1);
        end;
        
        v_dim = pop_sz/trn_sz;
        pop = reshape(pop, [trn_sz, v_dim]);
        eva = reshape(fit(:,2), [trn_sz, v_dim]);
        if conv
            eva_c = reshape(fit_c(:,2), [trn_sz, v_dim]);
        else
            eva_c = [];
        end;
        
        if mod(chrom_sz, 4) == 0 && quant
            chrom_sz = chrom_sz/4;
        elseif mod(chrom_sz, 2) == 0 && ~quant
            chrom_sz = chrom_sz/2;
        end;
        
        if iter > max_iter
            off = iter-max_iter;
            mult = off_ratio*max_iter;
            off = mult*ceil(off/mult);
        else
            off = 0;
        end;
        
        dial.setMainString('Old data parsed');
            
    else
        error('MyERR:InvalidInputArg', 'Invalid input argument DIST');
    end;
    
    % Store max & mean for progress tracking
    progress_data(iter, :) = [max(eva(:)) sum(sum(eva))/numel(eva)];
    
    if isempty(eva_c)
        no_conv = 1;
        eva_c = eva;
    else
        no_conv = 0;
    end;

    condition = 0;
    while iter <= max_iter + off
        dial.setMainString(['Working on population ' num2str(iter)...
            ' of ' num2str(max_iter + off) '...']);

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
            tmp_val = eva_c(i);
            eva_c(i) = eva_c(pos+i-1);
            eva_c(pos+i-1) = tmp_val;
        end;
        if no_conv
            fname = save_state(pop, eva, iter);
        else
            fname = save_state(pop, eva, eva_c, iter);
        end;
        
        pop_tmp = pop;  % Tmp variable needed for parfor-loop
        eva_tmp = eva;
        eva_c_tmp = eva_c;
        
        dial.terminate();

        parfor i=fit_sz/trn_sz+1:v_dim
            % Selection
            % ---------
            inds = pop(:, i);
            vals = eva(:, i);
            vals_c = eva_c(:, i);
            
            for j=1:trn_sz
                rand_index = randi(pop_sz, 1);
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
                % clearvars swap_part spouse
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
            % % Pure random mutation positions
            % for j=1:trn_sz
            %     mask = rand(chrom_sz);
            %     mask = (mask <= mut_prob);
            %     mut_chrom = chroms{j};
            %     mut_chrom(mask == 1) = ...
            %         abs(mut_chrom(mask == 1) - 1);
            %     chroms{j} = mut_chrom;
            % end;
            
            % Mutation positions in clusters
            for j=1:trn_sz
                mask = genMask(chroms{j}, round(chrom_sz/16));
                mask2 = rand(1, numel(mask(mask == 1)));
                mask2 = (mask2 <= mut_prob);
                mask(mask == 1) = mask2;
                mut_chrom = chroms{j};
                mut_chrom(mask == 1) = ...
                    abs(mut_chrom(mask == 1) - 1);
                chroms{j} = mut_chrom;
            end;
            
            % % Mutation positions in the neighbourhood of the chromosome
            % for j=1:trn_sz
            %     mask = rand(chrom_sz);
            %     mask = (mask <= mut_prob);
            %     tot_els = numel(find(mask));
            %     if tot_els > 1
            %         nb_els = tot_els;
            %         mask = zeros(chrom_sz);
            %         while nb_els > 0
            %             tpmask = zeros(chrom_sz);
            %             % Gen start pos
            %             xpos = randi([1 chrom_sz]);
            %             ypos = randi([1 chrom_sz]);
            %
            %             % Gen height
            %             maxh = min(nb_els, chrom_sz-ypos+1);
            %             height = randi([1 maxh]);
            % 
            %             % Gen width
            %             maxw = min(floor(nb_els/height), chrom_sz-xpos+1);
            %             width = randi([1 maxw]);
            % 
            %             % Rotate
            %             tpmask(ypos:ypos+height-1, xpos:xpos+width-1) = 1;
            %             if rand < .5
            %                 tpmask = tpmask';
            %             end;
            % 
            %             % Adapt
            %             mask(tpmask == 1) = 1;
            %             nb_els = tot_els - numel(find(mask));
            %         end;
            %     end                  
            %     mut_chrom = chroms{j};
            %     mut_chrom(mask == 1) = ...
            %         abs(mut_chrom(mask == 1) - 1);
            %     chroms{j} = mut_chrom;
            % end;

            % Save new pop & evaluate
            % -----------------------
            for j=1:trn_sz
                inds{j} = AntArray(chrom2mat(chroms{j}, quant), ...
                                    freq, [], elsp, [], 0);
                [vals(j), vals_c(j)] = fitness(inds{j}, dist, mode);
            end;
            eva_tmp(:, i) = vals;
            pop_tmp(:, i) = inds;
            eva_c_tmp(:, i) = vals_c;
                
        end;

        eva = eva_tmp;
        pop = pop_tmp;
        eva_c = eva_c_tmp;

        iter = iter+1;
        
        % Store max & mean for progress tracking
        progress_data(iter, :) = [max(eva(:)) sum(sum(eva))/numel(eva)];
        
        % Evaluate condition
        % ------------------
        if iter > mut_iter
            last_data = progress_data(iter-mut_iter:iter, :);
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
            mut_prob = 0.1;
            condition = 0;
        else
            mut_prob = mut_prob_df;
        end;
        
        % Add iterations if no good convergence
        % -------------------------------------
        if iter > max_iter + off
            off_iter = round(off_ratio*max_iter);
            % Test enough identical samples
            last_data = progress_data(iter-off_iter:iter, :);
            last_max = last_data(:, 1);
            if ~isempty(last_max(abs(last_max - last_max(end))>10^-6))
                off = off + off_iter;
            end;
        end;
    end;
    if no_conv
        save_state(pop, eva, iter);
    else
        save_state(pop, eva, eva_c, iter);
    end;
    dial.terminate();
    
    delete(dial);
    
catch ME
    switch ME.identifier
        case 'MyERR:Terminated'
            warning('MyWARN:Terminated', 'Operation terminated by user');
        otherwise
            rethrow(ME);
    end;
end;

% Stop parallel pool
parallel_pool('stop');

% Plot progress
% -------------
if ~exist('fname', 'var')
    warning('MyWARN:NoPlot', 'Fitness plot not generated');
    return;
elseif iter > max_iter
    [optim_val, pos] = max(eva(:));
    optim_sol = pop{pos};
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
if mode == 0
    label = 'Fitness [Vm]';
else
    label = 'Fitness [$m^2$]';
end;
ylabel(label, 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'FontSize', 16);
hold off;

print_plots(gcf, ['fitness_' fname]);
close all
    
end

%% Function to generate the mask for crossover
function patrn = genCrossoverPattern(side_lg, quant)
    %GENCROSSOVERPATTERN generates mask for chromosome 1-point crossover
    %
    % This function generates a mask that will be used for chromosome
    % crossover. The mask correspond to 1-point crossover in 2D: when the
    % full matrix is recontructed from the chromosome, the mask will 
    % delimit 2 contiguous regions.
    %
    % PATTERN = GENCROSSOVERPATTERN(SIDE_LEN, QUANT)
    % INPUT:
    %     SIDE_LEN:   side length of the 2D chromosome
    %     QUANT:      (optional) is the chromosome quantized? [default = 1]
    % OUTPUT:
    %     PATTERN:    matrix containing the mask for 2D 1-point crossover
    
    if nargin < 2
        quant = 1;
    end;

    if mod(side_lg, 2) ~= 0
        error 'side_lg is not an even number';
    elseif side_lg < 5
        error 'side_lg must be greater than 4';
    end;

    % Theoretical size of matrix (effective elements generated)
    if quant
        step = 2;
    else
        step = 1;
    end;

    th_side = side_lg/step;
    th_diff = ceil(th_side/8);
    mat = zeros(th_side);
    gen_max_off = th_side;
    
    maxit = 20;

    % Generate theoretical pattern
    prev_off = randi([2 th_side]);
    prev_lgt = th_side - prev_off + 1;
    mat(end, prev_off:end) = 1;
    count_l = 1;
    count_r = 1;
    it1 = 0;
    it2 = 0;
    
    maxi = randi([th_diff th_side-2]);
    i = 1;
    while i <= maxi
        min_off = prev_off - th_diff;
        min_off = max(min_off, 2);
        
        max_off = prev_off + th_diff;
        if prev_lgt <= th_diff
            max_off = min(max_off, prev_off + prev_lgt - 1);
        end;
        max_off = min(max_off, gen_max_off);
        
        off = randi([min_off, max_off]);
        cond = (count_l == th_diff && off == prev_off) ...
            || (prev_lgt == 1 && count_r == th_diff && off == prev_off + prev_lgt) ...
            || (off == gen_max_off ...
                && (count_r >= th_diff-1 || count_l >= th_diff-1)) ...
            || (off == prev_off && prev_lgt == 2 ...
                && count_r >= th_diff-1 && count_l >= th_diff-1 ...
                && prev_off+prev_lgt == gen_max_off+1);
        while cond && it1 < maxit
            off = randi([min_off, max_off]);
            cond = (count_l == th_diff && off == prev_off) ...
            || (prev_lgt == 1 && count_r == th_diff && off == prev_off + prev_lgt) ...
            || (off == gen_max_off ...
                && (count_r >= th_diff-1 || count_l >= th_diff-1)) ...
            || (off == prev_off && prev_lgt == 2 ...
                && count_r >= th_diff-1 && count_l >= th_diff-1 ...
                && prev_off+prev_lgt == gen_max_off+1);
            it1 = it1 + 1;
        end;
        
        max_lgt = prev_off + prev_lgt - off + th_diff;
        max_lgt = min(gen_max_off - off + 1, max_lgt);
        
        min_lgt = prev_off + prev_lgt - th_diff - off;
        if off < prev_off
            min_lgt = max(prev_off - off + 1, min_lgt);
        end;
        min_lgt = max(1, min_lgt);
        
        lgt = randi([min_lgt, max_lgt]);
        cond = (count_r == th_diff && off + lgt == prev_off + prev_lgt) ...
                || (count_l >= th_diff-1 && lgt == 1 && off == prev_off);
        while gen_max_off ~= th_side && cond && it2 < maxit
            lgt = randi([min_lgt, max_lgt]);
            cond = (count_r == th_diff && off + lgt == prev_off + prev_lgt) ...
                || (count_l >= th_diff-1 && lgt == 1 && off == prev_off);
            it2 = it2 + 1;
        end;

        mat(end-i, off:off+lgt-1) = 1;
        
        if it1 >= maxit || it2 >= maxit
            error 'Infinite loop when generating crossover pattern';
        else
            it1 = 0;
            it2 = 0;
        end;
        
        if off == prev_off
            count_l = count_l + 1;
        else
            count_l = 1;
        end;
        if off+lgt == prev_off+prev_lgt && off+lgt-1 ~= th_side
            count_r = count_r + 1;
        else
            count_r = 1;
        end;

        prev_off = off;
        prev_lgt = lgt;

        % If there is a blank on both sides of the line, reduce max length
        if gen_max_off == th_side && off+lgt-1 < th_side
            gen_max_off = gen_max_off - 1;
        end;
        
        i = i+1;
    end;

    % Quantize
    if quant
        patrn = zeros(side_lg);
        mat = reshape(mat, 1, numel(mat));
        mat = repmat(mat, 2, 1);
        mat = reshape(mat, side_lg, side_lg/step);
        patrn(:, 1:2:end) = mat;
        patrn(:, 2:2:end) = mat;
    else
        patrn = mat;
    end;

    % Vertical or horizontal
    if rand < 0.5
        patrn = patrn';
    end;

end

%% Function to generate clusters
function ptrns = genLayout(side, nb)
    %GENLAYOUT Create random grid matrices
    %
    %   ptrns = GENLAYOUT(side, num)
    %
    %   INPUT:
    %       side:   side length of the matrix
    %       num:    number of matrices to generate
    %   OUTPUT:
    %       ptrns:  cell array containing the generated matrices
    
    if nargin < 2 || isempty(nb)
        nb = 1;
    end;
    if side < 1 || nb < 1
        error('MyERR:Invalid', 'Invalid input argument');
    end;
    
    ptrns = cell(1, nb);
    
    side1 = side/2;
    side2 = side/4;
    
    for i=1:nb
        tmp = zeros(side);
        
        ypos = randi([1 side2]);
        while ypos <= side
            xpos = randi([1 side1]);
            while xpos <= side
                width = min(randi([1 side2]), side-xpos+1);
                height = min(randi([1 side2]), side-ypos+1);

                tmp(ypos:ypos+height-1, xpos:xpos+width-1)=1;
                xpos = xpos+width+randi([1 side1]);
            end;
            ypos = ypos+randi([1 side2]);
        end;
        
        ptrns{i} = tmp;
    end;
end

%% Function to save current state of the algorithm
function fname = save_state(pop, eva, eva_c, iter)
    % Save current state of the optimization algorithm
    %
    % INPUT:
    %   pop:    population cell array OR config data
    %   eva:    fitness vector
    %   eva_c:  fitness vector (of the other fitness mode)
    %   iter:   iteration
    % OUTPUT
    %   fname:  name of subfolder containing the saved data

    if nargin == 1 && length(pop) == numel(pop)
        config = 1;
        iter = 0;
    elseif isempty(eva) && length(pop) == numel(pop)
        config = 1;
    else
        config = 0;
    end;
    if nargin == 3
        iter = eva_c;
        eva_c = [];
    end;

    if ~config
        if numel(pop) ~= numel(eva)
            error 'pop and eva have different size';
        end;
        if numel(pop) ~= length(pop)
            pop = reshape(pop, 1, numel(pop));
        end;
        if numel(eva) ~= length(eva)
            eva = reshape(eva, 1, numel(eva));
        end;
        if ~isempty(eva_c) && numel(eva_c) ~= length(eva_c)
            eva_c = reshape(eva_c, 1, numel(eva_c));
        end;

        if size(eva, 1) ~= 1
            eva = eva';
        end;
        if ~isempty(eva_c) && size(eva_c, 1) ~= 1
            eva_c = eva_c';
        end;
    end;

    persistent prefix_folder_name

    % Determine sav folder
    if iter < 0
        error 'Invalid iter argument';
    elseif ispref('ga_2D', 'save_folder')
        prefix_folder_name = getpref('ga_2D', 'save_folder');
        rmpref('ga_2D', 'save_folder');
    elseif config && iter == 0
        prefix_folder_name = datestr(now, 'yyyymmdd/HHMMSS/');
    end;

    if ~config
        dir = [prefix_folder_name num2str(iter) '/'];
        if ~exist(dir, 'dir')
            mkdir(dir);
        end;

        % Generate individual-fitness array
        pairs = [1:length(pop); eva]';
        save([dir 'fitness.dat'], 'pairs', '-ASCII');
        clearvars pairs

        if ~isempty(eva_c)
            pairs_c = [1:length(pop); eva_c]';
            [tmp_val, pos] = max(pairs_c(:,2));
            pairs_c(pos,:) = pairs_c(2,:);
            pairs_c(2, :) = pairs_c(1, :);
            pairs_c(1, :) = [pos tmp_val];

            save([dir 'fitness_conv.dat'], 'pairs_c', '-ASCII');
            clearvars pairs_c
        end;

        % Save all array arrangements
        for i=1:length(pop)
            arrangmt = AntArray.quantize(pop{i}.M, 1);
            save([dir 'arrangement_' num2str(i) '.dat'], 'arrangmt', '-ASCII');
        end;

        fname = prefix_folder_name(end-6:end-1);
    else
        if ~exist(prefix_folder_name, 'dir')
            mkdir(prefix_folder_name);
        end;
        save([prefix_folder_name 'config.ga_2D'], 'pop', '-ASCII');
    end;
end
