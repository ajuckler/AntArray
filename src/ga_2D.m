%GA_2D Optimize the array arrangement using a genetic algorithm
%   [OPTIM_SOL, OPTIM_VAL] = GA_2D(DIST, START_POP, QUANT)
%
%   The algorithm parameters are hard-coded as follow:
%   chromosome size:    
%       if QUANT:   1/4 of matrix size in START_POP or of AntArray default
%       if ~QUANT:  1/2 of matrix size in START_POP or of AntArray default
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
%   INPUT:
%       DIST:       distance from the array where the fitness function will
%                   be evaluated [mm]
%       START_POP:  (optional) individuals to be included in the population
%                   [cell array of AntArray elements]
%       QUANT:      (optional) is the array arrangement quantized?
%                   [boolean]
%
%   OUTPUT:
%       OPTIM_SOL:  the optimal array arrangement [AntArray]
%       OPTIM_VAL:  fitness of the optimal arrangement [double]
%
%   See also FITNESS GENCROSSOVERPATTERN CHROM2MAT MAT2CHROM

%   Copyright 2016, Antoine Juckler. All rights reserved.

function [optim_sol, optim_val] = ga_2D(dist, start_pop, quant)

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

mut_prob = 0.001;   % mutation probability

max_iter = 50;      % max number of iterations

if mod(pop_sz, trn_sz) ~= 0
    pop_sz = pop_sz + mod(pop_sz, trn_sz);
end;

try
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
                                   AntArray.quantize(start_pop{i}.M, 2, 1));
                end;
            end;
        else
            if ~isa(start_pop, 'AntArray')
                error('START_POP is not of type AntArray');
            elseif quant
                tmp = AntArray(...
                               AntArray.quantize(start_pop.M, 2, 1));
                start_pop = cell(1,1);
                start_pop{1} = tmp;
            else
                tmp = start_pop;
                start_pop = cell(1,1);
                start_pop{1} = tmp;
            end;
        end;
        chrom_sz = size(start_pop{1}.M, 1);
    end;

    % Init wait window
    dial = WaitDialog();
    dial.setMainString('Starting genetic algorithm')

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
            eva_ln(j) = fitness(pop_ln{j}, dist);
        end;
        eva(:, i) = eva_ln;
    end;
    
    dial.terminate();
    dial.setMainString('Initial population evaluated');

    iter = 1;
    progress_data = zeros(max_iter+1, 2);
    save_state(pop, eva, iter);
    dial.setSubString(['Generation ' num2str(iter) ' data saved']);
    dial.terminate();

    while iter <= max_iter
        dial.setMainString(['Working on generation ' num2str(iter) '...']);

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

        % Store max & mean for progress tracking
        progress_data(iter, :) = [eva(1) sum(sum(eva))/numel(eva)];

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
                vals(j) = fitness(inds{j}, dist);
            end;
            eva_tmp(:, i) = vals;
            pop_tmp(:, i) = inds;
        end;

        eva = eva_tmp;
        pop = pop_tmp;

        iter = iter+1;

        save_state(pop, eva, iter);
    end;
    dial.terminate();

    [optim_val, pos] = max(eva(:));
    optim_sol = pop{pos};

    % Store max & mean for progress tracking
    progress_data(iter, :) = [optim_val sum(sum(eva))/numel(eva)];
    delete(dial);
    
catch ME
    switch ME.identifier
        case 'MyERR:Terminated'
            warning('Operation terminated by user');
        otherwise
            rethrow(ME);
    end;
end;

% Plot progress
% -------------
figure(1);
plot(1:length(progress_data), progress_data(:,1), '-b', ...
    'LineWidth', 2, 'DisplayName', 'max');
hold on
plot(1:length(progress_data), progress_data(:,2), '-r', ...
    'LineWidth', 2, 'DisplayName', 'mean');

xlim([1 length(progress_data)]);

L = legend('Location', 'northwest');
set(L, 'Interpreter', 'latex');

xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Fitness', 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'FontSize', 16);

print_plots(gcf, 'fitness');
close all
    
end