function [optim_sol, optim_val] = ga_2D(start_pop)
% Performs optimization of antenna array arrangement using a genetic
% algorithm
%
% INPUT:
%

if nargin < 1
    start_pop = [];
end;

% Parameters
% ----------
chrom_sz = [];      % chromosome size
                    %   empty = 1/4 of start_pop or of AntArray default

pop_sz = 50;        % population size
trn_sz = 2;         % tournament selection size
fit_sz = trn_size;  % number of elements passed through
                    %   must be a multiple of trn_sz

mut_prob = 0.01;    % mutation probability

max_iter = 50;      % max number of iterations

if mod(pop_sz, trn_sz) ~= 0
    pop_sz = pop_sz + mod(pop_sz, trn_sz);
end;

% Generate population
% -------------------
pop = cell(1, pop_sz);
if ~isempty(start_pop) % if some elements already given, check
    for i=1:length(start_pop)
        if ~isa(start_pop{i}, 'AntArray')
            error('START_POP is not of type AntArray');
        end;
    end;
    chrom_sz = size(start_pop{1}.M, 1);
end;

if isempty(chrom_sz)
    temp = AntArray();
    chrom_sz = size(temp.M, 1);
    clearvars temp;
end;

if mod(chrom_sz, 4) == 0
    chrom_sz = chrom_sz/4;
end;

end_index = min(length(start_pop), pop_sz);
pop(1:end_index) = start_pop(1:end_index);
for i=end_index:pop_sz
    chrom = randi([0 1], 1, chrom_sz*chrom_sz);
    pop{i} = AntArray(chrom2mat(chrom));
end;

% Reshape for future parallel implementation
v_dim = pop_sz/trn_sz;
pop = reshape(pop, [v_dim, trn_sz]);
eva = zeros(v_dim, trn_sz);

% First evaluation
% ----------------
parfor i=1:v_dim
    pop_ln = pop(i,:);
    eva_ln = eva(i,:);
    for j=1:trn_sz
        eva_ln(j) = evaluate(pop_ln{j});
    end;
    eva(i,:) = eva_ln;
end;

% PRINT AND SAVE CURRENT STATE
% gen folder per iteration with random_name
% save arr with number-val pair, and array arrangements
iter = 1;
save_state(pop, eva, iter);

while iter <= max_iter
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
    
    parfor i=fit_sz/trn_sz+1:v_dim
        % Selection
        % ---------
        inds = pop(i,:);
        vals = eva(i,:);
        for j=1:trn_sz
            rand_index = start_index+j;
            while rand_index == (i-1)*trn_sz+j
                rand_index = randi(pop_sz, 1);
            end;
            
            if vals(j) < eva(rand_index)
                inds{j} = pop{rand_index};
            end;
            
            % Convert to chromosomes
            inds{j} = mat2chrom(inds{j}.M);
        end;
        
        % Crossover
        % ---------
        cross_patrn = genCrossoverPattern(chrom_sz);
        chroms = cell(1, trn_sz);
        if trn_sz ~= 2
            for j=1:trn_sz
                temp_chrom = inds{j};   % Parent 1
                swap_part = randi([0 1],1);
                spouse = j;
                while spouse == j
                    spouse = randi(trn_sz, 1);  % Parent 2
                end;
                temp_chrom(cross_patrn == swap_part) = ...
                    spouse(cross_patrn == swap_part);   % Child
                chroms{j} = temp_chrom;
            end;
            clearvars swap_part spouse
        else
            par1 = inds{1};     % Parents
            par2 = inds{2};
            temp_chrom = par1;
            temp_chrom(cross_patrn == swap_part) = ...
                par2(cross_patrn == swap_part);     % Child 1
            chroms{1} = temp_chrom;
            temp_chrom = par2;
            temp_chrom(cross_patrn == swap_part) = ...
                par1(cross_patrn == swap_part);     % Child 2
            chroms{2} = temp_chrom;
        end;
        clearvars temp_chrom cross_patrn;
        
        % Mutation
        % --------
        if rand <= mut_prob
            pos = randi(chrom_sz*chrom_sz, 1);  % Mutation position
            chrom_nb = randi([1 trn_sz], 1);    % Affected chromosome
            mut_chrom = chroms{chrom_nb};
            mut_chrom(pos) = ~mut_chrom(pos);
            chroms{chrom_nb} = mut_chrom;
            clearvars mut_chrom chrom_nb pos
        end;
        
        % Save new pop & evaluate
        % -----------------------
        for j=1:trn_sz
            inds{j} = AntArray(chrom2mat(chroms{j}));
            vals(j) = evaluate(inds{j});
        end;
        eva_tmp(i,:) = vals;
        pop_tmp(i,:) = inds;
    end;
    
    eva = eva_tmp;
    pop = pop_tmp;
    
    iter = iter+1;
    
    save_state(pop, eva, iter)
end;

[optim_val, pos] = max(eva(:));
optim_sol = pop{pos};

    
end