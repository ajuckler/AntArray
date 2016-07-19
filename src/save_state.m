function fname = save_state(pop, eva, eva_c, iter)
% Save current state of the optimization algorithm
%
% INPUT:
%   pop:    population cell array OR config data
%   eva:    fitness vector
%   eva_c:  fitness vector (of the other fitness mode)
%   iter:   iteration
%

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
        pairs_c(pos,:) = pairs_c(1,:);
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