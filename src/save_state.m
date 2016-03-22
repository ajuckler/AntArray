function save_state(pop, eva, iter)
% Save current state of the optimization algorithm
%
% INPUT:
%   pop:    population cell array
%   eva:    fitness vector
%   iter:   iteration
%

if numel(pop) ~= numel(eva)
    error 'pop and eva have different size';
end;
if numel(pop) ~= length(pop)
    pop = reshape(pop, 1, numel(pop));
end;
if numel(eva) ~= length(eva)
    eva = reshape(eva, 1, numel(eva));
end;

if size(eva, 1) ~= 1
    eva = eva';
end;

persistent prefix_folder_name

% Determine sav folder
if iter < 1
    error 'Invalid iter argument';
elseif iter == 1
    prefix_folder_name = datestr(now, 'yyyymmdd/HHMMSS/');
end;

dir = [prefix_folder_name num2str(iter) '/'];
if ~exist(dir, 'dir')
    mkdir(dir);
end;

% Generate individual-fitness array
pairs = [1:length(pop); eva]';
save([dir 'fitness.dat'], 'pairs', '-ASCII');
clearvars pairs

% Save all array arrangements
for i=1:length(pop)
    arrangmt = AntArray.quantize(pop{i}.M, 1);
    save([dir 'arrangement_' num2str(i) '.dat'], 'arrangmt', '-ASCII');
end;

end