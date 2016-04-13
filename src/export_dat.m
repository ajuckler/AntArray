%
% Copyright 2015, Antoine JUCKLER
%

function export_dat(h, path, mode)
    if nargin < 3
        mode = 0;
    else
        mode = (mode > 0);
    end;
    
    date_v = datevec(date);
    dirpath = cell(1,3);
    for i=1:3
        dirpath{i} = mat2str(date_v(i));
        if numel(dirpath{i}) < 2
            dirpath{i} = ['0' dirpath{i}];
        end;
    end;
    dirpath = sprintf('%s', dirpath{:});
    
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end;
    
    dirpath = [dirpath '/dat'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end;
    
    if mode == 0
        save([dirpath '/' path '.dat'], 'h', '-ASCII');
    else
        save([dirpath '/' path '.mat'], 'h', '-MAT', '-v7.3');
    end;
end