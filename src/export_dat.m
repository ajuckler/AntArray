function export_dat(h, path, mode)
    %EXPORT_DAT export data from plots or variable to file
    %
    % [ ] = EXPORT_DAT(h, path, mode)
    %
    % INPUT
    %   h:      handler to plot or variable to save
    %   path:   save name
    %   mode:   (optional) type of save: 0 for plot, 1 for variable
    %           [default = 0]
    
    % Copyright 2015-2016, Antoine JUCKLER. All rights reserved.
    
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