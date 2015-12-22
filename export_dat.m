%
% Copyright 2015, Antoine JUCKLER
%

function export_dat(h, path)
    date_v = datevec(date);
    dirpath = cell(1,3);
    for i=1:3
        dirpath{i} = mat2str(date_v(i));
    end;
    dirpath = sprintf('%s', dirpath{:});
    
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end;
    
    dirpath = [dirpath '/dat'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end;
    
    save([dirpath '/' path '.dat'], 'h', '-ASCII');
end