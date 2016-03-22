%
% Copyright 2015, Antoine JUCKLER
%

function print_plots(h, path) % h must be gcf from plot
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
    
    set(h,'PaperOrientation','landscape');
    set(h,'PaperPosition',[0.25 0.25 28.41 19.72]);	%this is for A4
    
    if verLessThan('matlab','8.4')
        print([dirpath '/' path '.pdf'], '-dpdf', '-painters');
    else
        print([dirpath '/' path '.pdf'], '-dpdf', '-opengl', '-r300');
    end
    
    
    if ~strcmp(path(end-1:end), 'BW')
        dirpath = [dirpath '/fig'];
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end;
        
        saveas(h, [dirpath '/' path '.fig']);
    end;
end