function print_plots(h, path) % h must be gcf from plot
    %PRINT_PLOTS exports the plots to pdf
    %
    % The function will create both a pdf and a fig file.
    %
    % [ ] = PRINT_PLOTS(h, path)
    %
    % INPUT
    %   h:      handler to plot object
    %   path:   desired save name
    
    % Copyright 2015-2016, Antoine JUCKLER. All rights reserved.
    
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
        print(h, [dirpath '/' path '.pdf'], '-dpdf', '-painters');
    else
        print(h, [dirpath '/' path '.pdf'], '-dpdf', '-opengl', '-r300');
    end
    
    
    if ~strcmp(path(end-1:end), 'BW')
        dirpath = [dirpath '/fig'];
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end;
        
        saveas(h, [dirpath '/' path '.fig']);
    end;
end