function res = parallel_pool(mode)
    %PARALLEL_POOL starts, stops or get the number of workers in a parallel
    %pool
    %
    %   RES = PARALLEL_POOL(MODE)
    %   INPUT:
    %       MODE:   'start' to start the pool
    %               'stop' to stop the pool
    %               'status' to get number of workers
    %               If the pool was already started/stopped, nothing will
    %               happen
    %   OUTPUT:
    %       RES:    depends on the inputted mode
    %               if mode is 'start' or 'stop':   0 for success
    %                                               -1 for error
    %               if mode is 'status':            number of workers

    % Copyright 2015-2016, Antoine JUCKLER. All rights reserved
    
    persistent poolobj;

    if strcmp(mode, 'start')
        if verLessThan('matlab','8.2')
            if ~matlabpool('size')
                matlabpool open
            end;
        else
            if isempty(gcp('nocreate'))
                poolobj = parpool;
            elseif isempty(poolobj)
                poolobj = gcp;
            end;
        end;
        res = 0;
    elseif strcmp(mode, 'stop')
        if verLessThan('matlab','8.2') 
            if matlabpool('size')
                matlabpool close
            end;
        else
            if ~isempty(gcp('nocreate'))
                delete(poolobj);
                poolobj = [];
            end;
        end;
        res = 0;
    elseif strcmp(mode, 'status')
        if verLessThan('matlab', '8.2')
            res = matlabpool('size');
        else
            if ~isempty(gcp('nocreate'))
                res = poolobj.NumWorkers;
            else
                res = 1;
            end;
        end;
    else
        poolobj = [];
        res = -1;
        error 'Invalid mode parameter';
    end;

end