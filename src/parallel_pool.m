%PARALLEL_POOL starts or stops parallel pool
%
%   [] = PARALLEL_POOL(MODE)
%   INPUT:
%       MODE:   'start' to start the pool
%               'stop' to stop the pool
%               If the pool was already started/stopped, nothing will
%               happen

% Copyright 2016, Antoine Juckler. All rights reserved

function parallel_pool(mode)
    persistent poolobj;

    if strcmp(mode, 'start')
        if verLessThan('matlab','8.2')
            if ~matlabpool('size')
                matlabpool open
            end;
        else
            if isempty(gcp('nocreate'))
                poolobj = parpool;
            end;
        end;
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
    else
        poolobj = [];
        error 'Invalid mode parameter';
    end;

end