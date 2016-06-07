%REPEAT_ARRAY generates field patterns of repeated array at different array
%separations
%
%   [] = REPEAT_ARRAY(SRC, MIN, MAX)
%   INPUT:
%       SRC     Source file name
%       MIN     (optional) minimum separation in terms of elements
%               [default=0]
%       MAX     (optional) maximum separation in terms of elements
%               [default=16]
%

% Copyright 2016 Antoine Juckler. All rights reserved.

function repeat_array(src, min, max)
    if nargin < 3 || isempty(max)
        max = 16;
    else
        max = 2*round(max/2);
    end;
    if nargin < 2 || isempty(min)
        min = 0;
    else
        min = 2*round(min/2);
    end;
    
    parallel_pool('start');
    
    ant = AntArray(src);
    src_arr = ant.M;
    
    pos = strfind(src, '_');
    name = src(pos(end)+1:end);
    
%     ant = AntArray(src_arr, 60500, [], .84);
%     ant = ant.setName(name);
%     ant.plotAntArray(1);
%     ant = ant.genPattern(2000, 3000, 'YZ', 30);
%     ant = ant.genPattern(2000, 3000, 'YZ-BW');
%     ant = ant.genPattern(2000, 3000, 'YZ-main', 30);
%     ant = ant.genPattern(11000, 3000, 'XY', 30);
%     ant = ant.genPattern(11000, 3000, 'XY-BW', 30);
%     ant.E_strength(11000);
    
    for i=min:2:max
        if i < 0
            len = length(src_arr)+i/2;
            arr = [src_arr(1:len, 1:len) src_arr(1:len, end-len+1:end);
                  src_arr(end-len+1:end, 1:len) src_arr(end-len+1:end, end-len+1:end)];
        else
            len = length(src_arr);
            arr = [src_arr zeros(len, i) src_arr;
                zeros(i, i+2*len);
                src_arr zeros(len, i) src_arr];
        end;
        
        ant = AntArray(arr, 60500, [], .84);
        ant = ant.setName([name '_rep_' num2str(i)]);
        ant.plotAntArray(1);
        ant = ant.genPattern(2000, 3000, 'YZ', 30);
        ant = ant.genPattern(2000, 3000, 'YZ-BW');
        ant = ant.genPattern(2000, 3000, 'YZ-main', 30);
        ant = ant.genPattern(11000, 3000, 'XY', 30);
        ant = ant.genPattern(11000, 3000, 'XY-BW', 30);
        ant.E_strength(11000);
    end;
    
    parallel_pool('stop');
end