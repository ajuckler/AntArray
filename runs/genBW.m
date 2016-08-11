function genBW(name, opt)
    if nargin < 2 || isempty(opt)
        opt = 'YZ';
    end;
    ant = AntArray(['dat/elements_' name], 60500, [], .84, [], 0);
    ant = ant.setName(name);
    if strcmpi(opt, 'YZ')
        ant.genPattern(2000, 3000, 'YZ-BW');
    elseif strcmpi(opt, 'XY')
        ant.genPattern(11000, 3000, 'XY-BW');
    elseif strcmpi(opt, 'theta')
        ant.genPattern(11000, 3000, 'theta-BW', [], 0);
    else
        error 'Opt parameter not recognized';
    end;
    % ant.genPattern(11000, 3000, 'XY-BW');
end