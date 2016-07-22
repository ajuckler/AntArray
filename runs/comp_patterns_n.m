clear

parallel_pool('start');

% Array containing the weights
wArr = zeros(32,1);

fsize = 64;
% Uniform
    % Init
    arr = AntArray(zeros(fsize), 60500, [], 0.84, [], 0);
    arr = arr.setName(['fullx_' mat2str(fsize)]);
    arr = arr.setMax('XY', 30);
    arr = arr.setMax('YZ', 30);
    arr = arr.setMax('E', 25);
    arr = arr.setMin('XY', -60);
    arr = arr.setMin('YZ', -60);
    arr = arr.setMin('E', -15);

    % Create elements' pattern
    arr = arr.adaptArray(ones(fsize), 100000, 0, 0);
    
    % Plots
    arr.plotAntArray();
    arr = arr.genPattern(11000, 3000, 'theta', 30, pi/4);
    arr = arr.genPattern([], [], 'theta-BW', [], pi/4);
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    wArr(2,1) = arr.weight('XY');
    wArr(2,2) = arr.weight('theta', pi/4);
    arr = arr.E_strength(15000, 0, 0, 500);
    arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
    arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');

% Squares
for rem_els=0:15
    disp(['Squares, iteration ' num2str(rem_els) ' of 15']);
    % Init
    arr = AntArray(zeros(fsize), 60500, [], 0.84, [], 0);
    arr = arr.setName(['sq2_' mat2str(rem_els)]);
    arr = arr.setMax('XY', 30);
    arr = arr.setMax('YZ', 30);
    arr = arr.setMax('E', 25);
    arr = arr.setMin('XY', -60);
    arr = arr.setMin('YZ', -60);
    arr = arr.setMin('E', -15);

    % Create elements' pattern
    tmp = gen_sq_arrgt(64, rem_els);
    arr = arr.adaptArray(tmp, 100000, 0, 0);

    % Plots
    arr.plotAntArray();
    arr = arr.genPattern(11000, 3000, 'theta', 30, pi/4);
    arr = arr.genPattern([], [], 'theta-BW', [], pi/4);
    wArr(rem_els+1,3) = arr.weight('theta', pi/4);
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    wArr(rem_els+1,4) = arr.weight('XY');
    arr = arr.E_strength(15000, 0, 0, 500);
    arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
    arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
end;

% Triangles
for rem_els=0:15
    disp(['Triangles, iteration ' num2str(rem_els) ' of 15']);
    % Init
    arr = AntArray(zeros(fsize), 60500, [], 0.84, [], 0);
    arr = arr.setName(['tr_' num2str(rem_els)]);
    arr = arr.setMax('XY', 30);
    arr = arr.setMax('YZ', 30);
    arr = arr.setMax('E', 25);
    arr = arr.setMin('XY', -60);
    arr = arr.setMin('YZ', -60);
    arr = arr.setMin('E', -15);

    % Create elements' pattern
    tmp = gen_tr_arrgt(64, 2, rem_els);
    arr = arr.adaptArray(tmp, 100000, 0, 0);
    
    % Plots
    arr.plotAntArray();
    arr = arr.genPattern(11000, 3000, 'theta', 30, 0);
    arr = arr.genPattern([], [], 'theta-BW', [], 0);
    wArr(rem_els+1,5) = arr.weight('theta', 0);
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    wArr(rem_els+1,6) = arr.weight('XY');
    arr = arr.E_strength(15000, 0, 0, 500);
    arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
    arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
end;

% Triangles 2
for rem_els=0:15
    disp(['Triangles 2, iteration ' num2str(rem_els) ' of 15']);
    % Init
    arr = AntArray(zeros(fsize), 60500, [], 0.84, [], 0);
    arr = arr.setName(['tr2_' num2str(rem_els)]);
    arr = arr.setMax('XY', 30);
    arr = arr.setMax('YZ', 30);
    arr = arr.setMax('E', 25);
    arr = arr.setMin('XY', -60);
    arr = arr.setMin('YZ', -60);
    arr = arr.setMin('E', -15);

    % Create elements' pattern
    tmp = gen_tr_arrgt(64, 4, rem_els);
    arr = arr.adaptArray(tmp, 100000, 0, 0);

    % Plots
    arr.plotAntArray();
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    wArr(rem_els+1,7) = arr.weight('XY');
    arr = arr.E_strength(15000, 0, 0, 500);
    arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
    arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
end;

% Circles
for rem_els=0:31
    disp(['Circles, iteration ' num2str(rem_els) ' of 31']);
    % Init
    arr = AntArray(zeros(fsize), 60500, [], 0.84, [], 0);
    arr = arr.setName(['circ_' num2str(rem_els)]);
    arr = arr.setMax('XY', 30);
    arr = arr.setMax('YZ', 30);
    arr = arr.setMax('E', 25);
    arr = arr.setMin('XY', -60);
    arr = arr.setMin('YZ', -60);
    arr = arr.setMin('E', -15);

    % Create elements' pattern
    tmp = drawCircle(64, 32, rem_els);
    arr = arr.adaptArray(tmp, 100000, 0, 0);
    
    % Plots
    arr.plotAntArray();
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    wArr(rem_els+1,1) = arr.weight('XY');
    arr = arr.E_strength(15000, 0, 0, 500);
    arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
    arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
    arr = arr.genPattern(10*1000, 3000, 'YZ', 30);
    arr = arr.genPattern(10*1000, [], 'YZ-BW');
end;

export_dat(wArr, ['weights_' datestr(now, 'yyyymmdd')]);

parallel_pool('stop');