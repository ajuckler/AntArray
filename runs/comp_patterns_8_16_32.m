clear

parallel_pool('start');

for sz = [8 16 32]
    % Array containing the weights
    wArr = cell(round(sz/2),8);
    % Uniform
    arr = AntArray(zeros(sz), 60500, [], 0.84, [], 0);
    arr = arr.setName(['fullx_' mat2str(sz)]);
    arr = arr.setMax('XY', 30);
    arr = arr.setMax('YZ', 30);
    arr = arr.setMax('E', 25);
    arr = arr.setMin('XY', -60);
    arr = arr.setMin('YZ', -60);
    arr = arr.setMin('E', -15);

    % Create elements' pattern
    arr = arr.adaptArray(ones(sz), 100000, 0, 0);
    
    arr = arr.setComments(sprintf('Elements spacing: 0.84$\\lambda$'));
    
    % Plots
    arr.plotAntArray();
    arr = arr.genPattern(11000, 3000, 'theta', 30, pi/4);
    arr = arr.genPattern([], [], 'theta-BW', [], pi/4);
    wArr{1, 1} = 'Uniform theta';
    wArr{2, 1} = arr.weight('theta', pi/4);
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    wArr{1, 2} = 'Uniform XY';
    wArr{2, 2} = arr.weight('XY');
    arr = arr.E_strength(15000, 0, 0, 500);
    arr = arr.genPattern(2000, 3000, 'YZ', 30);
    arr = arr.genPattern(2000, [], 'YZ-BW');

    for rem_els=0:round(sz/4)-1
        disp(['Squares, iteration ' num2str(rem_els) ' of ' num2str(round(sz/4)-1)]);
        % Init
        arr = AntArray(zeros(sz), 60500, [], 0.84, [], 0);
        arr = arr.setName(['sq2_' mat2str(sz) '_' mat2str(rem_els)]);
        arr = arr.setMax('XY', 30);
        arr = arr.setMax('YZ', 30);
        arr = arr.setMax('E', 25);
        arr = arr.setMin('XY', -60);
        arr = arr.setMin('YZ', -60);
        arr = arr.setMin('E', -15);

        % Create elements' pattern
        tmp = gen_sq_arrgt(sz, rem_els);
        arr = arr.adaptArray(tmp, 100000, 0, 0);

        el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
        arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
            'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));

        % Plots
        arr.plotAntArray();
        arr = arr.genPattern(11000, 3000, 'theta', 30, pi/4);
        arr = arr.genPattern([], [], 'theta-BW', [], pi/4);
        wArr{1, 3} = 'SQ theta';
        wArr{rem_els+2, 3} = arr.weight('theta', pi/4);
        arr = arr.genPattern(11000, 3000, 'XY', 30);
        arr = arr.genPattern([], [], 'XY-BW');
        wArr{1, 4} = 'SQ XY';
        wArr{rem_els+2, 4} = arr.weight('XY');
        arr = arr.E_strength(15000, 0, 0, 500);
        arr = arr.genPattern(2000, 3000, 'YZ', 30);
        arr = arr.genPattern(2000, [], 'YZ-BW');
    end;

    % Triangles
    for rem_els=0:round(sz/4)-1
        disp(['Triangles, iteration ' num2str(rem_els) ' of ' num2str(round(sz/4)-1)]);
        % Init
        arr = AntArray(zeros(sz), 60500, [], 0.84, [], 0);
        arr = arr.setName(['tr_' num2str(sz) '_' num2str(rem_els)]);
        arr = arr.setMax('XY', 30);
        arr = arr.setMax('YZ', 30);
        arr = arr.setMax('E', 25);
        arr = arr.setMin('XY', -60);
        arr = arr.setMin('YZ', -60);
        arr = arr.setMin('E', -15);

        % Create elements' pattern
        tmp = gen_tr_arrgt(sz, 2, rem_els);
        arr = arr.adaptArray(tmp, 100000, 0, 0);

        el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
        arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
            'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));

        % Plots
        arr.plotAntArray();
        arr = arr.genPattern(11000, 3000, 'theta', 30, 0);
        arr = arr.genPattern([], [], 'theta-BW', [], 0);
        wArr{1, 5} = 'TR theta';
        wArr{rem_els+2, 5} = arr.weight('theta', 0);
        arr = arr.genPattern(11000, 3000, 'XY', 30);
        arr = arr.genPattern([], [], 'XY-BW');
        wArr{1, 6} = 'TR XY';
        wArr{rem_els+2, 6} = arr.weight('XY');
        arr = arr.E_strength(15000, 0, 0, 500);
        arr = arr.genPattern(2000, 3000, 'YZ', 30);
        arr = arr.genPattern(2000, [], 'YZ-BW');
    end;

    % Triangles 2
    for rem_els=0:round(sz/4)-1
        disp(['Triangles 2, iteration ' num2str(rem_els) ' of ' num2str(round(sz/4)-1)]);
        % Init
        arr = AntArray(zeros(sz), 60500, [], 0.84, [], 0);
        arr = arr.setName(['tr2_' num2str(sz) '_' num2str(rem_els)]);
        arr = arr.setMax('XY', 30);
        arr = arr.setMax('YZ', 30);
        arr = arr.setMax('E', 25);
        arr = arr.setMin('XY', -60);
        arr = arr.setMin('YZ', -60);
        arr = arr.setMin('E', -15);

        % Create elements' pattern
        tmp = gen_tr_arrgt(sz, 4, rem_els);
        arr = arr.adaptArray(tmp, 100000, 0, 0);

        el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
        arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
            'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));

        % Plots
        arr.plotAntArray();
        arr = arr.genPattern(11000, 3000, 'XY', 30);
        arr = arr.genPattern([], [], 'XY-BW');
        wArr{1, 7} = 'TR2 XY';
        wArr{rem_els+2, 7} = arr.weight('XY');
        arr = arr.E_strength(15000, 0, 0, 500);
        arr = arr.genPattern(2000, 3000, 'YZ', 30);
        arr = arr.genPattern(2000, [], 'YZ-BW');
    end;

    % Circles
    last_it = round(sz/2) - 1;
    for rem_els=0:last_it
        disp(['Circles, iteration ' num2str(rem_els) ' of ' num2str(last_it)]);
        % Init
        arr = AntArray(zeros(sz), 60500, [], 0.84, [], 0);
        arr = arr.setName(['circ_' num2str(sz) '_' num2str(rem_els)]);
        arr = arr.setMax('XY', 30);
        arr = arr.setMax('YZ', 30);
        arr = arr.setMax('E', 25);
        arr = arr.setMin('XY', -60);
        arr = arr.setMin('YZ', -60);
        arr = arr.setMin('E', -15);
    
        % Create elements' pattern
        tmp = drawCircle(sz, round(sz/2), rem_els);
        arr = arr.adaptArray(tmp, 100000, 0, 0);
        
        el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
        arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
            'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));
    
        % Plots
        arr.plotAntArray();
        arr = arr.genPattern(11000, 3000, 'XY', 30);
        arr = arr.genPattern([], [], 'XY-BW');
        wArr{1, 8} = 'Circ XY';
        wArr{rem_els+2, 8} = arr.weight('XY');
        arr = arr.E_strength(15000, 0, 0, 500);
        arr = arr.genPattern(2000, 3000, 'YZ', 30);
        arr = arr.genPattern(2000, [], 'YZ-BW');
    end;
    mat2lat(wArr, ['weights_' num2str(sz)], 'string');
end;

parallel_pool('stop');