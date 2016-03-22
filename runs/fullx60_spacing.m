% Element spacing influence
for space=0.84:0.02:1
    % Init
    arr = AntArray(zeros(60), 60500, [], space);
    arr = arr.setName(['fullx60_' mat2str(space*100)]);
    arr = arr.setMax('XY', 30);
    arr = arr.setMax('YZ', 30);
    arr = arr.setMax('E', 25);
    arr = arr.setMin('XY', -60);
    arr = arr.setMin('YZ', -60);
    arr = arr.setMin('E', -15);

    % Create elements' pattern
    arr = arr.adaptArray(ones(60), 90000, 0, 0);
    
    arr = arr.setComments(sprintf(['Elements spacing: ' mat2str(space) ...
        '$\\lambda$']));

    % Plots
    iter = round((space - 0.84)/0.02 + 1);
    
    arr.plotAntArray();
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    arr = arr.genPattern(10*1000, 3000, 'YZ', 30);
    arr = arr.genPattern(10*1000, [], 'YZ-BW');
end;