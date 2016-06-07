clear

% Squares 2
rem_els = 9;
    arr = AntArray(zeros(64), 60500, [], 0.84);
    arr = arr.setName(['sq2_' mat2str(rem_els)]);

    % Create elements' pattern
    tmp = gen_sq_arrgt(64, rem_els);

    arr = arr.adaptArray(tmp, 100000, 0, 0);
    arr.plotAntArray();

% Triangles
rem_els = 4;
    arr = AntArray(zeros(64), 60500, [], 0.84);
    arr = arr.setName(['tr_64_' num2str(rem_els)]);

    % Create elements' pattern
    mat = gen_tr_arrgt(64, 2, rem_els);

    arr = arr.adaptArray(mat, 100000, 0, 0);

    arr.plotAntArray();

% Triangles 2
rem_els = 9;
    arr = AntArray(zeros(64), 60500, [], 0.84);
    arr = arr.setName(['tr2_64_' num2str(rem_els)]);

    % Create elements' pattern
    mat = gen_tr_arrgt(64, 4, rem_els);

    arr = arr.adaptArray(mat, 100000, 0, 0);
    arr.plotAntArray();


% Circles
rem_els = 22;
    arr = AntArray(zeros(64), 60500, [], 0.84);
    arr = arr.setName(['circ_' num2str(rem_els)]);

    % Create elements' pattern
    tmp = drawCircle(64, 32, rem_els);
    arr = arr.adaptArray(tmp, 100000, 0, 0);
    arr.plotAntArray();