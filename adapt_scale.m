clear
arr = AntArray(zeros(16));
arr = arr.setName('scale');
arr = arr.adaptArray(ones(16), 50000, 0, 0);
arr.genPattern(1000, 2000, 'YZ', 30);