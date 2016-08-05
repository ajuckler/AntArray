%GENCROSSOVERPATTERN generates mask for chromosome 1-point crossover
%
%   This function generates a mask that will be used for chromosome
%   crossover. The mask correspond to 1-point crossover in 2D: when the
%   full matrix is recontructed from the chromosome, the mask will delimit
%   2 contiguous regions.
%
%   PATTERN = GENCROSSOVERPATTERN(SIDE_LEN, QUANT)
%   INPUT:
%       SIDE_LEN:   side length of the 2D chromosome
%       QUANT:      (optional) is the chromosome quantized?
%   OUTPUT:
%       PATTERN:    matrix containing the mask for 2D 1-point crossover
%
%   See also GA_2D

%   Copyright 2016, Antoine Juckler. All rights reserved.

function patrn = genCrossoverPattern(side_lg, quant)
    if nargin < 2
        quant = 1;
    end;

    if mod(side_lg, 2) ~= 0
        error 'side_lg is not an even number';
    elseif side_lg < 5
        error 'side_lg must be greater than 4';
    end;

    % Theoretical size of matrix (effective elements generated)
    if quant
        step = 2;
    else
        step = 1;
    end;

    th_side = side_lg/step;
    th_diff = ceil(th_side/8);
    mat = zeros(th_side);
    gen_max_off = th_side;
    
    maxit = 20;

    % Generate theoretical pattern
    prev_off = randi([2 th_side]);
    prev_lgt = th_side - prev_off + 1;
    mat(end, prev_off:end) = 1;
    count_l = 1;
    count_r = 1;
    it1 = 0;
    it2 = 0;
    
    maxi = randi([th_diff th_side-2]);
    i = 1;
    while i <= maxi
        min_off = prev_off - th_diff;
        min_off = max(min_off, 2);
        
        max_off = prev_off + th_diff;
        if prev_lgt <= th_diff
            max_off = min(max_off, prev_off + prev_lgt - 1);
        end;
        max_off = min(max_off, gen_max_off);
        
        off = randi([min_off, max_off]);
        cond = (count_l == th_diff && off == prev_off) ...
            || (prev_lgt == 1 && count_r == th_diff && off == prev_off + prev_lgt) ...
            || (off == gen_max_off ...
                && (count_r >= th_diff-1 || count_l >= th_diff-1)) ...
            || (off == prev_off && prev_lgt == 2 ...
                && count_r >= th_diff-1 && count_l >= th_diff-1 ...
                && prev_off+prev_lgt == gen_max_off+1);
        while cond && it1 < maxit
            off = randi([min_off, max_off]);
            cond = (count_l == th_diff && off == prev_off) ...
            || (prev_lgt == 1 && count_r == th_diff && off == prev_off + prev_lgt) ...
            || (off == gen_max_off ...
                && (count_r >= th_diff-1 || count_l >= th_diff-1)) ...
            || (off == prev_off && prev_lgt == 2 ...
                && count_r >= th_diff-1 && count_l >= th_diff-1 ...
                && prev_off+prev_lgt == gen_max_off+1);
            it1 = it1 + 1;
        end;
        
        max_lgt = prev_off + prev_lgt - off + th_diff;
        max_lgt = min(gen_max_off - off + 1, max_lgt);
        
        min_lgt = prev_off + prev_lgt - th_diff - off;
        if off < prev_off
            min_lgt = max(prev_off - off + 1, min_lgt);
        end;
        min_lgt = max(1, min_lgt);
        
        lgt = randi([min_lgt, max_lgt]);
%         cond = (count_r == th_diff && off + lgt == prev_off + prev_lgt) ...
%             || (th_side ~= gen_max_off && count_r == th_diff-1 && ...
%                 lgt == 1 && off == gen_max_off);
        cond = (count_r == th_diff && off + lgt == prev_off + prev_lgt) ...
                || (count_l >= th_diff-1 && lgt == 1 && off == prev_off);
        while gen_max_off ~= th_side && cond && it2 < maxit
            lgt = randi([min_lgt, max_lgt]);
%             cond = (count_r == th_diff && off + lgt == prev_off + prev_lgt) ...
%                 || (th_side ~= gen_max_off && count_r == th_diff-1 && ...
%                     lgt == 1 && off == gen_max_off);
            cond = (count_r == th_diff && off + lgt == prev_off + prev_lgt) ...
                || (count_l >= th_diff-1 && lgt == 1 && off == prev_off);
            it2 = it2 + 1;
        end;
        
%         lgt = randi(max_lgt-1);             % Length
% 
%         max_off = prev_off + prev_lgt - 1;  % Max offset for contiguous pattern
%         max_off = min(max_lgt-lgt+1, max_off);
% 
%         min_off = prev_off - lgt + 1;       % Min offset for contiguous pattern
%         min_off = max(min_off, 2);
% 
%         off = randi([min_off max_off]);     % Offset

        mat(end-i, off:off+lgt-1) = 1;
        
        if it1 >= maxit || it2 >= maxit
            disp(['Iteration left: ' num2str(it1) '    Iteration right: ' num2str(it2)]);
            disp(mat);
            error 'Infinite loop when generating crossover pattern';
        else
            it1 = 0;
            it2 = 0;
        end;
        
        if off == prev_off
            count_l = count_l + 1;
        else
            count_l = 1;
        end;
        if off+lgt == prev_off+prev_lgt && off+lgt-1 ~= th_side
            count_r = count_r + 1;
        else
            count_r = 1;
        end;

        prev_off = off;
        prev_lgt = lgt;

        % If there is a blank on both sides of the line, reduce max length
        if gen_max_off == th_side && off+lgt-1 < th_side
            gen_max_off = gen_max_off - 1;
        end;
        
        i = i+1;
    end;

    % Quantize
    if quant
        patrn = zeros(side_lg);
        mat = reshape(mat, 1, numel(mat));
        mat = repmat(mat, 2, 1);
        mat = reshape(mat, side_lg, side_lg/step);
        patrn(:, 1:2:end) = mat;
        patrn(:, 2:2:end) = mat;
    else
        patrn = mat;
    end;

    % Vertical or horizontal
    if rand < 0.5
        patrn = patrn';
    end;

end