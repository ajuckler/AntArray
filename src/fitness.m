function weight = fitness(ant, dist)
    % Gen 2D pattern at certain distance
    % Apply weighting as on draft
    counter_th = 2;
    step = 30;
    const = step*step;
    
    ant = ant.disableWaitbars();
    
    [~, ptrn] = ant.genPattern(dist, 3000, 'YZ', step);
    
    half = floor(length(ptrn)/2);
    ptrn = ptrn(half:end, half:end);
    mask = zeros(length(ptrn));
    
    for i=1:length(ptrn)
        counter = 0;
        ptrn_ln = ptrn(i, :);
        vect = zeros(1, length(ptrn));
        
        if i > 1 && mask(i-1, 1) == 0
            break;
        end;
        
        for j=1:length(vect)
            if ptrn_ln(j) < AntArray.min_E
                counter = counter + 1;
                if counter >= counter_th
                    break;
                end;
            elseif counter ~= 0
                counter = 0;
                vect(max(1, j-counter_th):j) = 1;
            else
                vect(j) = 1;
            end;
        end;
        mask(i,:) = vect;
    end;
    
    ptrn = ptrn - AntArray.min_E;
    ptrn(mask == 0) = 0;
    
    if isempty(ptrn)
        weight = 0;
    else 
        % Compute on symmetric part
        weight = 4*const*sum(sum(ptrn(2:end, 2:end)));

        % Compute on symmetry lines
        weight = weight + const*ptrn(1,1);
        weight = weight + 2*const*sum(ptrn(1, 2:end));
        weight = weight + 2*const*sum(ptrn(2:end, 1));
    end;
end