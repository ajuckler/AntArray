%FITNESS Compute the fitness based on the beam volume
%
%   This function is used by GA_2D to associate a fitness on each array
%   arrangement. The fitness is equal to the volume of the beam that is
%   above the threshold electric field (as defined in ANTARRAY).
%
%   WEIGHT = FITNESS(ANT, DIST)
%   INPUT:
%       ANT:    ANTARRAY object, its YZ pattern will be computed
%       DIST:   distance from the array plane at which the fitness should
%               be evaluated [mm]
%   OUTPUT:
%       WEIGHT: the fitness [Vm]
%
%   See also GA_2D ANTARRAY

%   Copyright 2016, Antoine Juckler. All rights reserved.

function weight = fitness(ant, dist)
    counter_th = 2;
    step = 30;
    const = step*step/1000/1000; % converted to m²
    
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