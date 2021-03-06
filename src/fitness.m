%FITNESS Compute the fitness based on the beam volume
%
%   This function is used by GA_2D to associate a fitness to each array
%   arrangement. 
%   The way of computing the fitness depends on the MODE parameter:
%   For MODE = 0:   The fitness is equal to the volume of the beam that is
%                   above the threshold electric field (as defined in
%                   ANTARRAY).
%   For MODE = 1:   The fitness is equal to the surface formed by the
%                   intersection of the beam and a plane parallel to the
%                   array at the given distance
%
%   WEIGHT = FITNESS(ANT, DIST, MODE)
%   INPUT:
%       ANT:    ANTARRAY object, its YZ pattern will be computed
%       DIST:   distance from the array plane at which the fitness should
%               be evaluated [mm]
%       MODE:   (optional) if set to 1, the fitness will be computed from
%               the surface of the cut plane through the beam; if set to 0
%               it will be computed from the volume of the beam
%               [default = 0]
%   OUTPUT:
%       WEIGHT:     the fitness [Vm | m�]
%       WEIGHT_ALT: fitness using the other mode [m� | Vm]
%
%   See also GA_2D ANTARRAY LOCAL_OPT

%   Copyright 2015-2016, Antoine JUCKLER. All rights reserved.

function [weight, weight_alt] = fitness(ant, dist, mode)
    if nargin < 3
        mode = 0;
    else
        mode = (mode > 0);
    end;
    
    counter_th = 1;
    step = 30;
    const = step*step/1000/1000; % converted to m�
    
    mat = ant.M;
    if numel(mat(mat == 0)) == numel(mat)
        weight = 0;
        weight_alt = 0;
        return;
    end;
    
    ant = ant.waitbars(0);
    
    [~, ptrn] = ant.genPattern(dist, 3000, 'YZ-main', step);
    
    half = ceil(length(ptrn)/2);
    ptrn = ptrn(half:end, half:end);
    mask = zeros(length(ptrn));
    
    for i=1:length(ptrn)
        counter = 0;
        ptrn_ln = ptrn(i, :);
        vect = zeros(1, length(ptrn));
        
        if i > 2 && mask(i-counter_th:i-1, 1) == 0
            break;
        end;
        
        for j=1:length(vect)
            if ptrn_ln(j) < ant.min_E
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
    
    ptrn = ptrn - ant.min_E;
    ptrn = 10.^(ptrn./20);
    ptrn(mask == 0) = 0;
    
    if isempty(ptrn)
        weight = 0;
        weight_alt = 0;
    else
        if mode == 1
            ptrn_alt = ptrn;
            ptrn = mask;
        else
            ptrn_alt = mask;
        end;
        % Compute on symmetric part
        weight = 4*const*sum(sum(ptrn(2:end, 2:end)));
        weight_alt = 4*const*sum(sum(ptrn_alt(2:end, 2:end)));

        % Compute on symmetry lines
        weight = weight + const*ptrn(1,1);
        weight = weight + 2*const*sum(ptrn(1, 2:end));
        weight = weight + 2*const*sum(ptrn(2:end, 1));
        weight_alt = weight_alt + const*ptrn_alt(1,1);
        weight_alt = weight_alt + 2*const*sum(ptrn_alt(1, 2:end));
        weight_alt = weight_alt + 2*const*sum(ptrn_alt(2:end, 1));
    end;
end