%STABILITY Controls the stability of an antenna array arrangement
%
%   This function permits to evaluate the stability of an antenna array
%   arrangement. The array will be subject to perturbations that will turn
%   on/off elements with a certain probability. The fitness of the
%   perturbed array is then evaluated and compared to the one of the
%   original array.
%
%   [] = STABILITY(ANT, PROB, DIST, MODE)
%   INPUT:
%       ANT:    ANTARRAY object, its stability will be evaluated
%       PROB:   probability of the occurence of a perturbation in the array
%               arrangement
%       DIST:   distance from the array plane at which the fitness should
%               be evaluated [mm]
%       MODE:   (optional) if set to 1, the fitness will be computed from
%               the surface of the cut plane through the beam; if set to 0
%               it will be computed from the volume of the beam
%               [default = 0]
%
%   See also FITNESS ANTARRAY

%   Copyright 2016, Antoine Juckler. All rights reserved.

function stability(ant, prob, dist, mode)
if prob >= 1
    error 'PROB should be < 1';
end;
if nargin < 4
    mode = 0;
else
    mode = (mode > 0);
end;

max_iter = 25;

% Start parallel pool
persistent poolobj

if verLessThan('matlab','8.2')
    if ~matlabpool('size')
        matlabpool open
    end;
else
    if isempty(gcp('nocreate'))
        poolobj = parpool;
    end;
end;

plotdata = zeros(1, max_iter+1);
plotdata(end) = fitness(ant, dist, mode);

mat = ant.M;

parfor i = 1:max_iter
    % Generate noise matrix
    perturbation = rand(size(mat,1));
    perturbation = (perturbation <= prob);
    
    tmp = mat;
    tmp(perturbation == 1) = abs(tmp(perturbation == 1) - 1);
    
    plotdata(i) = fitness(AntArray(tmp), dist, mode);
end;

plotdata = round(plotdata.*10^4)./10^4;

% Stop parallel pool
if verLessThan('matlab','8.2') 
    if matlabpool('size')
        matlabpool close
    end;
else
    if ~isempty(gcp('nocreate'))
        delete(poolobj);
        poolobj = [];
    end;
end;

figure(1);
histogram(plotdata);
hold on;
fig_dim = axis;
plot([plotdata(end) plotdata(end)], [fig_dim(3) fig_dim(4)], ...
    '--r', 'Linewidth', 1);
title('Stability');
ylabel('Count');
xlabel('Fitness');
hold off;

savname = ['stability_' ant.name];
print_plots(gcf, savname);
export_dat(plotdata, savname);

close all;

avg = mean(plotdata);
std_dev = std(plotdata);
med = median(plotdata);
disp(['Mean: ' avg]);
disp(['Standard deviation: ' std_dev]);
disp(['Median: ' med]);

end
