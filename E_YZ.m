function [] = E_YZ(antarray, d, L, ss)
% @brief
% Compute the electric field of an antenna array at a certain distance
%
% Compute the electric field of a dipole antenna array operated at
% frequency f, over a plane located at a distance d and with a side-size
% L. The antenna elements are spaced by a distance s and have a length l.
%
% @param    M   Array containing the amplitude and phase of the current
%               on the constituting antenna elements
% @param    f   Frequency of operation in MHz
% @param    l   Antenna element length in mm
% @param    s   Inter-element spacing in mm
% @param    d   Plane distance to the array in mm
% @param    L   Plane side size in mm
% @param    ss  Step-size in mm

% dipole length 2.1826855509101mm

dim = round(L/ss)+1;
plotdata = zeros(dim+1);    % Larger to be able to plot evth
plane = zeros(dim, dim, 3);

progress = waitbar(0, 'Computations in progress...',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(progress, 'canceling', 0);

for i=1:dim
    z = L/2 - (i-1)*ss;
    for j=1:dim
        y = (j-1)*ss - L/2;
        [E_x, E_y, E_z] = E_array(antarray, d, y, z);
        plane(i,j,:) = [E_x, E_y, E_z];
        plotdata(i,j) = 20*log10(sqrt(sum(abs(plane(i,j,:)))));
%         plotdata(i,j) = abs(plane(i,j,3));
    end;
    if getappdata(progress, 'canceling')
        close all;
        delete(progress);
        return;
    end;
    waitbar(i/dim);
end;

waitbar(1, progress, 'Generating plots...');

% ========================================================================
% Plot
% Generate axes
if L/2 < 1/100
    fact = 1000;
elseif L/2 < 1
    fact = 100;
else
    fact = 1;
end;

range = -L/2:ss:L/2';
range = range.*fact;
[absc, oord] = meshgrid([range range(end)+ss*fact]);  % Larger to be able to plot evth

% Plot field
figure(1);
surf(absc, oord, plotdata, 'EdgeColor', 'none');
view(2);
colormap jet;
cbar_h = colorbar('eastoutside');

% Adapt color map
title(cbar_h, 'dB\,V/m', 'Interpreter', 'latex', 'FontSize', 16);
cmax = max(max(plotdata(:,:)));
caxis([0 cmax+5-mod(cmax,5)]);

% Adapt ticks
if mod(L*fact,4) == 0
    tick_fact = 4;
elseif mod(L*fact,6) == 0
    tick_fact = 6;
else
    tick_fact = 2;
end;

spe_ticks = zeros(tick_fact+1,1);
for ii=1:length(spe_ticks)
    spe_ticks(ii) = range(1) + (ii-1)*L*fact/tick_fact;
end;
spe_ticks_pos = spe_ticks+ss*fact/2;

xlim([range(1), range(end)+ss*fact]);
ylim([range(1), range(end)+ss*fact]);
set(gca, 'XTick', spe_ticks_pos, ...
	'XTickLabel', spe_ticks);
set(gca, 'YTick', spe_ticks_pos, ...
    'YTickLabel', spe_ticks);

% Set labels and title
title(['\textbf{Electric field at ', mat2str(d), 'm}'], ...
    'Interpreter', 'latex', 'FontSize', 24);

switch fact
    case 1000
        unit = 'mm';
    case 100
        unit = 'cm';
    otherwise
        unit = 'm';
end;
xlabel(['${\rm y}_{\rm pos}$ [' unit ']'], ...
    'Interpreter', 'latex', 'FontSize', 22);
ylabel(['${\rm z}_{\rm pos}$ [' unit ']'], ...
    'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'FontSize', 16);

if isfield(antarray, 'name')
    name = ['_' antarray.name];
else
    name = '';
end;

print_plots(gcf, ['pattern_' mat2str(d) name]);

delete(progress);

close all
end
