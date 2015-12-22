function [] = E_XY(antarray, d, L, ss)
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
% @param    d   Plane depth in mm
% @param    L   Plane side size in mm
% @param    ss  Step-size in mm

% dipole length 2.1826855509101mm

dim1 = round(L/ss)+1;
dim2 = round(d/ss)+1;
plotdata = zeros(dim2+1, dim1+1);    % Larger to be able to plot evth
plane = zeros(dim2, dim1, 3);

progress = waitbar(0, 'Computations in progress...',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(progress, 'canceling', 0);

for i=1:dim2
    x = (i-1)*ss;
    for j=1:dim1
        y = (j-1)*ss - L/2;
        [E_x, E_y, E_z] = E_array(antarray, x, y, 0);
        plane(i,j,:) = [E_x, E_y, E_z];
        plotdata(i,j) = 20*log10(sqrt(sum(abs(plane(i,j,:)))));
%         plotdata(i,j) = abs(plane(i,j,3));
    end;
    if getappdata(progress, 'canceling')
        close all;
        delete(progress);
        return;
    end;
    waitbar(i/dim2);
end;

waitbar(1, progress, 'Generating plots...');

% ========================================================================
% Plot 1 (Field 2D)
% Generates axes
if L/2 < 1/100
    fact_L = 1000;
elseif L/2 < 1
    fact_L = 100;
else
    fact_L = 1;
end;

if d < 1/100
    fact_d = 1000;
elseif d < 1
    fact_d = 100;
else
    fact_d = 1;
end;

range_L = -L/2:ss:L/2';
range_L = range_L.*fact_L;
range_d = 0:ss:round(d/ss)*ss;
range_d = range_d.*fact_d;
[absc, oord] = meshgrid([range_L range_L(end)+ss*fact_L], ...
    [range_d range_d(end)+ss*fact_d]);  % Larger to be able to plot evth

% Plot field
figure(1);
surf(absc, oord, plotdata, 'EdgeColor', 'none');
view(2);
colormap jet;
cbar_h = colorbar('eastoutside');

% Adapt color map
title(cbar_h, 'dB\,V/m', 'Interpreter', 'latex', 'FontSize', 16);
ln_nb = ceil(10e-2/ss);
max_val = max(plotdata(ln_nb,:));
caxis([0 max_val+5-mod(max_val,5)]);

% Adapt ticks
if mod(L*fact_L,4) == 0
    tick_fact_L = 4;
elseif mod(L*fact_L,6) == 0
    tick_fact_L = 6;
else
    tick_fact_L = 2;
end;

if mod(d*fact_d,4) == 0
    tick_fact_d = 4;
elseif mod(d*fact_d,6) == 0
    tick_fact_d = 6;
else
    tick_fact_d = 2;
end;

spe_ticks_L = zeros(tick_fact_L+1,1);
for ii=1:length(spe_ticks_L)
    spe_ticks_L(ii) = range_L(1) + (ii-1)*L*fact_L/tick_fact_L;
end;
spe_ticks_L_pos = spe_ticks_L+ss*fact_L/2;

spe_ticks_d = zeros(tick_fact_d+1,1);
for ii=1:length(spe_ticks_d)
    spe_ticks_d(ii) = range_d(1) + (ii-1)*d*fact_d/tick_fact_d;
end;
spe_ticks_d_pos = spe_ticks_d+ss*fact_d/2;

xlim([range_L(1), range_L(end)+ss*fact_L]);
ylim([range_d(1), range_d(end)+ss*fact_d]);
set(gca, 'XTick', spe_ticks_L_pos, ...
	'XTickLabel', spe_ticks_L);
set(gca, 'YTick', spe_ticks_d_pos, ...
    'YTickLabel', spe_ticks_d);

% Set labels and title
title(['\textbf{Electric field in the XY plane}'], ...
    'Interpreter', 'latex', 'FontSize', 24);

switch fact_L
    case 1000
        unit_L = 'mm';
    case 100
        unit_L = 'cm';
    otherwise
        unit_L = 'm';
end;
switch fact_d
    case 1000
        unit_d = 'mm';
    case 100
        unit_d = 'cm';
    otherwise
        unit_d = 'm';
end;
xlabel(['${\rm y}_{\rm pos}$ [' unit_L ']'], 'Interpreter', ...
    'latex', 'FontSize', 22);
ylabel(['${\rm x}_{\rm pos}$ [' unit_d ']'], 'Interpreter', ...
    'latex', 'FontSize', 22);
set(gca, 'FontSize', 16);

if isfield(antarray, 'name')
    name = ['_' antarray.name];
else
    name = '';
end;
print_plots(gcf, ['pattern_XY' name]);

close(gcf);

% ========================================================================
% Plot 2 (Field strength at beam centre)
figure(2);
plot(range_d, plotdata(1:end-1,round((dim2+1)/2)), 'LineWidth', 2);
xlim([0 d*fact_d]);
view_axis = axis;
ylim([view_axis(3) max_val+5-mod(max_val,5)]);

title('\textbf{Field strength at the beam centre}', ...
    'Interpreter', 'latex', 'FontSize', 24);
xlabel(['Distance to the array centre [' unit_d ']'], ...
    'Interpreter', 'latex', 'FontSize', 22);
ylabel(['Field strength [dB\,V/m]'], ...
    'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'FontSize', 16);

print_plots(gcf, ['E_strength' name]);

delete(progress);

close all;
end
