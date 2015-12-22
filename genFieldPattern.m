function [] = genFieldPattern(antarray, d, L, mode, ss)

% @param    M   Array containing the amplitude and phase of the current
%               on the constituting antenna elements
% @param    f   Frequency of operation in MHz
% @param    l   Antenna element length in mm
% @param    s   Inter-element spacing in mm
% @param    d   Plane distance to the array in mm
% @param    L   Plane side size in mm
% @param    ss  Step-size in mm

L = L/1000;
d = d/1000;

if nargin < 5
    ss = antarray.spacing;
else
    ss = ss/1000;
    step_num = round(L/ss);
    if mod(step_num, 2) ~= 0
        step_num = step_num + 1;
    end;
    ss = L/step_num;
end;

disp(['Step size: ' mat2str(ss*1000) 'mm']);
disp(['Array size: ' mat2str(length(antarray.M)*antarray.spacing*1000) 'mm']);

if ~isfield(antarray, 'M')
    error('Missing array field: M');
elseif ~isfield(antarray, 'freq')
    error('Missing array field: freq');
elseif ~isfield(antarray, 'el_len')
    error('Missing array field: el_len');
elseif ~isfield(antarray, 'spacing')
    error('Missing array field: spacing');
end;

if strcmp(mode, 'YZ')
    E_YZ(antarray, d, L, ss);
elseif strcmp(mode, 'XY')
    E_XY(antarray, d, L, ss);
else
    error('Unhandled mode');
end;


