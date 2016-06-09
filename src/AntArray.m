%AntArray  Class for handling custom array of dipole antennas
%
%This class permits to create custom arrays of dipole antennas, all placed
%along the Z-direction on a grid with constant spacing centred at the
%origin and located in the YZ-plane.
%Elements of the array can be tuned as to focus the beam on a certain point
%of space, and their amplitude can be further tuned as to realise a
%specific input power over the entire array.
%Several functions are provided to compute the electric field distribution
%on different cut planes or along a line parallel to the X-axis.
%
%AntArray Properties
%   M               matrix containing element's excitations
%   freq            frequency of operation [Hz]
%   norm_freq       frequency used for power normalization [Hz]
%   el_len          dipole's length
%   spacing         inter-elements spacing
%   name            name for the figures
%   min_XY          min dB scale value for XY-patterns
%   max_XY          max dB scale value for XY-patterns
%   min_YZ          min dB scale value for YZ-patterns
%   max_YZ          max dB scale value for YZ-patterns
%   min_E_strength  min dB scale value for E-strength patterns
%   max_E_strength  max dB scale value for E-strength patterns
%   dire            matrix of element's groupings
%   dire_str        cell array of element's groupings
%   comments        comment string, will be printed on elements' plot
%   pwr             input power [W]
%   weight_ang      weight angle [rad]
%   normalized      has the input power been normalized (bool)
%   plotres         save the resulting plots (bool)
%   dispwait        display the waitbars (bool)
%
%AntArray Methods
%   AntArray        constructor
%   setName         set the name property
%   setMin          set the min value for XY, YZ or E-strength plots
%   setMax          set the max value for XY, YZ or E-strength plots
%   setWeightAngle  set the weight angle
%   setComments     set the comments
%   setNormFreq     set the frequency used at normalization
%   setNormPwr      set the input power
%   waitbars        turn waitbars on/off
%   adaptArray      add elements to the array, focused at a specific point
%   adaptAmp        adapt elements' amplitude according to given profile
%   rstArray        clear the elements' array
%   genPattern      generate electric field distribution along a plane
%   E_strength      generate the electric field strength on a line
%   weight          compute the weight
%   plotWeight      plot the weight pattern
%   directivity     compute the directivity
%   directivity_alt compute the directivity (alternative function)
%   plotAntArray    plot the elements' pattern
%   getVal          [static] get field value at given position
%   quantize        [static] quantize a matrix
%
%Use the DOC command for detailed explanations

% Copyright 2015-2016, Antoine JUCKLER. All rights reserved

classdef AntArray
    properties (Constant, Access='public')
        c0 = 299792458; % Speed of light in free space
        Z0 = 119.9169832*pi; % Impedance of free space
        opt_pts = 51;	% Number of discretization steps in one dimension
        min_E = -1.3;   % Minimal electric field for reception [dB V/m]
    end
    properties (GetAccess='public', SetAccess='private')
        M;              % matrix containing element's excitations
        freq;           % frequency of operation [Hz]
        norm_freq;      % frequency used for power normalization [Hz]
        el_len;         % dipole's length
        spacing;        % inter-elements spacing
        opt_win;        % side length and depth of optimization window
        name;           % name used for the figures
        min_XY;         % min dB scale value for XY-patterns
        max_XY;         % max dB scale value for XY-patterns
        min_YZ;         % min dB scale value for YZ-patterns
        max_YZ;         % max dB scale value for YZ-patterns
        min_E_strength; % min dB scale value for E-strength patterns
        max_E_strength; % max dB scale value for E-strength patterns
        dir;            % matrix of element's groupings
        dir_str;        % cell array of element's groupings
        comments;       % comment string, will be printed on elements' plot
        pwr;            % input power [W]
        weight_ang      % weight angle [rad]
        normalized      % has the input power been normalized (bool)
        plotres = 1     % save the resulting plots (bool)
        dispwait = 1    % display the waitbars (bool)
    end
    methods
        %% Constructor
        function obj = AntArray(M, f, l, s)
            %ANTARRAY constructor
            %
            % obj = ANTARRAY(M, F, L, S)
            % obj = ANTARRAY(path, F, L, S)
            %
            % INPUT
            %   M:      square matrix of antenna elements' excitation
            %   path:   path to matrix file
            %   F:      frequency of operation [MHz]
            %   L:      length of dipole element [mm]
            %   S:      inter-element spacing [fraction of wavelength]
            % OUTPUT
            %   obj:    ANTARRAY object
            %
            
            if nargin == 0 || isempty(M)
                M = zeros(64);
            elseif isa(M, 'char')
                if strcmp(M(end-3:end), '.mat')
                    M = M(1:end-4);
                end;
                if ~strcmp(M(end-3:end), '.dat')
                    M = [M '.dat'];
                end;
                if ~exist(M, 'file')
                    M(end-2) = 'm';
                    if ~exist(M, 'file')
                        error('MyErr:FileNotFound', ['file not found: ' M]);
                    else
                        M = load(M);
                        fl = fieldnames(M);
                        M = M.(fl{1});
                    end;
                else
                    M = dlmread(M);
                end;
            end;
                
            if size(M,1) ~= size(M,2)
                error('The matrix containing the excitations must be square');
            end;
            
            if mod(obj.opt_pts, 2) == 0
                error('Invalid number of discretization steps (should be odd)');
            end;
            
            if nargin < 2
                f = 60500;
                l = 2.1826855509101;
                s = 1;
            end;
            if isempty(f)
                f = 60500;
            end;
            if isempty(l)
                l = 2.1826855509101;
            end;
            if isempty(s)
                s = 1;
            end;
            
            
            obj.M = M;
            obj.freq = f*10^6;
            obj.el_len = l/1000;
            obj.spacing = obj.c0/obj.freq*s;
            
            obj.name = '';
            obj.comments = '';
            obj.opt_win = [0,0];   
            obj.max_XY = 0;
            obj.max_YZ = 0;
            obj.min_XY = 0;
            obj.min_YZ = 0;
            obj.max_E_strength = 0;
            obj.min_E_strength = 0;
            obj.normalized = 0;
            obj.pwr = 10*10^-3;
            obj.weight_ang = pi/18; % 10°
            
            obj.dir = zeros(length(M));
            obj.dir(M~=0) = 1;
            obj.dir_str{1} = 'Unknown';
            obj.norm_freq = obj.freq;
        end
               
        %% Function to set the optimization window dimension
        function obj = setOptWin(obj, side, dist)
           % INPUT
           %    obj:    AntArray object
           %    side:   side length of the optimization window [mm]
           %    dist:   distance to the array plane [mm]
           
           obj.opt_win = [side/1000, dist/1000];
        end
        
        %% Function to set the output file name
        function obj = setName(obj, name)
            %SETNAME set the name for file saving
            %
            % obj = SETNAME(obj, name)
            %
            % INPUT
            %   obj:    AntArray object
            %   name:   desired output file name
            % OUTPUT
            %   obj:    AntArray object
            
            if name(1) ~= '_'
                name = ['_' name];
            end;
            obj.name = name;
        end
        
        %% Function to set the maximal scale value
        function obj = setMax(obj, pattern, val)
            %SETMAX set the max value for plots
            %
            % obj = SETMAX(obj, patrn, val)
            %   possible patrn values:
            %       XY  for XY-plane
            %       YZ  for YZ-plane
            %       E   for the E-strength
            %
            % INPUT
            %   obj:    AntArray object
            %   patrn:  name of the plane where the value should apply
            %   val:    maximal value [dB]
            % OUTPUT
            %   obj:    ANTARRAY object
            %
            % See also SETMIN GENPATTERN E_STRENGTH
            
            if strcmp(pattern, 'XY')
                obj.max_XY = val;
            elseif strcmp(pattern, 'YZ')
                obj.max_YZ = val;
            elseif strcmp(pattern, 'E')
                obj.max_E_strength = val;
            else
                error('Unhandled pattern plane');
            end;
        end
        
        %% Function to set the minimal scale value
        function obj = setMin(obj, pattern, val)
            %SETMIN set the min value for plots
            %
            % obj = SETMIN(obj, patrn, val)
            %   possible patrn values:
            %       XY  for XY-plane
            %       YZ  for YZ-plane
            %       E   for the E-strength
            %
            % INPUT
            %   obj:    AntArray object
            %   patrn:  name of the plane where the value should apply
            %   val:    minimal value [dB]
            % OUTPUT
            %   obj:    AntArray object
            %
            % See also SETMAX GENPATTERN E_STRENGTH
            
            if strcmp(pattern, 'XY')
                obj.min_XY = val;
            elseif strcmp(pattern, 'YZ')
                obj.min_YZ = val;
            elseif strcmp(pattern, 'E')
                obj.min_E_strength = val;
            else
                error('Unhandled pattern plane');
            end;
        end
        
        %% Function to set the aperture angle used for weighting
        function obj = setWeightAngle(obj, val)
            %SETWEIGHTANGLE set the weight angle
            %
            % obj = SETWEIGHTANGLE(obj, val)
            %
            % INPUT
            %   obj:    AntArray object
            %   val:    angle [rad]
            % OUTPUT
            %   obj:    AntArray object
            %
            % See also WEIGHT PLOTWEIGHT

            if isempty(val)
               return;
            end;

            if val > pi/2
               val = pi/2;
            elseif val < 0
               val = 0;
            end;

            obj.weight_ang = val;
        end
        
        %% Function to set comments to be printed on the elements' plot
        function obj = setComments(obj, comm)
            %SETCOMMENTS set the comments string
            %
            % The comment string will be printed on the elements' plot
            %
            % obj = SETCOMMENTS(obj, comm)
            %
            % INPUT
            %   obj:    AntArray object
            %   comm:   comments to be displayed
            % OUTPUT
            %   obj:    AntArray object
            % 
            % See also PLOTANTARRAY
            
            obj.comments = comm;
        end
        
        %% Function to set the frequency used for normalization
        function obj = setNormFreq(obj, n_freq)
            %SETNORMFREQ set the frequency used at normalization
            %
            % This frequency will be used for computations during
            % normalization of the input power to the specified level.
            %
            % obj = SETNORMFREQ(obj, n_freq)
            %
            % INPUT
            %   obj:    AntArray object
            %   n_freq: frequency used at normalization [MHz]
            % OUTPUT
            %   obj:    AntArray object
            % 
            % See also SETNORMPWR

            obj.norm_freq = n_freq*10^6;
            obj.normalized = 0;
        end
        
        %% Function to set the power used for normalization
        function obj = setNormPwr(obj, n_pwr)
            %SETNORMPWR set the input power for normalization
            %
            % Before generating an electric field distribution, the
            % amplitude of the antenna elements will be scaled as to
            % correspond to the specified normalization power.
            %
            % obj = SETNORMPWR(obj, n_pwr)
            %
            % INPUT
            %   obj:    AntArray object
            %   n_pwr:  desired input power [W]
            % OUTPUT
            %   obj:    AntArray object
            % 
            % See also SETNORMFREQ
            
            obj.pwr = n_pwr;
            obj.normalized = 0;
        end
        
        %% Turn waitbars on/off
        function obj = waitbars(obj, status)
            %WAITBARS turn waitbars on or off
            %
            % obj = WAITBARS(obj, status)
            %
            % INPUT
            %   obj:    AntArray object
            %   status: (optional) ON (1) or OFF (0) [default=1]
            % OUTPUT
            %   obj:    AntArray object

            if nargin < 2 || isempty(status)
                status = 1;
            else
                status = status > 0;
            end;
            obj.dispwait = status;
        end
        
        %% Function to add focused antenna pattern to the array
        function obj = adaptArray(obj, M, x, y, z)
            %ADAPTARRAY add elements to the antenna array
            %
            % The elements' excitation will be adapted as to focus to a
            % specific point of space.
            %
            % obj = ADAPTARRAY(obj, M, x, y, z)
            %
            % INPUT
            %   obj:    AntArray object
            %   M:      matrix of elements that need to be steered
            %   x,y,z:  desired focus point [mm]
            % OUTPUT
            %   obj:    AntArray object
            %
            % See also RSTARRAY ADAPTAMP
            
            if nargin < 3
                x = 100000;
                y = 0;
                z = 0;
            end;
            
            if size(M,1) < size(obj.M,1) || size(M,2) < size(obj.M, 2)
                error('Matrix sizes do not match');
            end;
            
            dir_index = max(max(obj.dir))+1;
            if dir_index == 1
                dir_index = 2;
            end;
            obj.dir_str{dir_index} = ['(' mat2str(x) ', ' mat2str(y) ', ' ...
                mat2str(z) ')'];
            obj.dir(M~=0) = dir_index;
            
            if ~isempty(M(M(M ~= 0)~=1))
                pos = find(M(M ~= 0)~=1);
                obj.M(pos) = M(pos);
                return
            end;
            
            
            lambda = obj.c0/obj.freq;
            s = obj.spacing;

            x = x/1000;
            y = y/1000;
            z = z/1000;

            rho_g = sqrt(x^2+y^2+z^2);
            rho = rho_g;

            L = length(obj.M)*s;
            gamma = rho_g/2/L^2*lambda;
            disp(['Gamma: ' mat2str(gamma)]);
            
            tmp_M = M;
            
            % Excite elements appropriately for focusing
            if mod(size(M,1), 2) == 0
                turnover = size(M,1)/2;

                for i=1:size(M,1)
                    z_el = (turnover-i)*s + s/2;
                    for j=1:size(M,2)
                        y_el = (j-turnover)*s - s/2;

                        if M(i,j) == 1
                            phi = 2*pi/lambda * (sqrt(x^2+(y_el-y)^2+(z_el-z)^2) - rho);
                            tmp_M(i,j) = exp(-1j*phi);
                        end;
                    end;
                end;
            else
                turnover = (size(M,1)-1)/2;

                for i=1:size(M,1)
                    z_el = (turnover-i+1)*s;
                    for j=1:size(M,2)
                        y_el = (j-turnover-1)*s;

                        if M(i,j) == 1
                            phi = 2*pi/lambda * (sqrt(x^2+(y_el-y)^2+(z_el-z)^2) - rho);
                            tmp_M(i,j) = exp(-1j*phi);
                        end;
                    end;
                end;
            end;
            
            % Assign final values
            obj.M = tmp_M;
            obj.normalized = 0;
            
        end
        
        %% Function to adapt element amplitude according to a profile
        function obj = adaptAmp(obj, norm_x, amps, els, mode)
            %ADAPTAMP adapt the elements' amplitude
            %
            % The elements' amplitude will be adapted according to the
            % specified profile.
            % The profile is defined by means of two vectors, the first
            % corresponding to the normalized positions from the centre of
            % the array, the latter corresponding to the normalized
            % amplitude factor.
            % There are two different modes of operation. The profile can
            % be applied according to the elements position from the array
            % centre, or to its effective distance to the centre. With the
            % first mode, all elements on the edge will have the same
            % factor applied. With the second mode, the corner elements are
            % located further from the centre than other edge elements and
            % will thus have a different factor applied.
            %
            % obj = ADAPTAMP(obj, norm_x, amps, els, mode)
            %
            % INPUT
            %   obj:    AntArray object   
            %   norm_x: vector of normalized x-positions
            %   amps:   vector of corresponding normalized amplitudes
            %   els:    (optional) Matrix containing the elements to be
            %           adapted
            %   mode:   method to calculate the distance to the array
            %           centre: element position (P) or effective distance
            %           (D)
            % OUTPUT
            %   obj:    AntArray object
            %
            % See also ADAPTARRAY RSTARRAY
            
            if nargin < 4 || isempty(els)
                els = obj.M;
            elseif size(els,1) ~= size(els,2) || size(els,1) ~= size(obj.M, 1)
                error('Matrix dimensions do not agree');
            end;
            
            if nargin < 5
                mode = 'P';
            elseif mode ~= 'P' && mode ~= 'D'
                error('Unhandled mode');
            end;
            
            if length(norm_x) ~= length(amps)
                error('norm_x and amps dimensions do not agree');
            elseif norm_x(1) ~= 0 || norm_x(end) ~= 1
                error('norm_x extreme values should be 0 and 1');
            else
                for i=1:length(norm_x)-1
                    if norm_x(i+1) <= norm_x(i)
                        error('norm_x must be strictly increasing');
                    end;
                end;
            end;
            
            [y, x] = find(els ~= 0);
            turnover = size(obj.M,1)/2 + 0.5;
            
            if mode == 'D'
                for i=1:length(x)
                    xx = abs(turnover - x(i))*obj.spacing;
                    yy = abs(turnover - y(i))*obj.spacing;
                    
                    els(y(i), x(i)) = sqrt(xx^2 + yy^2);
                end;
            elseif mode == 'P'
                for i=1:length(x)
                    d = max(abs(turnover-x(i)), abs(turnover-y(i)))*obj.spacing;
                    els(y(i), x(i)) = d;
                end;
            end;
            
            norm_x = norm_x.*max(max(els));
            
            for i=1:length(x)
                namp = interp1(norm_x, amps, els(y(i), x(i)));
                
                curramp = abs(obj.M(y(i),x(i)));
                obj.M(y(i), x(i)) = obj.M(y(i), x(i))/curramp*namp;
            end;
        end
        
        %% Function to reset the antenna array
        function obj = rstArray(obj, M)
            %RSTARRAY clear the elements' matrix
            %
            % obj = RSTARRAY(obj, M)
            %
            % If no matrix argument is given, a 1-element array will be
            % created
            %
            % INPUT
            %   obj:    AntArray object
            %   M:      (optional) replacement matrix
            % OUTPUT
            %   obj:    AntArray object
            %
            % See also ADAPTARRAY ADAPTAMP
            
            if nargin < 2
                M = 0;
            end;
            
            obj.M = M;
            obj.dir = M;
            obj.dir(M~=0) = 1;
            obj.dir_str = cell(1,1);
            obj.dir_str{1} = 'Unknown';
            obj.normalized = 0;
        end
        
        %% Function to compute and plot the field pattern
        function [obj, ptrn] = genPattern(obj, d, L, mode, ss, theta)
            %GENPATTERN generate electric field distribution
            %
            % Compute and plot the electric field pattern of the antenna
            % array on a plane. The plane can be parallel to the array or
            % perpendicular to it and passing through the propagation axis.
            % Before each generation, the amplitude of the antenna elements
            % will be scaled as to correspond to the specified input power.
            % The patterns are computed at the frequency given at
            % contruction of the AntArray object.
            %
            % [obj ptrn] = GENPATTERN(obj, d, L, mode, ss, theta)
            %
            % The behaviour is determined by the mode parameter, which can
            % take three different values:
            % If mode = 'YZ':
            %   Plot the fields on a square surface parallel to the YZ
            %   plane
            % If mode = 'YZ-MAIN':
            %   Plot only the main beam in the YZ-plane
            % If mode = 'XY':
            %   Plot the fields in the XY-plane
            % If mode = 'THETA':
            %   Plot the fields in the plane at a given angle from the
            %   Z-axis, through the X-axis
            %
            % Moreover, once a pattern has been generated with one of the
            % previous mode, a black and white pattern can be generated by
            % using the same command with 'YZ-BW', 'XY-BW', or 'THETA-BW'
            % as mode parameter instead. The white region in the generated
            % pattern is the region where the electric field intensity is
            % greater than the minimum reception level.
            %   
            %
            % INPUT
            %   obj:    AntArray object
            %   d:      distance to the array (YZ) or distance from the
            %           array (XY and THETA) [mm]
            %   L:      side length of the plot surface [mm]
            %   mode:   YZ, XY or THETA see above
            %   ss:     step size for the plot [mm]
            %   theta:  [if mode=theta] angle wrt Z-axis [rad]
            % OUTPUT
            %   obj:    AntArray object
            %   ptrn:   generated pattern
            %
            % See also SET_E_MIN SETNORMPWR ANTARRAY
            
            plot_val = obj.plotres;
            
            if nargout == 2
                obj.plotres = 0;
            end;
            mode = upper(mode);
            
            L = L/1000;
            d = d/1000;

            if nargin < 5 || isempty(ss)
                ss = obj.spacing;
            else
                ss = ss/1000;
                step_num = round(L/ss);
                if mod(step_num, 2) ~= 0
                    step_num = step_num + 1;
                end;
                ss = L/step_num;
            end;
            
            if obj.normalized == 0
                fprintf('Normalizing input power...');
                obj = obj.normalize();
                fprintf('\tdone\n');
            end;
            
            disp(['Generating ' mode ' pattern...']);

            if strcmp(mode, 'YZ')
                fprintf(['\tStep size: ' mat2str(ss*1000) 'mm\n']);
                fprintf(['\tArray size: ' mat2str(length(obj.M)*obj.spacing*1000) 'mm\n']);
                
                if length(d) > 1
                    warning('Only the pattern at the latest distance will be returned');
                end;

                for i=1:length(d)
                    fprintf(['\tGenerating pattern ' mat2str(i) ...
                        ' of ' mat2str(length(d)) '...']);
                    
                    ptrn = E_YZ(obj, d(i), L, ss);
                    fprintf('\tdone\n');
                end;
            elseif strcmp(mode, 'XY')
                fprintf(['\tStep size: ' mat2str(ss*1000) 'mm\n']);
                fprintf(['\tArray size: ' mat2str(length(obj.M)*obj.spacing*1000) 'mm\n']);

                ptrn = E_XY(obj, d, L, ss);
                fprintf('\tdone\n');
            elseif strcmp(mode, 'YZ-MAIN')
                fprintf(['\tStep size: ' mat2str(ss*1000) 'mm\n']);
                fprintf(['\tArray size: ' mat2str(length(obj.M)*obj.spacing*1000) 'mm\n']);
                
                 if length(d) > 1
                    warning('Only the pattern at the latest distance will be returned');
                end;

                for i=1:length(d)
                    fprintf(['\tGenerating pattern ' mat2str(i) ...
                        ' of ' mat2str(length(d)) '...']);
                    
                    ptrn = E_YZ_main(obj, d(i), L, ss);
                    fprintf('\tdone\n');
                end;
            elseif strcmp(mode, 'YZ-BW')
                for i=1:length(d)
                    fprintf(['\tGenerating BW pattern ' mat2str(i) ...
                        ' of ' mat2str(length(d)) '...']);
                    ptrn = E_BW(obj, d(i));
                    fprintf('\tdone\n');
                end;
            elseif strcmp(mode, 'XY-BW')
                ptrn = E_BW(obj);
                fprintf('\tdone\n');
            elseif strcmp(mode, 'THETA')
                if isempty(theta)
                    theta = pi/4;
                end;
                fprintf(['\tStep size: ' mat2str(ss*1000) 'mm\n']);
                fprintf(['\tArray size: ' mat2str(length(obj.M)*obj.spacing*1000) 'mm\n']);

                ptrn = E_theta(obj, d, L, ss, theta);
                fprintf('\tdone\n');
            elseif strcmp(mode, 'THETA-BW')
                if isempty(theta)
                    theta = pi/4;
                end;
                
                ptrn = E_BW(obj, [], theta);
                fprintf('\tdone\n');
            else
                error('Unhandled mode');
            end;
            
            obj.plotres = plot_val;
        end
        
        %% Function to compute the E-field strength along a line
        function obj = E_strength(obj, L, y, z, smp)
            %E_STRENGTH compute the electric field along a line
            %
            % The electric field intensity is computed along a line
            % parallel to the X-axis.
            % Before each generation, the amplitude of the antenna elements
            % will be scaled as to correspond to the specified input power.
            % The fields are computed at the frequency given at
            % contruction of the AntArray object.
            %
            % obj = E_STRENGTH(obj, L, y, z, smp)
            %
            % INPUT
            %   obj:    AntArray object
            %   L:      maximal distance to the array plane [mm]
            %   y:      (optional) Y-position of the line
            %   z:      (optional) Z-position of the line
            %   smp:    (optional) number of samples [default 100]
            % OUTPUT
            %   obj:    AntArray object
            %
            % See also SETNORMPWR ANTARRAY
            
            if nargin < 5
                smp = 100;
            end;
            if nargin < 4 || isempty(z)
                z = 0;
            end;
            if nargin < 3 || isempty(y)
                y = 0;
            end;
            if nargin < 2
                error('Not enought input arguments');
            end;
            
            if obj.normalized == 0
                fprintf('Normalizing input power...');
                obj = obj.normalize();
                fprintf('\tdone\n');
            end;
            
            fprintf('Computing E strength profile...');
            
            y = y/1000;
            z = z/1000;
            
            % Create waitbar
            if obj.dispwait
                progress = waitbar(0, 'Computations in progress...', ...
                    'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
                setappdata(progress, 'canceling', 0);
                posProg = get(progress, 'Position');
                uicontrol('Parent', progress, 'Style', 'pushbutton', ...
                    'Position', [posProg(3)*0.12, posProg(4)*0.22, 80, 23], ...
                    'String', 'Terminate', ...
                    'Callback', 'setappdata(gcbf, ''terminating'', 1)', ...
                    'Visible', 'on');
                setappdata(progress, 'terminating', 0);

                waitbar(0);
            end;

            absc = logspace(2, ceil(log10(L)), smp);
            oord = zeros(1, length(absc));
            
            % Partitions for parallel computing
            iterations = 1;
            s_samples = smp;
            adapted_samples = 0;
            if smp > 50
                divs = factor2(round(smp));
                if length(divs) <= 2
                    divs = factor2(round(smp+1));
                    adapted_samples = 1;
                end;
                for i=1:round(length(divs)/2)
                    if divs(end-i+1) > iterations && divs(i) < 25
                        iterations = divs(i);
                    end;
                end;
                if adapted_samples
                    s_samples = (smp+1)/iterations;
                else
                    s_samples = smp/iterations;
                end;
            end;
            
            % Computations
            for k=1:iterations
                if adapted_samples && k == iterations
                    cut_oord = zeros(1, s_samples-1);
                    cut_absc = absc(1, end-length(cut_oord):end);
                else
                    cut_oord = zeros(1, s_samples);
                    cut_absc = absc(1, (k-1)*s_samples+1:k*s_samples);
                end;
                if iterations == 1
                    for i=1:length(cut_oord)
                        E = zeros(1,3);
                        [E(1), E(2), E(3)] = E_array(obj, cut_absc(i)/1000, y, z);
                        cut_oord(i) = 20*log10(sqrt(sum(abs(E(:)).^2)));

                        waitbar(i/smp);
                        
                        if getappdata(progress, 'canceling')
                            close all;
                            delete(progress);
                            return;
                        elseif getappdata(progress, 'terminating')
                            close all;
                            delete(progress);
                            throw(MException('MyERR:Terminated', ...
                            'Program terminated by user'));
                        end;
                    end;
                else
                    parfor i=1:length(cut_oord)
                        E = zeros(1,3);
                        [E(1), E(2), E(3)] = E_array(obj, cut_absc(i)/1000, y, z);
                        cut_oord(i) = 20*log10(sqrt(sum(abs(E(:)).^2)));
                    end;
                    
                    if getappdata(progress, 'canceling')
                        close all;
                        delete(progress);
                        return;
                    elseif getappdata(progress, 'terminating')
                        close all;
                        delete(progress);
                        throw(MException('MyERR:Terminated', ...
                            'Program terminated by user'));
                    end;
                    
                    waitbar(k/iterations);
                end;
                
                if adapted_samples && k == iterations
                    oord(1, end-length(cut_oord):end) = cut_oord(:);
                else
                    oord(1, (k-1)*s_samples+1:k*s_samples) = cut_oord(:);
                end;
            end;
        
            % Plot
            waitbar(1, progress, 'Generating plots...');
            fig = figure();
            semilogx(absc, oord, 'LineWidth', 2);
            xlim([absc(1) absc(end)]);
            hold on;
            
            view_axis = axis;
            if obj.max_E_strength ~= 0
                max_val = obj.max_E_strength-1;
            else
                max_val = max(oord);
            end;
            if obj.min_E_strength ~= 0
                min_val = obj.min_E_strength+1;
            else
                min_val = view_axis(3);
            end;
            ylim([min_val-5+mod(min_val,5) max_val+5-mod(max_val,5)]);
            
            % Find fitline
            pos = length(oord);
            prev = oord(pos);
            curr = oord(pos-1);
            pos = pos-2;
            while curr >= prev && pos > 0
                prev = curr;
                curr = oord(pos);
                pos = pos - 1;
            end;
            pos = pos + 1;
            if pos > .05*length(oord)
                pmax = pos;
                prev = oord(pos);
                curr = oord(pos-1);
                pos = pmax-2;
                while pos > 0 && curr <= prev
                    prev = curr;
                    curr = oord(pos);
                    pos = pos - 1;
                end;
                if pos ~= 0
                    bias = mean(oord(1:pos));
                    parea = axis;
                    if bias > oord(pmax)
                        plot([parea(1) parea(2)], [bias bias], '--', ...
                            'Color', [1, 0.35, 0], 'LineWidth', 1.2);
                        vert = absc(pmax);
                    else
                        plot([parea(1) parea(2)], [bias bias], '--k', ...
                            'LineWidth', 1);
                        pos = length(oord);
                        prev = oord(pos);
                        curr = oord(pos-1);
                        pos = pos-2;
                        while prev <= bias && curr < bias
                            prev = curr;
                            curr = oord(pos);
                            pos = pos - 1;
                        end;
                        dx = log10(absc(pos+2)) - log10(absc(pos+1));
                        dy = prev - curr;
                        Dx = dx/dy*(bias-prev);
                        vert = 10^(log10(absc(pos+2))+Dx);
                    end;
                    
                    plot([vert vert], [parea(3) parea(4)], '--k', ...
                        'LineWidth', 1);
                    
                    text(absc(round(.96*length(absc))), ...
                        bias - (parea(4)-parea(3))*.035, ...
                        [num2str(bias,3) '\,dB\,V/m'], ...
                        'HorizontalAlignment', 'right', ...
                        'Interpreter', 'latex', 'FontSize', 18);
                    text(absc(round(.96*length(absc))), ...
                        bias - (parea(4)-parea(3))*.08, ...
                        [int2str(round(vert)) '\,mm'], ...
                        'HorizontalAlignment', 'right', ...
                        'Interpreter', 'latex', 'FontSize', 18);
                end;
            end;
            
            % Adapt title and axes
            title(['\textbf{Field strength along (' mat2str(z) ', ' ...
                mat2str(y) ')}'], 'Interpreter', 'latex', 'FontSize', 24);
            xlabel('Distance to the array centre [mm]', ...
                'Interpreter', 'latex', 'FontSize', 22);
            ylabel('Field strength [dB\,V/m]', ...
                'Interpreter', 'latex', 'FontSize', 22);
            set(get(fig, 'CurrentAxes'), 'FontSize', 16);
            hold off;

            print_plots(fig, ['E_strength' obj.name]);

            delete(progress);

            close all;
            fprintf('\tdone\n');
        end
        
        %% Weighting function
        function w = weight(obj, mode, d)
            %WEIGHT compute the weight of the realized pattern
            %
            % The weight consists of the sum of all cells with an electric
            % field intensity above the minimal level, that are located
            % inside the weight angle.
            %
            % Patterns must have been generated with the GENPATTERN
            % function before computing their weight.
            %
            % w = WEIGHT(obj, mode, d)
            %
            %
            % INPUT
            %   obj:    AntArray object
            %   mode:   XY, YZ or THETA, see explanation in GENPATTERN
            %   d:      (optional) distance of the YZ pattern [mm]
            %           OR theta angle [rad]
            % OUTPUT
            %   w:      corresponding weight
            %
            % See also GENPATTERN SET_MIN_E SETWEIGHTANGLE PLOTWEIGHT
            
            w = 0;
            
            if strcmp(mode, 'XY')
                savname = ['pattern_XY' obj.name];
                mode = 0;
            elseif strcmp(mode, 'YZ')
                if isempty(d)
                    error('Parameter d is required for this mode');
                end;
                savname = ['pattern' obj.name '_' mat2str(d/1000)];
                mode = 1;
            elseif strcmp(mode, 'theta')
                if isempty(d)
                    error('Parameter d is required for this mode');
                end;
                savname = ['pattern_theta_' mat2str(d, 3) obj.name];
                mode = 0;
            else
                error('Invalid mode parameter');
            end;
            
            filename = AntArray.openFile(savname);
            if filename == 0
                return;
            end;
            
            A = dlmread(filename);
            
            x_dev = A(1,2:end-1).*1000; % Convert to mm
            A = A(2:end-1, 2:end-1);
            
            ss = abs(x_dev(2) - x_dev(1));
            
            if mode==0
                elW = zeros(size(A,1), floor(size(A,2)/2));
                
                for i=1:size(elW,2)
                    elW(round(i/tan(obj.weight_ang/2)):end,i) = 1;
                end;
                
                hwidth = size(obj.M,1)/2*obj.spacing*1000; % Convert to mm
                hwidth = floor(hwidth/ss);
                
                elW(1:end,1:hwidth-1) = 1;
                
                elW = [elW(:,end:-1:1) ones(size(elW,1),1) elW];
            elseif mode==1
                elW = zeros(ceil(size(A,1)/2));
                
                r1 = floor(size(obj.M,1)*obj.spacing*1000/ss);
                intersect = r1/tan(obj.weight_ang/2);
                
                if d <= intersect
                    r = r1;
                else
                    r = d*tan(obj.weight_ang/2);
                end;
                
                elW(1, 1:min(end, round(r/ss)+1)) = 1;
                
                for i=1:min(size(elW,1), floor(r/ss))
                    maxval = round(sqrt(r^2-(i*ss)^2)/ss)+1;
                    if maxval < 1e-6
                        continue;
                    end;
                    if maxval > size(elW,1)
                        elW(i+1, 1:end) = 1;
                    else
                        elW(i+1, 1:maxval) = 1;
                    end;
                end;
                elW(elW'==1) = 1;

                elW = [elW(end:-1:2,end:-1:2) elW(end:-1:2,:);
                        elW(:,end:-1:2) elW];
            end;
            
            weightM = zeros(size(A,1), size(A,2));
            weightM(A>obj.min_E) = 1;
            weightM(elW ~= 1) = 0;
            w = sum(sum(weightM));
        end
        
        %% Plot weighting function
        function plotWeight(obj, mode, d)
            %PLOTWEIGHT plot the pattern realized by the weighting angle
            %
            % [ ] =PLOT WEIGHT(obj, mode, d)
            %
            % INPUT
            %   obj:    AntArray object
            %   mode:   XY, YZ or THETA, see explanation in GENPATTERN
            %   d:      (optional) distance of the YZ pattern [mm]
            %           OR theta angle [rad]
            %
            % See also GENPATTERN SET_MIN_E SETWEIGHTANGLE WEIGHT
            
            if nargin < 3
                d = [];
            end;
            
            if strcmp(mode, 'XY')
                savname = ['pattern_XY' obj.name];
                mode = 0;
            elseif strcmp(mode, 'YZ')
                if isempty(d)
                    error('Parameter d is required for this mode');
                end;
                savname = ['pattern' obj.name '_' mat2str(d/1000)];
                mode = 1;
            elseif strcmp(mode, 'theta')
                if isempty(d)
                    error('Parameter d is required for this mode');
                end;
                savname = ['pattern_theta_' mat2str(d, 3) obj.name];
                mode = 0;
            else
                error('Invalid mode parameter');
            end;
            
            filename = AntArray.openFile(savname);
            if filename == 0
                return;
            end;
            
            A = dlmread(filename);
            
            x_dev = A(1,2:end-1).*1000;
            y_dev = A(2:end-1,1)'.*1000; % Convert to mm
            A = A(2:end-1, 2:end-1);
            
            ss = abs(x_dev(2) - x_dev(1));
            
            if mode==0
                elW = zeros(size(A,1), floor(size(A,2)/2));
                
                for i=1:size(elW,2)
                    elW(round(i/tan(obj.weight_ang/2)):end,i) = 1;
                end;
                
                hwidth = size(obj.M,1)/2*obj.spacing*1000; % Convert to mm
                hwidth = floor(hwidth/ss);
                
                elW(1:end,1:hwidth-1) = 1;
                
                elW = [elW(:,end:-1:1) ones(size(elW,1),1) elW];
            elseif mode==1
                elW = zeros(ceil(size(A,1)/2));
                
                r = d*tan(obj.weight_ang/2);
                
                elW(1, 1:min(end, round(r/ss)+1)) = 1;
                
                for i=1:min(size(elW,1), floor(r/ss))
                    maxval = round(sqrt(r^2-(i*ss)^2)/ss)+1;
                    if maxval < 1e-6
                        continue;
                    end;
                    if maxval > size(elW,1)
                        elW(i+1, 1:end) = 1;
                    else
                        elW(i+1, 1:maxval) = 1;
                    end;
                end;
                elW(elW'==1) = 1;

                elW = [elW(end:-1:2,end:-1:2) elW(end:-1:2,:);
                        elW(:,end:-1:2) elW];
            end;
            
            % ========================================================================
            % Plot
            x_dev = x_dev./1000; % Reconvert to m
            y_dev = y_dev./1000;
            plotdata = [elW zeros(size(elW,1),1); zeros(1,size(elW,2)+1)];
            % Generate axes
            if x_dev(end) < 1/100
                fact_x = 1000;
            elseif x_dev(end) < 1
                fact_x = 100;
            else
                fact_x = 1;
            end;
            if y_dev(end) < 1/100
                fact_y = 1000;
            elseif y_dev(end) < 1
                fact_y = 100;
            else
                fact_y = 1;
            end;

            range_x = x_dev.*fact_x;
            range_y = y_dev.*fact_y;
            ss_x = x_dev(2) - x_dev(1);
            ss_y = y_dev(2) - y_dev(1);
            [absc, oord] = meshgrid([range_x range_x(end)+ss_x*fact_x], ...
                [range_y range_y(end)+ss_y*fact_y]);  % Larger to be able to plot evth

            % Plot field
            figure(1);
            surf(absc, oord, plotdata, 'EdgeColor', 'none', ...
                'LineStyle', 'none');
            view(2);
            hold on;
            colormap gray;
            colorbar('eastoutside');

            % Adapt color map
            caxis([0 1]);

            % Adapt ticks
            if mod(2*x_dev(end)*fact_x,4) == 0
                tick_fact_x = 4;
            elseif mod(2*x_dev(end)*fact_x,6) == 0
                tick_fact_x = 6;
            else
                tick_fact_x = 2;
            end;
            if mod(2*y_dev(end)*fact_y,4) == 0
                tick_fact_y = 4;
            elseif mod(2*y_dev(end)*fact_y,6) == 0
                tick_fact_y = 6;
            else
                tick_fact_y = 2;
            end;

            spe_ticks_x = zeros(tick_fact_x+1,1);
            for ii=1:length(spe_ticks_x)
                spe_ticks_x(ii) = range_x(1) + (ii-1)*2*x_dev(end)*fact_x/tick_fact_x;
            end;
            spe_ticks_pos_x = spe_ticks_x+ss_x*fact_x/2;
            spe_ticks_y = zeros(tick_fact_y+1,1);
            for ii=1:length(spe_ticks_y)
                spe_ticks_y(ii) = range_y(1) + (ii-1)*2*y_dev(end)*fact_y/tick_fact_y;
            end;
            spe_ticks_pos_y = spe_ticks_y+ss_y*fact_y/2;

            xlim([range_x(1), range_x(end)+ss_x*fact_x]);
            ylim([range_y(1), range_y(end)+ss_y*fact_y]);
            set(gca, 'XTick', spe_ticks_pos_x, ...
                'XTickLabel', spe_ticks_x);
            set(gca, 'YTick', spe_ticks_pos_y, ...
                'YTickLabel', spe_ticks_y);

            % Set labels and title
            if mode == 1
                title(['\textbf{Desired electric field at ', mat2str(d), 'm}'], ...
                    'Interpreter', 'latex', 'FontSize', 24);
            elseif mode == 0 && ~isempty(d)
                title(['\textbf{Desired electric field at $\theta=\pi\cdot' ...
                rats(d/pi) '$}'], ...
                'Interpreter', 'latex', 'FontSize', 24);
            else
                title('\textbf{Desired electric field in the XY plane}', ...
                    'Interpreter', 'latex', 'FontSize', 24);
            end;
            
            switch fact_x
                case 1000
                    unit_x = 'mm';
                case 100
                    unit_x = 'cm';
                otherwise
                    unit_x = 'm';
            end;
            switch fact_y
                case 1000
                    unit_y = 'mm';
                case 100
                    unit_y = 'cm';
                otherwise
                    unit_y = 'm';
            end;
            if mode == 1
                xlabel(['${\rm y}_{\rm pos}$ [' unit_x ']'], ...
                    'Interpreter', 'latex', 'FontSize', 22);
                ylabel(['${\rm z}_{\rm pos}$ [' unit_y ']'], ...
                    'Interpreter', 'latex', 'FontSize', 22);
            elseif mode == 0 && ~isempty(d)
                xlabel(['Distance to the array centre [' unit_x ']'], 'Interpreter', ...
                    'latex', 'FontSize', 22);
                ylabel(['${\rm x}_{\rm pos}$ [' unit_y ']'], 'Interpreter', ...
                    'latex', 'FontSize', 22);
            else
                xlabel(['${\rm y}_{\rm pos}$ [' unit_x ']'], ...
                    'Interpreter', 'latex', 'FontSize', 22);
                ylabel(['${\rm x}_{\rm pos}$ [' unit_y ']'], ...
                    'Interpreter', 'latex', 'FontSize', 22);
            end;
            set(gca, 'FontSize', 16);
            
            print_plots(gcf, [savname '_weight']);
            
            close all;
            
        end
        
        %% Function to compute the directivity
        function obj = directivity(obj)
            %DIRECTIVITY compute and plot the directivity
            %
            % This function uses the standard formula for computing the
            % directivity, which takes more time than the alternative
            % version
            %
            % obj = DIRECTIVITY(obj)
            %
            % INPUT
            %   obj:    AntArray object
            % OUTPUT
            %   obj:    AntArray object
            %
            % See also DIRECTIVITY_ALT

            r = 50000;
            fun = @(theta, phi) ...
                obj.dir_den(r, theta, phi);
            
            den = integral2(fun, 0, pi, 0, 2*pi)/4/pi/r^2;
            
            D = 0:180/5;
            theta0 = pi/3;
            for phi_i=1:length(D)
                phi0 = (D(phi_i)-90)*5*pi/180;
                
                D(phi_i) = 2*obj.E_tot(r, theta0, phi0)/den;
            end;
            
            disp(['mean = ' num2str(mean(D),4)]);
            
            % ========================================================================
            % Plot
            figure(1);
            plot(-90:5:90, D, 'LineWidth', 2);
            hold on
            
            xlim([-90 90]);
            
            title(['\textbf{Directivity at $\theta=' ...
                num2str(theta0*180/pi) '$}'], ...
                'Interpreter', 'latex', 'FontSize', 24);
            xlabel('$\phi$ [deg]', 'Interpreter', 'latex', 'FontSize', 22);
            ylabel('Directivity []', 'Interpreter', 'latex', 'FontSize', 22);
            set(gca, 'FontSize', 16);
            
            hold off
        end
        
        %% Function to compute the directivity
        function obj = directivity_alt(obj)
            %DIRECTIVITY_ALT compute and plot the directivity
            %
            % This function uses an altenative method to compute
            % directivity, which requires computation of the input power of
            % the array. Generally, it is faster than the standard
            % expression using integrals.
            %
            % obj = DIRECTIVITY_ALT(obj)
            %
            % INPUT
            %   obj:    AntArray object
            % OUTPUT
            %   obj:    AntArray object
            %
            % See also DIRECTIVITY SETNORMPWR
            
            r = 90000;            
            den = obj.inpower()/4/pi/r^2;
            
            D = 0:180/5;
            theta0 = pi/3;
            for phi_i=1:length(D)
                phi0 = (D(phi_i)-90)*5*pi/180;
                
                D(phi_i) = 2*obj.E_tot(r, theta0, phi0)/den/obj.Z0;
            end;
            
            disp(['mean = ' num2str(mean(D))]);
            
            % ========================================================================
            % Plot
            figure(1);
            plot(-90:5:90, D, 'LineWidth', 2);
            hold on
            
            xlim([-90 90]);
            
            title(['\textbf{Directivity at $\theta=' ...
                num2str(theta0*180/pi) '$}'], ...
                'Interpreter', 'latex', 'FontSize', 24);
            xlabel('$\phi$ [deg]', 'Interpreter', 'latex', 'FontSize', 22);
            ylabel('Directivity []', 'Interpreter', 'latex', 'FontSize', 22);
            set(gca, 'FontSize', 16);
            
            hold off
        end
        
        %% Function to show the antenna repartition
        function plotAntArray(obj, cont)
            %PLOTANTARRAY plot the antenna arrangement
            %
            % Plot the antenna arrangement of matrix M. Their point of
            % focus is also indicated, as well as the antenna size. If a
            % comment was defined, it will be printed too.
            %
            % [ ] = PLOTANTARRAY(obj, cont)
            %
            % INPUT
            %   obj:    AntArray object
            %   cont:   (optional) whether to plot the contour of the
            %           antenna array [default = 0]
            %
            % See also SETCOMMENTS ADAPTARRAY RSTARRAY
            
            if nargin < 2 || cont < 1
                cont = 0;
            else
                cont = 1;
            end;
            
            colors = {'m', 'c', 'r', 'g', 'b', 'k'};
            markers = {'o', '+', 'x', 's', 'd'};
            
            maxval = max(max(obj.dir));
            if maxval == 0
                error('No antenna element was initiated');
            end;
            
            fig = figure();
            hold on;
            for i=1:maxval
                if ~isempty(find(obj.dir,i))
                    [rl, cl] = find(obj.dir == i);
                    rl = length(obj.M) - rl + 1;
                    colorind = mod(i, length(colors));
                    markerind = mod(ceil(i/length(markers)),length(markers));
                    plot(cl, rl, 'LineStyle', 'none', ...
                        'DisplayName', [obj.dir_str{i} ' - ' mat2str(length(rl))], ...
                        'MarkerFaceColor', colors{colorind}, ...
                        'Marker', markers{markerind}, 'MarkerSize', 5);
                end;
            end;
            
            xlim([-1 length(obj.M)+1]);
            ylim([-1 length(obj.M)+1]);
            
            text(0, -length(obj.M)*0.03, ...
                [mat2str(length(obj.M)) 'x' mat2str(length(obj.M))], ...
                'Interpreter', 'latex', 'FontSize', 20);
            
            L = legend('Location', 'eastoutside');
            set(L, 'Interpreter', 'latex', 'FontSize', 20);
            
            if cont == 1
                cc = length(obj.M)+1;
                plot([0 cc], [0 0], '--k', 'Linewidth', 1);
                plot([0 cc], [cc cc], '--k', 'Linewidth', 1);
                plot([0 0], [0 cc], '--k', 'Linewidth', 1);
                plot([cc cc], [0 cc], '--k', 'Linewidth', 1);
            end;
            
            title('\textbf{Array arrangement}', ...
                'Interpreter', 'latex', 'FontSize', 24);
            axis off;
            
            if ~strcmp(obj.comments, '')
                drawnow
                text(length(obj.M)*1.04, -length(obj.M)*0.02, ...
                    obj.comments, ...
                    'Interpreter', 'latex', 'FontSize', 18);
            end

            savname = ['elements' obj.name];
            print_plots(fig, savname);
            export_dat(round(abs(obj.M)), savname, 1);

            close all
        end
    end
    methods (Access='private')
        %% Function to compute and plot the fields for the YZ-mode
        function ptrn = E_YZ(obj, d, L, ss)
            % INPUT
            %   obj:    AntArray object
            %   d:      Distance to the array [m]
            %   L:      Side length of the plot surface [m]
            %   ss:     Step size for the plot [m]
            
            dim = round(L/ss)+1;
            ext_dim = dim+1;
            plotdata = zeros(ext_dim); % Larger to be able to plot evth

            % Create waitbar
            if obj.dispwait
                progress = waitbar(0, 'Computations in progress...', ...
                    'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
                setappdata(progress, 'canceling', 0);
                posProg = get(progress, 'Position');
                uicontrol('Parent', progress, 'Style', 'pushbutton', ...
                    'Position', [posProg(3)*0.12, posProg(4)*0.22, 80, 23], ...
                    'String', 'Terminate', ...
                    'Callback', 'setappdata(gcbf, ''terminating'', 1)', ...
                    'Visible', 'on');
                setappdata(progress, 'terminating', 0);

                waitbar(0);
            end;
            
            % Search for symmetry
            dim_z = dim;
            dim_y = dim;

            if isempty(obj.M(obj.M ~= obj.M(end:-1:1,:)))
                dim_z = ceil(dim_z/2);
            end;
            if isempty(obj.M(obj.M ~= obj.M(:, end:-1:1)))
                dim_y = ceil(dim_y/2);
            end;
            
            % Compute
            for i=1:dim_z
                z = L/2 - (i-1)*ss;
                slice = plotdata(ext_dim-i,:);
                parfor j=1:dim_y
                    y = (j-1)*ss - L/2;
                    E = zeros(1,3);
                    [E(1), E(2), E(3)] = E_array(obj, d, y, z);
                    slice(j) = 20*log10(sqrt(sum(abs(E(:)).^2)));
                end;
                plotdata(ext_dim-i,:) = slice;
                
                if obj.dispwait
                    waitbar(i/dim_z);
                    if getappdata(progress, 'canceling')
                        close all;
                        delete(progress);
                        return;
                    elseif getappdata(progress, 'terminating')
                        close all;
                        delete(progress);
                        throw(MException('MyERR:Terminated', ...
                        'Program terminated by user'));
                    end;
                end;
            end;
            
            if dim_y ~= dim
                plotdata(:, end-1:-1:end/2) = plotdata(:, 1:end/2);
            end;
            if dim_z ~= dim
                plotdata(end/2:-1:1, :) = plotdata(end/2:end-1, :);
            end;
            
            ptrn = plotdata(1:end-1, 1:end-1);

            if ~obj.plotres
                if obj.dispwait
                    delete(progress);
                end;
                return
            end;
            if obj.dispwait
                waitbar(1, progress, 'Generating plots...');
            end;

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
            surf(absc, oord, plotdata, 'EdgeColor', 'none', ...
                'LineStyle', 'none');
            view(2);
            hold on;
            colormap jet;
            cbar_h = colorbar('eastoutside');
            
            % Plot antenna
            cmax = max(max(plotdata(:,:)));
            antsize = length(obj.M)*obj.spacing*fact;
            lpos = antsize/2;
            plot3([-lpos+ss/2 lpos+ss/2], [-lpos+ss/2 -lpos+ss/2], ...
                [cmax+5 cmax+5], '-k', 'LineWidth', 1);
            plot3([-lpos+ss/2 lpos+ss/2], [lpos+ss/2 lpos+ss/2],...
                [cmax+5 cmax+5], '-k', 'LineWidth', 1);
            plot3([-lpos+ss/2 -lpos+ss/2], [-lpos+ss/2 lpos+ss/2],...
                [cmax+5 cmax+5], '-k', 'LineWidth', 1);
            plot3([lpos+ss/2 lpos+ss/2], [-lpos+ss/2 lpos+ss/2],...
                [cmax+5 cmax+5], '-k', 'LineWidth', 1);
            plot3(ss/2, ss/2, cmax+5, '-k', 'MarkerFaceColor', 'k', ...
                'Marker', '+', 'MarkerSize', 5);

            % Adapt color map
            title(cbar_h, 'dB\,V/m', 'Interpreter', 'latex', 'FontSize', 16);
            if obj.max_YZ ~= 0
                cmax = obj.max_YZ-1;
            end;
            if obj.min_YZ ~= 0
                cmin = obj.min_YZ+1;
            else
                cmin = min(min(plotdata));
            end;
            caxis([cmin-5+mod(cmin,5) cmax+5-mod(cmax,5)]);

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
            title(['\textbf{Electric field at ', mat2str(d), '\,m}'], ...
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

            savname = ['pattern' obj.name '_' mat2str(d)];
            print_plots(gcf, savname);
            export_dat([0 absc(1,:)./fact; oord(:,1)./fact plotdata], savname);
            
            if obj.dispwait
                delete(progress);
            end;

            close all
        end
        
        %% Function to compute and plot the main fields for the YZ-mode
        function ptrn = E_YZ_main(obj, d, L, ss)
            % INPUT
            %   obj:    AntArray object
            %   d:      Distance to the array [m]
            %   L:      Side length of the plot surface [m]
            %   ss:     Step size for the plot [m]
            
            dim = round(L/ss)+1;
            ext_dim = dim+1;
            plotdata = ones(ext_dim); % Larger to be able to plot evth
            plotdata = plotdata.*(AntArray.min_E-40); % 40dB below threshold
                        
            % Search for symmetry
            % -------------------
            dim_z = dim;
            dim_y = dim;

            if isempty(obj.M(obj.M ~= obj.M(end:-1:1,:)))
                dim_z = ceil(dim_z/2);
            end;
            if isempty(obj.M(obj.M ~= obj.M(:, end:-1:1)))
                dim_y = ceil(dim_y/2);
            end;
            
            % Compute
            % -------
            counter_th = 1;
            parfor i=1:dim_z
                counter = 0;
                
                z = L/2 - (i-1)*ss;
                slice = plotdata(ext_dim-i,:);

                if dim_y ~= dim
                    for j=dim_y:-1:1
                        y = (j-1)*ss - L/2;
                        E = zeros(1,3);
                        [E(1), E(2), E(3)] = E_array(obj, d, y, z);
                        val = 20*log10(sqrt(sum(abs(E(:)).^2)));
                        slice(j) = val;
                        
                        if val < AntArray.min_E
                            counter = counter + 1;
                            if counter >= counter_th
                                break;
                            end;
                        elseif counter ~= 0
                            counter = 0;
                        end;                        
                    end;
                else
                    for j=floor(dim_y/2):-1:1
                        y = (j-1)*ss - L/2;
                        E = zeros(1,3);
                        [E(1), E(2), E(3)] = E_array(obj, d, y, z);
                        val = 20*log10(sqrt(sum(abs(E(:)).^2)));
                        slice(j) = val;
                        
                        if val < AntArray.min_E
                            counter = counter + 1;
                            if counter >= counter_th
                                break;
                            end;
                        elseif counter ~= 0
                            counter = 0;
                        end;   
                    end;
                    for j=ceil(dim_y/2):dim_y
                        y = (j-1)*ss - L/2;
                        E = zeros(1,3);
                        [E(1), E(2), E(3)] = E_array(obj, d, y, z);
                        val = 20*log10(sqrt(sum(abs(E(:)).^2)));
                        slice(j) = val;
                        
                        if val < AntArray.min_E
                            counter = counter + 1;
                            if counter >= counter_th
                                break;
                            end;
                        elseif counter ~= 0
                            counter = 0;
                        end;   
                    end;
                end;
                    
                plotdata(ext_dim-i,:) = slice;
            end;
            
            if dim_y ~= dim
                plotdata(:, end-1:-1:end/2) = plotdata(:, 1:end/2);
            end;
            if dim_z ~= dim
                plotdata(end/2:-1:1, :) = plotdata(end/2:end-1, :);
            end;
            
            ptrn = plotdata(1:end-1, 1:end-1);
            
            if ~obj.plotres
                return
            end;

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
            surf(absc, oord, plotdata, 'EdgeColor', 'none', ...
                'LineStyle', 'none');
            view(2);
            hold on;
            colormap jet;
            cbar_h = colorbar('eastoutside');
            
            % Plot antenna
            cmax = max(max(plotdata(:,:)));
            antsize = length(obj.M)*obj.spacing*fact;
            lpos = antsize/2;
            plot3([-lpos+ss/2 lpos+ss/2], [-lpos+ss/2 -lpos+ss/2], ...
                [cmax+5 cmax+5], '-k', 'LineWidth', 1);
            plot3([-lpos+ss/2 lpos+ss/2], [lpos+ss/2 lpos+ss/2],...
                [cmax+5 cmax+5], '-k', 'LineWidth', 1);
            plot3([-lpos+ss/2 -lpos+ss/2], [-lpos+ss/2 lpos+ss/2],...
                [cmax+5 cmax+5], '-k', 'LineWidth', 1);
            plot3([lpos+ss/2 lpos+ss/2], [-lpos+ss/2 lpos+ss/2],...
                [cmax+5 cmax+5], '-k', 'LineWidth', 1);
            plot3(ss/2, ss/2, cmax+5, '-k', 'MarkerFaceColor', 'k', ...
                'Marker', '+', 'MarkerSize', 5);

            % Adapt color map
            title(cbar_h, 'dB\,V/m', 'Interpreter', 'latex', 'FontSize', 16);
            if obj.max_YZ ~= 0
                cmax = obj.max_YZ-1;
            end;
            if obj.min_YZ ~= 0
                cmin = obj.min_YZ+1;
            else
                cmin = min(min(plotdata));
            end;
            caxis([cmin-5+mod(cmin,5) cmax+5-mod(cmax,5)]);

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
            title(['\textbf{Electric field at ', mat2str(d), '\,m}'], ...
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

            savname = ['pattern' obj.name '_' mat2str(d) '_main'];
            print_plots(gcf, savname);
            export_dat([0 absc(1,:)./fact; oord(:,1)./fact plotdata], savname);

            close all
        end
        
        %% Function to generate a BW plot from existing plot
        function ptrn = E_BW(obj, d, theta)
            %INPUT
            %   d:  Distance to the array [m]
            
            if nargin < 2
                d = [];
            end;
            if nargin < 3
                theta = [];
            end;
            
            if ~isempty(d)
                savname = ['pattern' obj.name '_' mat2str(d)];
            elseif ~isempty(theta)
                savname = ['pattern_theta_' mat2str(theta,3) obj.name];
            else
                savname = ['pattern_XY' obj.name];
            end;
            
            filename = AntArray.openFile(savname);
            if filename == 0
                return;
            end;
            
            A = dlmread(filename);
            
            x_dev = A(1,2:end-1);
            y_dev = A(2:end-1,1)';
            A = A(2:end, 2:end);
            
            plotdata = zeros(size(A,1), size(A,2));
            plotdata(A >= obj.min_E)=1;
            
            ptrn = plotdata(1:end-1, 1:end-1);
            
            if ~obj.plotres
                delete(progress);
                return
            end;
            
            % ========================================================================
            % Plot
            % Generate axes
            if x_dev(end) < 1/100
                fact_x = 1000;
            elseif x_dev(end) < 1
                fact_x = 100;
            else
                fact_x = 1;
            end;
            if y_dev(end) < 1/100
                fact_y = 1000;
            elseif y_dev(end) < 1
                fact_y = 100;
            else
                fact_y = 1;
            end;

            range_x = x_dev.*fact_x;
            range_y = y_dev.*fact_y;
            ss_x = x_dev(2) - x_dev(1);
            ss_y = y_dev(2) - y_dev(1);
            [absc, oord] = meshgrid([range_x range_x(end)+ss_x*fact_x], ...
                [range_y range_y(end)+ss_y*fact_y]);  % Larger to be able to plot evth

            % Plot field
            figure(1);
            surf(absc, oord, plotdata, 'EdgeColor', 'none', ...
                'LineStyle', 'none');
            view(2);
            hold on;
            colormap gray;
            colorbar('eastoutside');

            % Adapt color map
            caxis([0 1]);
            
            % If YZ, plot antenna
            if ~isempty(d)
                cmax = max(max(plotdata)) + 10;
                antsize_x = length(obj.M)*obj.spacing*fact_x;
                antsize_y = length(obj.M)*obj.spacing*fact_y;
                lpos_x = antsize_x/2;
                lpos_y = antsize_y/2;
                plot3([-lpos_x+ss_x/2 lpos_x+ss_x/2], ...
                    [-lpos_y+ss_y/2 -lpos_y+ss_y/2], ...
                    [cmax+5 cmax+5], '-r', 'LineWidth', 1);
                plot3([-lpos_x+ss_x/2 lpos_x+ss_x/2], ...
                    [lpos_y+ss_y/2 lpos_y+ss_y/2],...
                    [cmax+5 cmax+5], '-r', 'LineWidth', 1);
                plot3([-lpos_x+ss_x/2 -lpos_x+ss_x/2], ...
                    [-lpos_y+ss_y/2 lpos_y+ss_y/2], ...
                    [cmax+5 cmax+5], '-r', 'LineWidth', 1);
                plot3([lpos_x+ss_x/2 lpos_x+ss_x/2], ...
                    [-lpos_y+ss_y/2 lpos_y+ss_y/2], ...
                    [cmax+5 cmax+5], '-r', 'LineWidth', 1);
                plot3(ss_x/2, ss_y/2, cmax+5, '-r', ...
                    'MarkerFaceColor', 'r', 'Marker', '+', 'MarkerSize', 5);
            end;

            % Adapt ticks
            if mod(2*x_dev(end)*fact_x,4) == 0
                tick_fact_x = 4;
            elseif mod(2*x_dev(end)*fact_x,6) == 0
                tick_fact_x = 6;
            else
                tick_fact_x = 2;
            end;
            if mod(2*y_dev(end)*fact_y,4) == 0
                tick_fact_y = 4;
            elseif mod(2*y_dev(end)*fact_y,6) == 0
                tick_fact_y = 6;
            else
                tick_fact_y = 2;
            end;

            spe_ticks_x = zeros(tick_fact_x+1,1);
            for ii=1:length(spe_ticks_x)
                spe_ticks_x(ii) = range_x(1) + (ii-1)*2*x_dev(end)*fact_x/tick_fact_x;
            end;
            spe_ticks_pos_x = spe_ticks_x+ss_x*fact_x/2;
            spe_ticks_y = zeros(tick_fact_y+1,1);
            for ii=1:length(spe_ticks_y)
                spe_ticks_y(ii) = range_y(1) + (ii-1)*2*y_dev(end)*fact_y/tick_fact_y;
            end;
            spe_ticks_pos_y = spe_ticks_y+ss_y*fact_y/2;

            xlim([range_x(1), range_x(end)+ss_x*fact_x]);
            ylim([range_y(1), range_y(end)+ss_y*fact_y]);
            set(gca, 'XTick', spe_ticks_pos_x, ...
                'XTickLabel', spe_ticks_x);
            set(gca, 'YTick', spe_ticks_pos_y, ...
                'YTickLabel', spe_ticks_y);

            % Set labels and title
            if ~isempty(d)
                title(['\textbf{Electric field at ', mat2str(d), 'm}'], ...
                    'Interpreter', 'latex', 'FontSize', 24);
            elseif ~isempty(theta)
                title(['\textbf{Electric field at $\theta=\pi\cdot' ...
                rats(theta/pi) '$}'], ...
                'Interpreter', 'latex', 'FontSize', 24);
            else
                title('\textbf{Electric field in the XY plane}', ...
                    'Interpreter', 'latex', 'FontSize', 24);
            end;
            
            switch fact_x
                case 1000
                    unit_x = 'mm';
                case 100
                    unit_x = 'cm';
                otherwise
                    unit_x = 'm';
            end;
            switch fact_y
                case 1000
                    unit_y = 'mm';
                case 100
                    unit_y = 'cm';
                otherwise
                    unit_y = 'm';
            end;
            if ~isempty(d)
                xlabel(['${\rm y}_{\rm pos}$ [' unit_x ']'], ...
                    'Interpreter', 'latex', 'FontSize', 22);
                ylabel(['${\rm z}_{\rm pos}$ [' unit_y ']'], ...
                    'Interpreter', 'latex', 'FontSize', 22);
            elseif ~isempty(theta)
                xlabel(['Distance to the array centre [' unit_x ']'], 'Interpreter', ...
                    'latex', 'FontSize', 22);
                ylabel(['${\rm x}_{\rm pos}$ [' unit_y ']'], 'Interpreter', ...
                    'latex', 'FontSize', 22);
            else
                xlabel(['${\rm y}_{\rm pos}$ [' unit_x ']'], ...
                    'Interpreter', 'latex', 'FontSize', 22);
                ylabel(['${\rm x}_{\rm pos}$ [' unit_y ']'], ...
                    'Interpreter', 'latex', 'FontSize', 22);
            end;
            set(gca, 'FontSize', 16);
            
            print_plots(gcf, [savname '_BW']);
            
            close all;
        end
        
        %% Function to compute and plot the fields for the XY-mode
        function ptrn = E_XY(obj, d, L, ss)
            % INPUT
            %   obj:    AntArray object
            %   d:      Distance to the array (maximal plot distance) [m]
            %   L:      Side length of the plot surface [m]
            %   ss:     Step size for the plot [m]
            
            dim1 = round(L/ss)+1;
            dim2 = round(d/ss)+1;
            plotdata = zeros(dim2+1, dim1+1); % Larger to be able to plot evth

            % Create waitbar
            progress = waitbar(0, 'Computations in progress...', ...
                'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
            setappdata(progress, 'canceling', 0);
            posProg = get(progress, 'Position');
            uicontrol('Parent', progress, 'Style', 'pushbutton', ...
                'Position', [posProg(3)*0.12, posProg(4)*0.22, 80, 23], ...
                'String', 'Terminate', ...
                'Callback', 'setappdata(gcbf, ''terminating'', 1)', ...
                'Visible', 'on');
            setappdata(progress, 'terminating', 0);
            
            waitbar(0);
            
            % Symmetry
            dim_y = dim1;
            
            if isempty(obj.M(obj.M ~= obj.M(:, end:-1:1)))
                dim_y = ceil(dim1/2);
            end;
            
            % Compute
            for i=1:dim2
                x = (i-1)*ss;
                slice = plotdata(i,:);
                parfor j=1:dim_y
                    y = (j-1)*ss - L/2;
                    E = zeros(1,3);
                    [E(1), E(2), E(3)] = E_array(obj, x, y, 0);
                    slice(j) = 20*log10(sqrt(sum(abs(E(:)).^2)));
                end;
                plotdata(i,:) = slice;
                
                waitbar(i/dim2);
                if getappdata(progress, 'canceling')
                    close all;
                    delete(progress);
                    return;
                elseif getappdata(progress, 'terminating')
                    close all;
                    delete(progress);
                    throw(MException('MyERR:Terminated', ...
                    'Program terminated by user'));
                end;
            end;
            
            if dim_y ~= dim1
                plotdata(:, end/2:end-1) = plotdata(:, end/2:-1:1);
            end;
            
            ptrn = plotdata(1:end-1, 1:end-1);

            if ~obj.plotres
                delete(progress);
                return
            end;
            
            waitbar(1, progress, 'Generating plots...');

            % ==============================================================
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
            if obj.max_XY == 0
                ln_nb = ceil(10e-2/ss);
                max_val = max(max(plotdata(ln_nb:end,:)));
            else
                max_val = obj.max_XY-1;
            end;
            if obj.min_XY ~= 0
                min_val = obj.min_XY+1;
            else
                min_val = min(min(plotdata(:,:)));
            end;
            caxis([min_val-5+mod(min_val,5) max_val+5-mod(max_val,5)]);

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
            title('\textbf{Electric field in the XY plane}', ...
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

            savname = ['pattern_XY' obj.name];
            print_plots(gcf, savname);
            export_dat([0 absc(1,:)./fact_L; oord(:,1)./fact_d plotdata], savname);
            
            delete(progress);
            
            close all;
        end
        
        %% Function to compute and plot the fields for the theta-mode
        function ptrn = E_theta(obj, d, L, ss, theta)
            % INPUT
            %   obj:    AntArray object
            %   d:      Distance to the array (maximal plot distance) [m]
            %   L:      Side length of the plot surface [m]
            %   ss:     Step size for the plot [m]
            %   theta:  Angle of the plot surface wrt z-axis [radians]
            
            if theta > pi/2
                theta = mod(theta, pi/2);
            elseif theta < -pi/2
                theta = -mod(-theta, pi/2);
            end;
            
            dim1 = round(L/ss)+1;
            dim2 = round(d/ss)+1;
            plotdata = zeros(dim2+1, dim1+1);    % Larger to be able to plot evth

            % Create waitbar
            progress = waitbar(0, 'Computations in progress...', ...
                'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
            setappdata(progress, 'canceling', 0);
            posProg = get(progress, 'Position');
            uicontrol('Parent', progress, 'Style', 'pushbutton', ...
                'Position', [posProg(3)*0.12, posProg(4)*0.22, 80, 23], ...
                'String', 'Terminate', ...
                'Callback', 'setappdata(gcbf, ''terminating'', 1)', ...
                'Visible', 'on');
            setappdata(progress, 'terminating', 0);
            
            waitbar(0);

            % Compute
            for i=1:dim2
                x = (i-1)*ss;
                slice = plotdata(i,:);
                parfor j=1:dim1
                    y_plot = (j-1)*ss - L/2;
                    y = y_plot*sin(theta);
                    z = y_plot*cos(theta);
                    E = zeros(1,3);
                    [E(1), E(2), E(3)] = E_array(obj, x, y, z);
                    slice(j) = 20*log10(sqrt(sum(abs(E(:)).^2)));
                end;
                plotdata(i,:) = slice;
                
                waitbar(i/dim2);
                if getappdata(progress, 'canceling')
                    close all;
                    delete(progress);
                    return;
                elseif getappdata(progress, 'terminating')
                    close all;
                    delete(progress);
                    throw(MException('MyERR:Terminated', ...
                    'Program terminated by user'));
                end;
            end;
            
            ptrn = plotdata(1:end-1, 1:end-1);
            
            if ~obj.plotres
                delete(progress);
                return
            end;

            waitbar(1, progress, 'Generating plots...');

            % ========================================================================
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
            if obj.max_XY == 0
                ln_nb = ceil(10e-2/ss);
                max_val = max(max(plotdata(ln_nb:end,:)));
            else
                max_val = obj.max_XY-1;
            end;
            if obj.min_XY ~= 0
                min_val = obj.min_XY+1;
            else
                min_val = min(min(plotdata(:,:)));
            end;
            caxis([min_val-5+mod(min_val,5) max_val+5-mod(max_val,5)]);

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
            title(['\textbf{Electric field at $\theta=\pi\cdot' ...
                rats(theta/pi) '$}'], ...
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
            xlabel(['Distance to the array centre [' unit_L ']'], 'Interpreter', ...
                'latex', 'FontSize', 22);
            ylabel(['${\rm x}_{\rm pos}$ [' unit_d ']'], 'Interpreter', ...
                'latex', 'FontSize', 22);
            set(gca, 'FontSize', 16);

            savname = ['pattern_theta_' mat2str(theta, 3) obj.name];
            print_plots(gcf, savname);
            export_dat([0 absc(1,:)./fact_L; oord(:,1)./fact_d plotdata], savname);
            
            delete(progress);
            
            close all;
        end
        
        %% Function to compute the electric field components of the array
        function [E_x, E_y, E_z] = E_array(obj, x, y, z, M_opt)
            % Compute the electric field components of the array in one
            % single point
            %
            % INPUT
            %    obj:    AntArray object
            %    x,y,z:  Coordinates of the point [m]
            %    M_opt:  (optional) optimization antenna array
            % OUTPUT
            %    E_{k}:  Electric field components in the {k}-direction

            if nargin == 5
                M_in = M_opt;
            else
                M_in = obj.M;
            end;
            
            s_in = obj.spacing;
            f_in = obj.freq;
            l_in = obj.el_len;

            E_x = 0;
            E_y = 0;
            E_z = 0;

            if mod(size(M_in,1), 2) == 0
                turnover = size(M_in,1)/2 + 0.5;
            else
                turnover = (size(M_in,1)+1)/2;
            end;
            
            P = find(M_in ~= 0);
            
            for p=1:numel(P)
                i = mod(P(p), size(M_in,1));
                if i==0
                   i = size(M_in,1);
                end
                j = ceil(P(p)/size(M_in, 1));
            
                z_el = (turnover-i)*s_in;
                zz = z - z_el;
                
                y_el = (j-turnover)*s_in;
                yy = y - y_el;
                
                E = zeros(1,3);
                
                [E(1), E(2), E(3)] = ...
                    AntArray.E_dipole(l_in, M_in(i,j), f_in, x, yy, zz);

                E_x = E_x + E(1);
                E_y = E_y + E(2);
                E_z = E_z + E(3);
            end;           
        end
        
        %% Function to compute E-field strength
        function E_strength = dir_den(obj, r, theta, phi)
           E = zeros(1,3);
           E_strength = zeros(size(theta,1), size(theta,2));
           for i=1:numel(theta)
                [E(1), E(2), E(3)] = obj.E_array(r*sin(theta(i))*cos(phi(i)), ...
                    r*sin(theta(i))*sin(phi(i)), r*cos(theta(i)));
                E_strength(i) = sum(abs(E(:)).^2)*sin(theta(i))*r^2;
           end;
        end
        
        %% Function to compute E-field strength
        function E_strength = E_tot(obj, r, theta, phi)
           E = zeros(1,3);
           [E(1), E(2), E(3)] = obj.E_array(r*sin(theta)*cos(phi), ...
                    r*sin(theta)*sin(phi), r*cos(theta));
           E_strength = sum(abs(E(:)).^2);
        end
        
        %% Function to compute the array input power
        function input_pwr = inpower_old(obj)
            [elcol, elrow] = find(obj.M ~= 0);
            
            if obj.dispwait
                progress = waitbar(0, 'Power computations in progress...',...
                    'CreateCancelBtn',...
                    'setappdata(gcbf,''canceling'',1)');
                setappdata(progress, 'canceling', 0);
            end;
            
            input_pwr = zeros(size(obj.M,1), size(obj.M,1), numel(elcol));
            
            if mod(size(obj.M,1), 2) == 0
                turnover = size(obj.M,1)/2 + 0.5;
            else
                turnover = (size(obj.M,1)+1)/2;
            end;
            
            k0 = 2*pi/obj.c0*obj.norm_freq;
            ss = obj.spacing;
            Zc = obj.Z0;
            ellen = obj.el_len;
            
            M_tmp = obj.M(obj.M ~= 0);
            
            for i=1:numel(obj.M)
                if obj.M(i) == 0
                    continue;
                end;
                row = mod(i, size(obj.M,1));
                col = ceil(i/size(obj.M, 1));
                if row == 0
                    row = size(obj.M,1);
                end;
                
                I1 = obj.M(i);

                z_el1 = (turnover-row)*ss;
                y_el1 = (col-turnover)*ss;
                parfor j=1:numel(elcol)
                    z_el2 = (turnover-elcol(j))*ss;
                    y_el2 = (elrow(j)-turnover)*ss;
                    
                    z_dist = (z_el1-z_el2)^2;
                    y_dist = (y_el1-y_el2)^2;
                    dist = sqrt(z_dist + y_dist);

                    if dist < 1e-8
                        input_pwr(row,col,j) = Zc*k0^2*ellen^2/6/pi;
                        input_pwr(row,col,j) = input_pwr(row,col,j)*I1*conj(M_tmp(j));
                    else
                        ang = atan((row-elcol(j))/(elrow(j)-col));
                        psi = k0*dist;

                        input_pwr(row,col,j) = Zc*k0^2*ellen^2/4/pi * ...
                            (2*(sin(psi)-psi*cos(psi))/psi^3 * (cos(ang))^2 + ...
                            ((psi^2-1)*sin(psi)+psi*cos(psi))/psi^3 * (sin(ang))^2);
                        input_pwr(row,col,j) = input_pwr(row,col,j)*I1*conj(M_tmp(j));
                    end;                        
                end;
                if obj.dispwait
                    waitbar(i/numel(obj.M));
                    if getappdata(progress, 'canceling')
                        close all;
                        delete(progress);
                        input_pwr = -1;
                        return;
                    end;
                end;
            end;
            input_pwr = abs(sum(sum(sum(input_pwr(:,:,:)))));
            
            if obj.dispwait
                delete(progress);
            end;
        end
        
        %% Function to compute the array input power
        function input_pwr = inpower(obj)
            [elcol, elrow, elval] = find(obj.M);
            
            if islogical(elval)
                elval = ones(size(elval,1), size(elval,2));
            end;
            
            if obj.dispwait
                progress = waitbar(0, 'Power computations in progress...',...
                    'CreateCancelBtn', ...
                    'setappdata(gcbf,''canceling'',1)');
                setappdata(progress, 'canceling', 0);
            end;
            
            input_pwr = zeros(size(obj.M,1), size(obj.M,1), numel(elcol));
            
            k0 = 2*pi/obj.c0*obj.norm_freq;
            ss = obj.spacing;
            Zc = obj.Z0;
            ellen = obj.el_len;
            
            % create dist & ang matrix
            D = zeros(size(input_pwr, 1), size(input_pwr, 2));
            Ang = D;
            for i=1:size(D,1)
                for j=1:i
                    D(i,j) = sqrt(ss^2*((i-1)^2+(j-1)^2));
                    Ang(i,j) = atan((j-1)/(i-1));
                end;
            end;
            Ang(1,1:end) = 0;
            
            D2 = D';
            D2(D ~= 0) = 0;
            D = D + D2;
            Ang2 = 2*pi - Ang';
            Ang2(Ang ~= 0) = 0;
            Ang(Ang ~= 0) = Ang(Ang ~= 0) + 3*pi/2;
            Ang = Ang + Ang2;
            Ang(1,1) = 0;
            Ang(2:end,1) = 3*pi/2;
            Ang2 = Ang';
            
            for i=1:numel(elcol)
                y1 = elcol(i);
                x1 = elrow(i);
                tmp_D = [D(y1:-1:1, x1:-1:1) D(y1:-1:1, 2:end-x1+1);
                    D(2:end-y1+1, x1:-1:1) D(2:end-y1+1, 2:end-x1+1)];
                P1 = Ang(y1:-1:2, x1:-1:2) - pi;
                P2 = Ang2(y1:-1:2, 1:end-x1+1) - 3*pi/2;
                P3 = Ang2(1:end-y1+1, x1:-1:2) - pi/2;
                tmp_A = [P1 P2;
                    P3 Ang(1:end-y1+1, 1:end-x1+1)];
                
                for j=1:numel(elcol)
                    y2 = elcol(j);
                    x2 = elrow(j);
                    if tmp_D(y2, x2) < 1e-8
                        val = Zc*k0^2*ellen^2/6/pi;
                        val = val * elval(i)*conj(elval(j));
                    else
                        psi = k0*tmp_D(y2, x2);
                        ang = tmp_A(y2, x2);
                        
                        val = Zc*k0^2*ellen^2/4/pi * ...
                            (2*(sin(psi)-psi*cos(psi))/psi^3 * (cos(ang))^2 + ...
                            ((psi^2-1)*sin(psi)+psi*cos(psi))/psi^3 * (sin(ang))^2);
                        val = val * elval(i)*conj(elval(j));
                    end;
                    input_pwr(y1, x1, j) = val;
                end;
                if obj.dispwait
                    waitbar(i/numel(elcol));
                    if getappdata(progress, 'canceling')
                        close all;
                        delete(progress);
                        input_pwr = -1;
                        return;
                    end;
                end;
            end;
            
            input_pwr = abs(sum(sum(sum(input_pwr(:,:,:)))));
            
            if obj.dispwait
                delete(progress);
            end;
        end
        
        %% Function to normalize the currents to match desired input power
        function obj = normalize(obj)
            input_pwr = inpower(obj)/2; % Divide by 2 as only radiated in half-space
            
            if input_pwr < 0
                warning('MyWARN:norm_power', 'Power was not normalized');
            else
                fact = sqrt(input_pwr/obj.pwr);
                obj.M = obj.M(:,:)./fact;
                obj.normalized = 1;
            end;
        end
        
    end
    methods (Static, Access = 'private')
        %% Function to compute the electric field component of a dipole
        function [E_x, E_y, E_z] = E_dipole(l, I, f, x, y, z)
            % Compute the electric field components generated by a dipole
            % antenna at a certain position in space
            %
            % INPUT
            %   l:      length of the dipole [m]
            %   I:      constant current through the dipole [A]
            %   f:      frequency of operation [Hz]
            %   x,y,z:  position wrt the dipole centre [m]
            
            Zc = AntArray.Z0;    % Characteristic impedance of free space

            r = sqrt(x^2+y^2+z^2);  % Convert to polar coordinates
            if r < 1e-16
                r = 1e-16;
            end;
            theta = acos(z/r);
            phi = atan2(y,x);       % Or atan(y/x)

            lambda = AntArray.c0/f; % Wavelength
            k = 2*pi/lambda;        % Wave number

            E_r = Zc/2/pi * I*l*cos(theta)/(r^2) * (1+1/(1j*k*r)) * exp(-1j*k*r);
            E_theta = 1j*Zc*k/4/pi * I*l*sin(theta)/r * (1+1/(1j*k*r)-1/(k*r)^2) * exp(-1j*k*r);
            % H_phi = 1j*k/4/pi * I*l*sin(theta)/r * (1+1/(1j*k*r)) * exp(-1j*k*r);

            E_x = E_r*sin(theta)*cos(phi) + E_theta*cos(theta)*cos(phi);
            E_y = E_r*sin(theta)*sin(phi) + E_theta*cos(theta)*sin(phi);
            E_z = E_r*cos(theta) - E_theta*sin(theta); 
        end
        
        %% Open file
        function filename = openFile(savname)
            attempts = 30;
            
            date_n = datenum(date);
            date_v = datevec(date_n);
            dirpath = cell(1,3);
            for i=1:3
                dirpath{i} = mat2str(date_v(i));
                if numel(dirpath{i}) < 2
                    dirpath{i} = ['0' dirpath{i}];
                end;
            end;
            dirpath = sprintf('%s', dirpath{:});
            dirpath = [dirpath '/dat'];
            
            filename = [dirpath '/' savname '.dat'];
            
            while (~exist(dirpath, 'dir') || ~exist(filename, 'file')) && attempts > 0
                date_n = addtodate(date_n, -1, 'day');
                attempts = attempts-1;
                
                date_v = datevec(date_n);
                dirpath = cell(1,3);
                for i=1:3
                    dirpath{i} = mat2str(date_v(i));
                    if numel(dirpath{i}) < 2
                        dirpath{i} = ['0' dirpath{i}];
                    end;
                end;
                dirpath = sprintf('%s', dirpath{:});
                dirpath = [dirpath '/dat'];
                
                filename = [dirpath '/' savname '.dat'];
            end;
            
            if attempts == 0
                warning('MyWARN:file_not_found', ...
                    ['File ' savname ' not found, skipped']);
                filename = 0;
            end; 
        end
    end
    methods (Static, Access = 'public')
        %% Get field value at a given position
        function val = getVal(x, y, file, folder)
            %GETVAL get the field value at a given position
            %
            % Extract the value of the electric field at a given position
            % from the given 'dat' file
            %
            % val = GETVAL(x, y, file, folder)
            % 
            % INPUT
            %    x,y:    coordinates of interest [mm]
            %    file:   name of the source file [.dat extension]
            %    folder: (optional) Name of the source folder
            % OUTPUT
            %   val:    electric field intensity [dB V/m]
            %
            % See also GENPATTERN
           
            if nargin < 4
                date_v = datevec(date);
                dirpath = cell(1,3);
                for i=1:3
                    dirpath{i} = mat2str(date_v(i));
                end;
                folder = sprintf('%s', dirpath{:});
            end;
            if folder(end) == '/'
                folder = folder(1:end-1);
            end;
            if strcmp(folder(end-3:end), '/dat')
                folder = folder(1:end-4);
            end;
            if ~strcmp(file(end-3:end), '.dat')
                file = [file '.dat'];
            end;
            
            if ~exist([folder '/dat/' file], 'file')
                error('The specified file does not exist');
            end;
            
            A = dlmread([folder '/dat/' file]);
            A = A(1:end-1,1:end-1);
            
            [x_dev, x_ind] = min(abs(A(1,2:end)-x));
            [y_dev, y_ind] = min(abs(A(2:end,1)-y));
         
            val = A(x_ind+1, y_ind+1);
            
            if x_dev > abs(A(1,3)-A(1,2)) || y_dev > abs(A(3,1)-A(2,1))
                warning('MyWARN:out_of_plot_area', ...
                    'The specified point might be out of plotted area');
            end;
            
        end
        
        %% Quantize a matrix
        function M = quantize(M, quant, lvl)
            %QUANTIZE quantize a matrix
            %
            % Quantize a matrix, using a smaller matrix. When the number of
            % active elements inside the smaller matrix is higher than a
            % given level, all the elements will be turned on. Otherwise
            % they will all be turned off.
            %
            % M = quantize(M, quant, lvl)
            %            
            % INPUT:
            %   M:      matrix to be quantized
            %   quant:  dimension of quantization
            %   lvl:    quantization level
            % OUPUT
            %   M:      quantized matrix
            
            quant = abs(round(quant));
            if nargin < 3
                lvl = 0;
            else
                lvl = abs(round(lvl));
            end;
            
            if size(M, 1) ~= size(M, 2)
                error 'M should be a square matrix';
            elseif mod(length(M), 2) ~= 0
                error 'M should be a multiple of quant';
            end;
            
            M(M < 10^-6) = 0;
            M(M > 10^-6) = 1;
            
            if quant > 1
                half_dim = round(size(M,1)/quant);

                for i=1:half_dim
                    y = (i-1)*quant + 1;
                    for j=1:half_dim
                        x = (j-1)*quant + 1;
                        sq = M(y:y+quant-1, x:x+quant-1);
                        if sum(sum(sq)) < lvl
                            sq = zeros(quant);
                        else
                            sq = ones(quant);
                        end;
                        M(y:y+quant-1, x:x+quant-1) = sq(:,:);
                    end;
                end;
            end;
        end;
    end
end