%
% Copyright 2015-2016, Antoine JUCKLER
%

classdef AntArray
    properties (Constant, Access='public')
        c0 = 299792458; % Speed of light in free space
        Z0 = 119.9169832*pi; % Impedance of free space
        opt_pts = 51;	% Number of discretization steps in one dimension
        min_E = -1.3;   % Minimal electric field for reception [dB]
    end
    properties (GetAccess='public', SetAccess='private')
        M;              % matrix of elements' excitation
        freq;           % frequency of operation
        norm_freq;      % frequency used for normalization
        el_len;         % dipoles' length
        spacing;        % inter-element spacing
        opt_win;        % side length and depth of optimization window
        name;           % name for the plot files
        max_XY;         % max dB scale value for XY-patterns
        max_YZ;         % max dB scale value for YZ-patterns
        min_XY;         % min dB scale value for XY-patterns
        min_YZ;         % min dB scale value for YZ-patterns
        max_E_strength; % max dB scale value for E-strength plot
        min_E_strength; % min dB scale value for E-strength plot
        dir;            % matrix of element's groups
        dir_str;        % cell array of element's groups
        comments;       % string to be printed on elements' plot
        
        normalized;     % Has the array been normalized?
        pwr;            % Input power [W]
        
        weight_ang;     % Aperture angle for weighting [rad]
    end
    methods
        %% Constructor
        function obj = AntArray(M, f, l, s)
            % INPUT
            %   M:      square matrix of antenna elements' excitation
            %   f:      frequency of operation [MHz]
            %   l:      length of dipole element [mm]
            %   s:      inter-element spacing [fraction of wavelength]
            
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
            obj.pwr = 10e-3;
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
            % INPUT
            %   obj:    AntArray object
            %   name:   desired output file name
            
            if name(1) ~= '_'
                name = ['_' name];
            end;
            obj.name = name;
        end
        
        %% Function to set the maximal scale value
        function obj = setMax(obj, pattern, val)
            % INPUT
            %   obj:        AntArray object
            %   pattern:    name of the plane where the value should apply
            %   val:        maximal value on the dB-scale
            
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
            % INPUT
            %   obj:        AntArray object
            %   pattern:    name of the plane where the value should apply
            %   val:        minimal value on the dB-scale
            
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
            % INPUT
            %   obj:        AntArray object
            %   val:        Weighting angle [rad]
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
        function obj = setComments(obj, comments)
            % INPUT
            %   obj:        AntArray object
            %   comments:   Comments to be displayed on the plot
            obj.comments = comments;
        end
        
        %% Function to set the frequency used for normalization
        function obj = setNormFreq(obj, freq)
            % INPUT
            %   obj:        AntArray object
            %   freq:       Frequency used at normalization [MHz]
            obj.norm_freq = freq;
            obj.normalized = 0;
        end;
        
        %% Function to add focused antenna pattern to the array
        function obj = adaptArray(obj, M, x, y, z)
            % Adapt the excitations of the input matrix components to focus
            % the beam at the given location.
            % If the optimization window of the AntArray object was
            % previously defined, the Levenberg-Marquardt algorithm will be
            % used to optimize the beam focussing
            %
            % INPUT
            %   obj:    AntArray object
            %   M:      matrix of elements that need to be steered
            %   x,y,z:  desired focus point [mm]
            
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
            
            % Start values for optimization
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
            
            % Optimize
            pos = find(M);
            if obj.opt_win(1) ~= 0
                tmp_M = optArray(obj, pos, tmp_M, x, y, z);
            end;
            
            % Assign final values
            obj.M(pos) = tmp_M(pos);
            obj.normalized = 0;
            
        end
        
        %% Function to adapt element amplitude according to a profile
        function obj = adaptAmp(obj, norm_x, amps, els, mode)
        	% Adapt the amplitude of selected antenna elements according to
        	% the inputted profile
            %
            % INPUT
            %   obj:    AntArray object   
            %   norm_x: Vector of normalized x-positions
            %   amps:   Vector of corresponding normalized amplitudes
            %   els:    (optional) Matrix containing the elements to be
            %           adapted
            %   mode:   Method to calculate the distance to the array
            %           centre: element position (P) or effective distance
            %           (D)
            
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
            % Reset the matrix containing the element's excitations, if no
            % matrix is specified, a 1-element array will be created
            %
            % INPUT
            %   obj:    AntArray object
            %   M:      (optional) replacement matrix
            
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
        function obj = genPattern(obj, d, L, mode, ss, theta)
            % Compute and plot the field pattern of the antenna array
            % Behaviour depends on the "mode" parameter:
            % If mode = 'YZ':
            %    Plot the fields on a square surface parallel to the YZ-plane
            % If mode = 'XY':
            %    Plot the fields in the XY-plane and the field strength
            %    along the X-axis
            % If mode = 'theta':
            %    Plot the fields in the plane at a given angle from the
            %    Z-axis
            %
            % INPUT
            %    obj:    AntArray object
            %    d:      Distance to the array [mm]
            %    L:      Side length of the plot surface [mm]
            %    mode:   YZ, XY or theta see above
            %    ss:     Step size for the plot [mm]
            %    theta:  [if mode=theta] angle wrt Z-axis [radians]
           
            L = L/1000;
            d = d/1000;

            if nargin < 5
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
                fprintf('Computing array impedance...');
                obj = obj.normalize();
                fprintf('\tdone\n');
            end;
            
            disp(['Generating ' mode ' pattern...']);

            if strcmp(mode, 'YZ')
                fprintf(['\tStep size: ' mat2str(ss*1000) 'mm\n']);
                fprintf(['\tArray size: ' mat2str(length(obj.M)*obj.spacing*1000) 'mm\n']);

                for i=1:length(d)
                    fprintf(['\tGenerating pattern ' mat2str(i) ...
                        ' of ' mat2str(length(d)) '...']);
                    E_YZ(obj, d(i), L, ss);
                    fprintf('\tdone\n');
                end;
            elseif strcmp(mode, 'XY')
                fprintf(['\tStep size: ' mat2str(ss*1000) 'mm\n']);
                fprintf(['\tArray size: ' mat2str(length(obj.M)*obj.spacing*1000) 'mm\n']);

                E_XY(obj, d, L, ss);
                fprintf('\tdone\n');
            elseif strcmp(mode, 'YZ-BW')
                for i=1:length(d)
                    fprintf(['\tGenerating BW pattern ' mat2str(i) ...
                        ' of ' mat2str(length(d)) '...']);
                    E_BW(obj, d(i));
                    fprintf('\tdone\n');
                end;
            elseif strcmp(mode, 'XY-BW')
                E_BW(obj);
                fprintf('\tdone\n');
            elseif strcmp(mode, 'theta')
                if isempty(theta)
                    theta = pi/4;
                end;
                fprintf(['\tStep size: ' mat2str(ss*1000) 'mm\n']);
                fprintf(['\tArray size: ' mat2str(length(obj.M)*obj.spacing*1000) 'mm\n']);

                E_theta(obj, d, L, ss, theta);
                fprintf('\tdone\n');
            elseif strcmp(mode, 'theta-BW')
                if isempty(theta)
                    theta = pi/4;
                end;
                
                E_BW(obj, [], theta);
                fprintf('\tdone\n');
            else
                error('Unhandled mode');
            end;
        end
        
        %% Function to compute the E-field strength along a line
        function obj = E_strength(obj, L, y, z, samples)
            % Compute the electric field of the array along a line
            % starting at the array plane
            %
            % INPUT
            %   obj:    AntArray object
            %   L:      Maximal distance to the array plane [mm]
            %   y:      (optional) Y-position of the line
            %   z:      (optional) Z-position of the line
            %   samples:(optional) nb of samples [default 100]
            
            if nargin < 5
                samples = 100;
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
                fprintf('Computing array impedance...');
                obj = obj.normalize();
                fprintf('\tdone\n');
            end;
            
            fprintf('Computing E strength profile...');
            
            y = y/1000;
            z = z/1000;
            
            progress = waitbar(0, 'Computations in progress...', ...
                'CreateCancelBtn', ...
                'setappdata(gcbf,''canceling'',1)');
            setappdata(progress, 'canceling', 0);
            
            absc = logspace(2, ceil(log10(L)), samples);
            oord = zeros(1, length(absc));

            for i=1:samples
                E = zeros(1,3);
                [E(1), E(2), E(3)] = E_array(obj, absc(i)/1000, y, z);
                oord(i) = 20*log10(sqrt(sum(abs(E(:)).^2)));
                if getappdata(progress, 'canceling')
                    close all;
                    delete(progress);
                    return;
                end;
                waitbar(i/samples);
            end;

            waitbar(1, progress, 'Generating plots...');
            semilogx(absc, oord, 'LineWidth', 2);
            xlim([absc(1) absc(end)]);
            
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

            title(['\textbf{Field strength along (' mat2str(z) ', ' ...
                mat2str(y) ')}'], 'Interpreter', 'latex', 'FontSize', 24);
            xlabel('Distance to the array centre [mm]', ...
                'Interpreter', 'latex', 'FontSize', 22);
            ylabel('Field strength [dB\,V/m]', ...
                'Interpreter', 'latex', 'FontSize', 22);
            set(gca, 'FontSize', 16);

            print_plots(gcf, ['E_strength' obj.name]);

            delete(progress);

            close all;
            fprintf('\tdone\n');
        end
        
        %% Weighting function
        function w = weight(obj, mode, d)
            % Get the weight of the realized pattern
            %
            % INPUT
            %   obj:    AntArray object
            %   mode:   XY, YZ or theta, see explanation in genPattern()
            %   d:      (optional) distance of the YZ pattern [mm]
            %           OR theta angle [rad]
            
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
            % Plot the pattern of realized with the wieghting angle
            %
            % INPUT
            %   obj:    AntArray object
            %   mode:   XY, YZ or theta, see explanation in genPattern()
            %   d:      (optional) distance of the YZ pattern [mm]
            %           OR theta angle [rad]
            
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
        function obj = directivity_sph(obj)
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
        function obj = directivity_sph_alt(obj)
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
        function plotAntArray(obj)
            colors = {'m', 'c', 'r', 'g', 'b', 'k'};
            markers = {'o', '+', 'x', 's', 'd'};
            
            maxval = max(max(obj.dir));
            if maxval == 0
                error('No antenna element was initiated');
            end;
            
            figure(1);
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
            
            text(0, -length(obj.M)*0.02, ...
                [mat2str(length(obj.M)) 'x' mat2str(length(obj.M))], ...
                'Interpreter', 'latex', 'FontSize', 20);
            
            L = legend('Location', 'eastoutside');
            set(L, 'Interpreter', 'latex', 'FontSize', 20);
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
            print_plots(gcf, savname);

            close all
        end
    end
    methods (Access='public')
        %% Function to compute and plot the fields for the YZ-mode
        function E_YZ(obj, d, L, ss)
            % INPUT
            %   obj:    AntArray object
            %   d:      Distance to the array [m]
            %   L:      Side length of the plot surface [m]
            %   ss:     Step size for the plot [m]
            
            dim = round(L/ss)+1;
            ext_dim = dim+1;
            plotdata = zeros(ext_dim);    % Larger to be able to plot evth

            progress = waitbar(0, 'Computations in progress...',...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
            setappdata(progress, 'canceling', 0);

            for i=1:dim
                z = L/2 - (i-1)*ss;
                slice = plotdata(ext_dim-i,:);
                parfor j=1:dim
                    y = (j-1)*ss - L/2;
                    E = zeros(1,3);
                    [E(1), E(2), E(3)] = E_array(obj, d, y, z);
                    slice(j) = 20*log10(sqrt(sum(abs(E(:)).^2)));
                end;
                plotdata(ext_dim-i,:) = slice;
                waitbar(i/dim);
                if getappdata(progress, 'canceling')
                    close all;
                    delete(progress);
                    return;
                end;
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

            delete(progress);

            close all
        end
        
        %% Function to generate a BW plot from existing plot
        function E_BW(obj, d, theta)
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
        function E_XY(obj, d, L, ss)
            % INPUT
            %   obj:    AntArray object
            %   d:      Distance to the array (maximal plot distance) [m]
            %   L:      Side length of the plot surface [m]
            %   ss:     Step size for the plot [m]
            
            dim1 = round(L/ss)+1;
            dim2 = round(d/ss)+1;
            plotdata = zeros(dim2+1, dim1+1);    % Larger to be able to plot evth

            progress = waitbar(0, 'Computations in progress...',...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
            setappdata(progress, 'canceling', 0);

            for i=1:dim2
                x = (i-1)*ss;
                slice = plotdata(i,:);
                parfor j=1:dim1
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
                end;
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
        function E_theta(obj, d, L, ss, theta)
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

            progress = waitbar(0, 'Computations in progress...',...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
            setappdata(progress, 'canceling', 0);

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
                end;
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
        
        %% Function to optimize the antenna excitations for focusing
        function M = optArray(obj, pos, in_M, x, y, z)
            % Optimize the input array excitation using the
            % Levenberg-Marquardt algorithm for a focus at the desired
            % position with a minimum number of side lobes
            %
            % INPUT
            %   obj:    AntArray object
            %   pos:    vector with the position of the elements to tune
            %   in_M:   excitation matrix to optimize
            %   x,y,z:  desired focus position [m]
            % OUTPUT
            %   M:      matrix of optimized excitations
            
            phi0 = angle(in_M(pos));
            
            if verLessThan('matlab','8.1')
                optoptions = optimset('Display', 'iter', ...
                    'TolFun', 1e-5, 'TolX', 5e-6);
            else
                optoptions = optimset('lsqnonlin', ...
                  'Display', 'final-detailed', ...
                  'TolFun', 1e-5, 'TolX', 2.5e-6);
            end;
            
            hlen = floor(size(obj.M,1)/2);
            base = 1:hlen;
            pos_tmp = repmat(base, hlen, 1);
            pos_tmp = max(pos_tmp, pos_tmp');
            
            if mod(size(obj.M,1),2) == 0
                pos_tmp = [pos_tmp(end:-1:1, end:-1:1) pos_tmp(end:-1:1, :);
                        pos_tmp(:, end:-1:1) pos_tmp];
            else
                pos_tmp = [pos_tmp(end:-1:1, end:-1:1) base(end:-1:1)' pos_tmp(end:-1:1, :);
                        base(end:-1:1) 0 base;
                        pos_tmp(:, end:-1:1) base' pos_tmp];
            end;
            
            max_dist = max(pos_tmp(pos))*obj.spacing;
            
            E = zeros(1,3);
            [E(1), E(2), E(3)] = AntArray.E_dipole(obj.el_len, 1, obj.freq, x, max_dist, 0);
            ref = max(abs(E).^2);
            
            bound = ones(length(phi0),1)*pi;
            
            phi = lsqnonlin(...
                    @(phi)obj.optFunc(phi, pos, x, y, z, ref), ...
                    phi0, -bound, bound, optoptions); % Invoke optimizer
            
            exc = arrayfun(@(el)exp(-1j*el), phi);
            M(pos) = exc;
            
        end
        
        %% Function used for optimization
        function F = optFunc(obj, phis, pos, x, y, z, ref)
            % Compute the elements of the least-square sum
            %
            % INPUT
            %   obj:    AntArray object
            %   phis:   excitation angles
            %   pos:    excited antennas' positions
            %   x,y,z:  desired focal point
            %   ref:    reference value used to settle field bounds
            % OUTPUT
            %   F:      vector containing the least-square terms
            
            lambda = obj.c0/obj.freq;
            
            % Compute step sizes
            step = obj.opt_win(1)/obj.opt_pts;
            nb_step_x = ceil(obj.opt_win(2)/step);
            if mod(nb_step_x, 2) ~= 1
                nb_step_x = nb_step_x-1;
            end;
            step_x = obj.opt_win(2)/nb_step_x;
            
            % Init F
            F = zeros(obj.opt_pts, obj.opt_pts, nb_step_x);
            
            % Create antenna matrix
            M_opt = zeros(size(obj.M,1));
            M_opt(pos) = arrayfun(@(el)exp(-1j*el),phis);
            
            turn = (obj.opt_pts+1)/2;
            turn_x = (nb_step_x+1)/2;
            
            for i=1:nb_step_x
                xx = (i-turn_x)*step_x + x;
                for j=1:obj.opt_pts
                    zz = (turn-j)*step + z;
                    for k=1:obj.opt_pts
                        yy = (k-turn)*step + y;
                        
                        [E(1), E(2), E(3)] = obj.E_array(xx, yy, zz, M_opt);
                        En2 = sum(abs(E).^2);
                        
                        dn = (x-xx)^2 + (y-yy)^2 + (z-zz)^2;
%                         rho = xx^2 + yy^2 + zz^2;
                        
                        if dn < 4*lambda
                            gn2 = max([length(pos)-1,1.5])*ref;
                            Gn2 = 2*length(pos)*ref;
                        elseif dn > 8*lambda
                            gn2 = 0;
                            Gn2 = ceil(length(pos)/2)*ref;
                        else
                            gn2 = 0;
                            Gn2 = 2*length(pos)*ref;
                        end;
                        
                        if gn2<En2 && En2<Gn2
                            F(j, k, i) = 0;
                        else                   
                            F(j, k, i) = 2*(Gn2 - En2)*(gn2-En2)...
                                + abs(Gn2-En2)*abs(gn2-En2);
                        end;
                    end
                end
            end
            
            % define G, g, target point radius
            % Evaluate E_dipole for each defined element
            
            % Maple for Jacobian
            % Use atan of phase to avoid periodicity problems
            
            % OUTPUT MUST BE ARRAY
            
        end
        
        %% Function to compute the array input power
        function input_pwr = inpower(obj)
            [elcol, elrow] = find(obj.M ~= 0);
            
            progress = waitbar(0, 'Power computations in progress...',...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
            setappdata(progress, 'canceling', 0);
            
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
                        if z_el2 < z_el1
                            z_dist = -sqrt(z_dist);
                        else
                            z_dist = sqrt(z_dist);
                        end;
                        if y_el2 < y_el1
                            y_dist = -sqrt(y_dist);
                        else
                            y_dist = sqrt(y_dist);
                        end;
                        ang = atan2(z_dist, y_dist);

                        psi = k0*dist;

                        input_pwr(row,col,j) = Zc*k0^2*ellen^2/4/pi * ...
                            (2*(sin(psi)-psi*cos(psi))/psi^3 * (cos(ang))^2 + ...
                            ((psi^2-1)*sin(psi)+psi*cos(psi))/psi^3 * (sin(ang))^2);
                        input_pwr(row,col,j) = input_pwr(row,col,j)*I1*conj(M_tmp(j));
                    end;                        
                end;
                waitbar(i/numel(obj.M));
                if getappdata(progress, 'canceling')
                    close all;
                    delete(progress);
                    return;
                end;
            end;
            input_pwr = abs(sum(sum(sum(input_pwr(:,:,:)))));
            
            delete(progress);
        end
        
        %% Function to normalize the currents to match desired input power
        function obj = normalize(obj)
            input_pwr = inpower(obj)/2; % Divide by 2 as only radiated in half-space
            
            fact = sqrt(input_pwr/obj.pwr);
            obj.M = obj.M(:,:)./fact;
            obj.normalized = 1;
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
                warning(['File ' savname ' not found, skipped']);
                filename = 0;
            end; 
        end
    end
    methods (Static, Access = 'public')
        %% Get field value at a given position
        function val = getVal(x, y, file, folder)
            % Extract the value of the electric field at a given position
            % from the given 'dat' file
            % 
            % INPUT
            %    x,y:    Coordinates of interest [mm]
            %    file:   Name of the source file [.dat extension]
            %    folder: (optional) Name of the source folder
           
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
                warning('Specified point might be out of plotted area');
            end;
            
        end
    end
end