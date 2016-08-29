%WaitDialog  Class for displaying current optimisation status
%
%WaitDialog Properties
%   handler         handler of the main window
%   main_string     main string to display
%   sub_string      additional information
%
%WaitDialog Methods
%   WaitDialog      constructor 
%   delete          delete the object
%   terminate   	check whether abort button has been pressed
%   setMainString   set main string
%   setSubString    set additional information string
%   abort           (static) action when pressing the abort button
%
%Use the DOC command for detailed explanations

% Copyright 2015-2016, Antoine JUCKLER. All rights reserved

classdef WaitDialog < handle
    
    properties (GetAccess='public', SetAccess='private')
        handler;        % Handler of the main window
        main_string;    % main string
        sub_string;     % additional information string
    end;
    methods
        %% Constructor
        function obj = WaitDialog()            
            %WAITDIALOG constructor
            %
            % obj = WAITDIALOG()
            %
            % OUTPUT
            %   obj:    WaitDialog object
            
            sizeW = [300 170];
            screensize = get(0, 'ScreenSize');
            xpos = ceil((screensize(3)-sizeW(1))/2);
            ypos = ceil((screensize(4)-sizeW(2))/2);
            obj.handler = dialog('Name', '', ...
                'Position', [xpos ypos sizeW], ...
                'Visible', 'off', ...
                'WindowStyle', 'normal');

            txt_width = sizeW(1) - 40;
            uicontrol('Parent', obj.handler, ...
               'Style', 'text', ...
               'Position', [20 135 txt_width 25], ...
               'FontSize', 11, ...
               'String', 'Optimization in progress...');

            btn_width = 70;
            btn_left = (sizeW(1) - btn_width)/2;
            uicontrol('Parent', obj.handler, ...
               'Position', [btn_left 15 btn_width 25], ...
               'String', 'Abort', ...
               'FontSize', 9, ...
               'Callback', {@WaitDialog.abort, obj});
           setappdata(obj.handler, 'abort', 0);

           uicontrol('Parent', obj.handler, ...
                'Style', 'text', ...
                'BackgroundColor', 'black', ...
                'Visible', 'on', ...
                'Position', [20 50 txt_width 80]);
            uicontrol('Parent', obj.handler, ...
                'Style', 'text', ...
                'BackgroundColor', 'white', ...
                'Visible', 'on', ...
                'Position', [21 51 txt_width-2 78]);

            obj.main_string = uicontrol('Parent', obj.handler, ...
                'Style', 'text', ...
                'Position', [21 110 txt_width-2 19], ...
                'BackgroundColor', 'white', ...
                'HorizontalAlignment', 'left', ...
                'FontSize', 9, ...
                'String', 'Starting...');

            sub_width = txt_width-2-(35-21);
            obj.sub_string = uicontrol('Parent', obj.handler, ...
                'Style', 'text', ...
                'Position', [35 90 sub_width 20], ...
                'BackgroundColor', 'white', ...
                'HorizontalAlignment', 'left', ...
                'FontSize', 9, ...
                'String', '');
            
            set(obj.handler, 'Visible', 'on');
            drawnow;
        end;
        
        %% Delete object
        function delete(obj)
            %DELETE delete the window object
            %
            % [ ] = DELETE(obj)
            %
            % INPUT
            %   obj:    WaitDialog object
            
            delete(obj.handler);
        end;
        
        %% Check whether abort button has been pressed
        function terminate(obj)
            %TERMINATE check whether the abort button has been pressed,
            %so yes, the object is deleted and an exception is thrown
            %
            % [ ] = TERMINATE(obj)
            %
            % INPUT
            %   obj:    WaitDialog object
            
            if ~getappdata(obj.handler, 'abort')
                return;
            end;
            delete(obj);
            throw(MException('MyERR:Terminated', 'Operation aborted by user'));
        end;
        
        %% Set main string
        function setMainString(obj, str)
            %SETMAINSTRING set the main information string, the additional
            %information string is cleared in the process
            %
            % [ ] = SETMAINSTRING(obj, str)
            %
            % INPUT
            %   obj:    WaitDialog object
            %   str:    main information string
            
            set(obj.main_string, 'String', str);
            set(obj.sub_string, 'String', '');
            drawnow;
        end;
        
        %% Set substring
        function setSubString(obj, str)
            %SETSUBSTRING set the additional information string
            %
            % [ ] = SETSUBSTRING(obj, str)
            %
            % INPUT
            %   obj:    WaitDialog object
            %   str:    additional information string
            
            set(obj.sub_string, 'String', str);
            drawnow;
        end;
    end;
    methods (Static, Access='public')
        %% Action when pressing the abort button
        function abort(~, ~, obj)
            %ABORT set the abort flag high
            %
            % [ ] = ABORT(~, ~, obj)
            %
            % INPUT
            %   obj:    WaitDialog object
            
            setappdata(obj.handler, 'abort', 1);
        end;
    end;
end