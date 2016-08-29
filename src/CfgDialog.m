%CfgDialog  Class for displaying configuration dialog
%
%The configuration dialog permits to input the parameters inherent to a
%class, for example AntArray.
%
%CfgDialog Properties
%   handler         handler of the main configuration window
%   ptr             handler of an invisible configuration window
%   fields          list of parameters to be set
%
%CfgDialog Methods
%   CfgDialog       constructor 
%   waitClosing     wait for the closing of the dialog
%   setDft          set default values for the fields
%   getFields       get input data
%   cancelAction    (static) action when pressing the cancel button
%   confirmAction   (static) action when pressing the confirm button
%   closeAction     (static) action when closing the dialog
%
%Use the DOC command for detailed explanations

% Copyright 2015-2016, Antoine JUCKLER. All rights reserved

classdef CfgDialog < handle
    
    properties (GetAccess='public', SetAccess='private')
        handler;    % Handler of the main dialog
        ptr;        % Handler of an invisible dialog used to store data
        fields;     % list of parameters to be set
    end;
    methods
        %% Constructor
        function obj = CfgDialog(list)
            %CFGDIALOG constructor
            %
            % obj = CFGDIALOG(list)
            %
            % INPUT
            %   list:   list of field names
            % OUTPUT
            %   obj:    CfgDialog object
            
            sizeW = [300 20*length(list)+50];
            screensize = get(0, 'ScreenSize');
            xpos = ceil((screensize(3)-sizeW(1))/2);
            ypos = ceil((screensize(4)-sizeW(2))/2);
            obj.ptr = figure('Name', 'Ptr', ...
                'Visible', 'off');
            obj.handler = figure('Name', 'Array parameters', ...
                'Units', 'pixels', ...
                'Position', [xpos ypos sizeW], ...
                'Visible', 'off', ...
                'DeleteFcn', {@CfgDialog.closeAction, obj});
            obj.fields = struct();

            for i=1:length(list)
                uicontrol('Parent', obj.handler, ...
                    'Style', 'text', ...
                    'Units', 'pixels', ...
                    'Position', [0 sizeW(2)-20*i sizeW(1)/2 20], ...
                    'FontSize', 9, ...
                    'HorizontalAlignment', 'right', ...
                    'String', [list{i} ' ']);
                uicontrol('Parent', obj.handler, ...
                    'Style', 'edit', ...
                    'Units', 'pixels', ...
                    'Position', [sizeW(1)/2 sizeW(2)-20*i sizeW(1)/2 20], ...
                    'FontSize', 9, ...
                    'tag', list{i});
                obj.fields.(list{i}) = [];
            end;

            btn_width = sizeW(1)*.4;
            btn_left = (sizeW(1)/2 - btn_width)/2;
            uicontrol('Parent', obj.handler, ...
               'Position', [btn_left 15 btn_width 25], ...
               'String', 'Cancel', ...
               'FontSize', 9, ...
               'BackgroundColor', [.99 .18 .18], ...
               'Callback', {@CfgDialog.cancelAction, obj});
            uicontrol('Parent', obj.handler, ...
               'Position', [btn_left+sizeW(1)/2 15 btn_width 25], ...
               'String', 'Confirm', ...
               'FontSize', 9, ...
               'BackgroundColor', [.09 .79 .06], ...
               'Callback', {@CfgDialog.confirmAction, obj});
            
            set(obj.handler, 'MenuBar', 'none');
            set(obj.handler, 'Visible', 'on');
            drawnow;
            
            handles = guihandles(obj.handler);
            fl = fieldnames(handles);
            for i=1:length(fl)
                if ~isfield(obj.fields, fl{i})
                    handles = rmfield(handles, fl{i});
                end;
            end;
            guidata(obj.handler, handles);
        end
        
        %% Wait for action on main dialog window
        function waitClosing(obj)
            %WAITCLOSING wait for an user action on the main dialog window
            %
            % [ ] = WAITCLOSING(obj)
            %
            % INPUT
            %   obj:    CfgDialog object
            %
            % See also CANCELACTION CONFIRMACTION CLOSEACTION
            
            waitfor(obj.handler);
        end
        
        %% Set field default values
        function setDft(obj, dft)
            %SETDFT set default values of the fields to be input on the
            %dialog
            %
            % [ ] = SETDTF(obj, dft)
            %
            % INPUT
            %   obj:    CfgDialog object
            %   dft:    structure of field name-value pairs
            
            fl = fieldnames(dft);
            for i=1:length(fl)
                if ~isfield(obj.fields, fl{i})
                    dft = rmfield(dft, fl{i});
                end;
            end;

            els = guidata(obj.handler);
            fl = fieldnames(dft);
            for i=1:length(fl)
                if isfield(els, fl{i})
                    set(els.(fl{i}), 'Value', dft.(fl{i}));
                end;
            end;
            drawnow;
        end
        
        %% Get input values
        function vals = getFields(obj)
            %GETFIELDS get the input values
            %
            % The invisible figure containing the input data will be
            % deleted in the process.
            %
            % vals = getFields(obj)
            %
            % INPUT
            %   obj:    CfgDialog object
            % OUTPUT
            %   vals:   struct of field name-value pairs
            
            vals = get(obj.ptr, 'UserData');
            close(obj.ptr);
        end;
    end;
    methods (Static)
        %% Action when pressing the cancel button
        function cancelAction(hObject, callbackdata, obj)
            %CANCELACTION cancel all the input data and close the window
            %
            % [ ] = CANCELACTION(hObject, callbackdata, obj)
            %
            % INPUT
            %   hObject:        handler of object calling the function
            %   callbackdata:   additional callback information
            %   obj:            CfgDialog object
            
            fl = fieldnames(obj.fields);
            for i=1:length(fl)
                obj.fields.(fl{i}) = [];
            end;
            
            set(obj.ptr, 'UserData', obj.fields);
            delete(obj.handler);
        end
        
        %% Action when pressing the confirm button
        function confirmAction(hObject, callbackdata, obj)
            %CONFIRMACTION copy the input data to the invisible handle and
            %close the window
            %
            % [ ] = CONFIRMACTION(hObject, callbackdata, obj)
            %
            % INPUT
            %   hObject:        handler of object calling the function
            %   callbackdata:   additional callback information
            %   obj:            CfgDialog object
            
            set(obj.handler, 'UserData', 1);
            els = guidata(obj.handler);
            fl = fieldnames(els);

            % Remove wrong fields from fl
            for i=1:length(fl)
                if ~isfield(obj.fields, fl{i})
                    els = rmfield(els, fl{i});
                end;
            end;
            
            fl = fieldnames(els);
            % Parse inputted values
            for i=1:length(fl)
                inval = get(els.(fl{i}), 'string');
                if isempty(inval)
                    if ~isempty(get(els.(fl{i}), 'value'))
                        obj.fields.(fl{i}) = get(els.(fl{i}), 'value');
                    else
                        obj.fields.(fl{i}) = [];
                    end;
                elseif strcmp(fl{i}, 'name')
                    obj.fields.(fl{i}) = inval;
                else
                    obj.fields.(fl{i}) = str2double(inval);
                end;
            end;
            set(obj.ptr, 'UserData', obj.fields);

            delete(obj.handler);
        end
        
        %% Action when closing the window
        function closeAction(hObject, callbackdata, obj)
            %CLOSEACTION calls CANCELACTION and close the window
            %
            % [ ] = CLOSEACTION(hObject, callbackdata, obj)
            %
            % INPUT
            %   hObject:        handler of object calling the function
            %   callbackdata:   additional callback information
            %   obj:            CfgDialog object

            if isempty(get(obj.handler, 'UserData'))
                CfgDialog.cancelAction(hObject, callbackdata, obj);
            end;
            
            if ishandle(obj.handler)
                delete(obj.handler);
            end;
        end
    end;
end