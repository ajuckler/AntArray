classdef CfgDialog < handle
    
    properties (GetAccess='public', SetAccess='private')
        handler;
        ptr;
        fields;
    end;
    methods
        function obj = CfgDialog(list)
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
                'DeleteFcn', {@close, obj});
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
               'Callback', {@cancel, obj});
            uicontrol('Parent', obj.handler, ...
               'Position', [btn_left+sizeW(1)/2 15 btn_width 25], ...
               'String', 'Confirm', ...
               'FontSize', 9, ...
               'BackgroundColor', [.09 .79 .06], ...
               'Callback', {@confirm, obj});
            
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
        
        function waitClosing(obj)
            waitfor(obj.handler);
        end
        
        function setDft(obj, dft)
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
        
        function vals = getFields(obj)
            vals = get(obj.ptr, 'UserData');
            close(obj.ptr);
        end;
    end;
    methods (Access='private')
        function cancel(hObject, callbackdata, obj)
            fl = fieldnames(obj.fields);
            for i=1:length(fl)
                obj.fields.(fl{i}) = [];
            end;
            
            set(obj.ptr, 'UserData', obj.fields);
            obj.close();
        end

        function confirm(hObject, callbackdata, obj)
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

            obj.close();
        end
        
        function close(obj)
            if get(obj.handler, 'UserData') ~= 1
                obj.cancel();
            end;
            delete(obj.handler);
        end
    end;
end