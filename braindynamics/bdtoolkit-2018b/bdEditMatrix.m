function outdata = bdEditMatrix(data,name)
    %bdEditMatrix Dialog box for editing matrix data.
    %Usage:
    %    outdata = bdEditMatrix(indata,name)
    %where
    %    'indata' is the initial matrix data (nxn)
    %    'name' is the title of the dialog box (text)
    %Returns the edited data in 'outdata'. If the user cancels the operation
    %then the returned data is identical to the initial data. 

    % init return data
    outdata = data;

    % construct dialog box
    fig = figure('Units','pixels', ...
        'Position',[randi(300,1,1) randi(300,1,1), 600, 270], ...
        'MenuBar','none', ...
        'Name',name, ...
        'NumberTitle','off', ...
        'ToolBar', 'none', ...
        'Resize','off');
    
    % set data cursor mode
    dcm = datacursormode(fig);
    set(dcm, 'DisplayStyle','datatip', 'Enable','on','UpdateFcn',@UpdateFcn);
    
    % axis for image
    ax = axes('parent',fig, 'Units','Pixels', 'Position',[390 60 200 200]);
   
    % image
    img = imagesc('CData',data, 'parent',ax);
    axis image ij
 
    % data table
    tbl = uitable(fig,'Position',[10 10 350 250], ...
        'Data',data, ...
        'ColumnWidth',{50}, ...
        'ColumnEditable',true, ...
        'CellEditCallback', @(src,~) CellEditCallback(src,img));
    
    % 'Cancel' button
    uicontrol('Style','pushbutton', ...
        'String','Cancel', ...
        'HorizontalAlignment','center', ...
        'FontUnits','pixels', ...
        'FontSize',12, ...
        'Parent', fig, ...
        'Callback', @(~,~) delete(fig), ...
        'Position',[460 10 60 20]);

    % 'OK' button
    uicontrol('Style','pushbutton', ...
        'String','OK', ...
        'HorizontalAlignment','center', ...
        'FontUnits','pixels', ...
        'FontSize',12, ...
        'Parent', fig, ...
        'Callback', @(~,~) OKCallback, ...     % delete the dialog box
        'Position',[530 10 60 20]);
    
    % wait for dialog box to be deleted
    uiwait(fig);
    
    % Callback for datatip
    function txt = UpdateFcn(~,evnt)
        pos = evnt.Position;
        val = data(pos(2),pos(1));
        txt = {num2str(pos,'X=%d, Y=%d'),num2str(val,'%0.4g')};
    end
    
    % Callback for uitable
    function CellEditCallback(tbl,img)
        %disp('CellEditCallback');
        data = get(tbl,'data');                 % Retrieve the new data from uitable
        set(img,'CData',data);                  % Update the image data
        %set(dcm, 'DisplayStyle','datatip');     % Hack to force the datatip to update (doesnt help)
    end

    % Callback for OK button
    function OKCallback()
        outdata = data;     % return the new data
        delete(fig);        % close the dialog box
    end
    
end
        

