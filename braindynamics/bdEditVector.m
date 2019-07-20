function outdata = bdEditVector(data,name,columnName)
    %bdEditVector Dialog box for editing vector data.
    %Usage:
    %    outdata = bdEditVector(indata,name,columnName)
    %where
    %    'indata' is the initial vector data (nx1)
    %    'name' is the title of the dialog box (text)
    %    'columnName' is the title of the data column (text)
    %Returns the edited data in 'outdata'. If the user cancels the operation
    %then the returned data is identical to the initial data. 

    % init return data (in case user cancels)
    outdata = data;

    % remember the orginal shape of our data
    orgsize = size(data);
    
    % ensure incoming data is a column vector
    data = reshape(data,[],1);
    
    % length of vector
    n = numel(data);
    
    % construct dialog box
    fig = figure('Units','pixels', ...
        'Position',[randi(300,1,1) randi(300,1,1), 400, 300], ...
        'MenuBar','none', ...
        'Name',name, ...
        'NumberTitle','off', ...
        'ToolBar', 'none', ...
        'Resize','off');

    % axes
    ax1 = axes('parent',fig, 'Units','pixels', 'Position',[185 205 200 70]);      % axes for bar graph
    ax2 = axes('parent',fig, 'Units','pixels', 'Position',[185  80 200 70]);      % axes for histogram

    % bar graph
    barg = bar(data, 'parent',ax1);
    xlim(ax1,[0.5 n+0.5]);
    xlabel('parent',ax1,'index');
    ylabel('parent',ax1,'value');
    title(ax1,['Values of ', columnName]);
    
    % histogram
    hist = histogram(data,'parent',ax2, 'Normalization','probability');
    xlabel('parent',ax2,'value');
    ylabel('parent',ax2,'proportion');
    title(ax2,['Histogram of ', columnName]);
    
    % data table
    tbl = uitable(fig,'Position',[10 10 125, 280], ...
        'Data',data, ...
        'ColumnName',{columnName}, ...
...        'ColumnWidth',{75}, ...
        'ColumnEditable',true, ...
        'CellEditCallback', @(src,~) CellEditCallback(src,barg,hist));
    
    % 'Cancel' button
    uicontrol('Style','pushbutton', ...
        'String','Cancel', ...
        'HorizontalAlignment','center', ...
        'FontUnits','pixels', ...
        'FontSize',12, ...
        'Parent', fig, ...
        'Callback', @(~,~) delete(fig), ...
        'Position',[255 10 60 20]);
    
    % 'OK' button
    uicontrol('Style','pushbutton', ...
        'String','OK', ...
        'HorizontalAlignment','center', ...
        'FontUnits','pixels', ...
        'FontSize',12, ...
        'Parent', fig, ...
        'Callback', @(~,~) OKCallback(), ...
        'Position',[325 10 60 20]);
    
    % wait for close button to delete the figure
    uiwait(fig);
    
    function CellEditCallback(tbl,barg,hist)
        %disp('CellEditCallback');
        data = get(tbl,'data');
        % update bar graph
        barg.YData = data;
        % update histogram
        hist.Data = data;
        hist.BinLimitsMode='auto';
    end

    % Callback for OK button
    function OKCallback()
        outdata = reshape(data,orgsize);    % return the new data in it original shape
        delete(fig);                        % close the dialog box
    end

%     % Callback for edit box
%     function EditCallback(hObj,~)
%         % get the incoming value
%         val = str2double(hObj.String);
%         if isnan(val)
%             dlg = errordlg(['''', hObj.String, ''' is not a valid number'],'Invalid Number','modal');
%             val = hObj.Value;           % restore the previous value                
%             uiwait(dlg);                % wait for dialog box to close
%         else
%             hObj.Value = val;           % remember the new value
%         end            
% 
%         % update the edit box string
%         hObj.String = num2str(val,'%g');
%     end

end
        

