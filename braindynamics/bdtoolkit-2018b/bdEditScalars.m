%bdEditScalars  Dialog box for editing scalar values.
%Usage:
%    outdata = bdEditScalars(pardef,name,descr)
%where
%    "pardeg" is a cell array of {value,"text"} pairs
%    "name" is the title of the dialog box (text)
%    "descr" is the description of the content (text)
%Returns the edited data in "outdata". If the user cancels the operation
%then the returned as empty.

function outdata = bdEditScalars(pardef,name,descr)
    outdata = [];

    % number of rows in pardef
    npardef = size(pardef,1);
    
    % array of edit box widgets
    editbox = gobjects(npardef);

    % row geometry
    rowh = 22;
    boxh = 20;
    boxw = 50;
    
    % open dialog box
    height = npardef*rowh + 100; 
    dlg = dialog("Position",[300 300 250 height],"Name",name);

    % first line
    ypos = height - 1.5*rowh;
    
    % heading
    uicontrol("Parent",dlg,...
        "Style","text",...
        "Position",[20 ypos 210 rowh],...
        "String",descr, ...
        "FontSize",12, ...
        "HorizontalAlignment","left", ...
    ... "BackgroundColor","r", ...
        "FontWeight","bold");

    % for each pardef extry
    for indx=1:npardef
        % next line
        ypos = ypos - rowh;

        % edit box
        editbox(indx) = uicontrol("Parent",dlg,...
            "Style","edit",...
            "Position",[20 ypos boxw boxh],...
            "String",num2str(pardef{indx,1}),...
            "FontSize",12, ...
            "HorizontalAlignment","center", ...
        ... "BackgroundColor","r", ...           
            "Callback",@editbox_callback);

        % number of neurons (text label)
        uicontrol("Parent",dlg,...
            "Style","text",...
            "Position",[75 ypos 155 boxh],...
            "String",pardef{indx,2}, ...
            "FontSize",12, ...
            "HorizontalAlignment","left", ...
        ... "BackgroundColor","r", ...
            "FontWeight","normal");
    end

    % next line
    ypos = ypos - 1.5*rowh;
    
    % syntax error (text)    
    err = uicontrol("Parent",dlg,...
        "Style","text",...
        "Position",[20 ypos 200 25],...
        "String","Syntax Error", ...
        "FontSize",12, ...
        "HorizontalAlignment","left", ...
        "ForegroundColor","r", ...
        "Visible","off", ...
        "FontWeight","normal");
    
    % next line
    ypos = ypos - rowh;

    % CANCEL button
    uicontrol("Parent",dlg,...
        "Position",[20 ypos 75 25],...
        "String","Cancel",...
        "Callback",@(~,~) cancel_callback );

    % CONTINUE button
    btn = uicontrol("Parent",dlg,...
        "Position",[250-95 ypos 75 25],...
        "String","Continue",...
        "Callback",@(~,~) delete(dlg) );
    
    % force the editbox values into outdata
    editbox_callback([],[]);
    
    % Wait for dialog to close
    uiwait(dlg);
   
    function editbox_callback(~,~)
        for indx = 1:npardef
            val = str2num(editbox(indx).String);
            if isempty(val)
                % invalid number
                btn.Enable = "off";      % disable the CONTINUE button
                err.Visible = "on";      % show the syntax error text
                return
            else
                outdata(indx) = val;
            end
        end
        btn.Enable = "on";              % enable the CONTINUE button
        err.Visible = "off";            % hide the syntax error text
    end

    function cancel_callback(~,~)
        outdata = [];
        delete(dlg);
    end
    
end
