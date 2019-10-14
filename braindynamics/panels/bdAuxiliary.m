classdef bdAuxiliary < bdPanel
    %bdAuxiliary Display panel for plotting model-specific functions.
    properties (Constant)
        title = "Auxiliary"
    end
    
    properties
        ax              % Handle to the plot axes
        UserData        % User data returned by the auxiliary function     
    end
    
    properties (Access=private)
        holdmenu        % handle to HOLD menu item
        submenu = []    % handle to subpanel selector menu item
        listener        % handle to our listener object
    end
    
    methods
        
        function this = bdAuxiliary(tabgroup,control)
            % Construct a new Auxiliary Panel in the given tabgroup

            % initialise the base class (specifically this.menu and this.tab)
            this@bdPanel(tabgroup);
            
            % assign default values to missing options in sys.panels.bdTimePortrait
            control.sys.panels.bdAuxiliary = bdAuxiliary.syscheck(control.sys);

            % configure the pull-down menu
            this.menu.Label = control.sys.panels.bdAuxiliary.title;
            this.InitHoldMenu(control);
            this.InitExportMenu(control);
            this.InitCloseMenu(control);

            % configure the panel graphics
            this.tab.Title = control.sys.panels.bdAuxiliary.title;
            this.InitSubpanel(control);
            
            % listen to the control panel for redraw events
            this.listener = addlistener(control,"redraw",@(~,~) this.redraw(control));    
        end
        
        function delete(this)
            % Destructor
            delete(this.listener)
        end
         
    end
    
    methods (Access=private)
        
        % Initialise the HOLD menu item
        function InitHoldMenu(this,control)
             % get the hold menu setting from sys.panels options
            if control.sys.panels.bdAuxiliary.hold
                holdcheck = "on";
            else
                holdcheck = "off";
            end
            
            % construct the menu item
            this.holdmenu = uimenu(this.menu, ...
                "Label","Hold", ...
                "Checked",holdcheck, ...
                "Callback", @HoldMenuCallback );

            % Menu callback function
            function HoldMenuCallback(menuitem,~)
                switch menuitem.Checked
                    case "on"
                        menuitem.Checked="off";
                    case "off"
                        menuitem.Checked="on";
                end
                % redraw this panel
                this.redraw(control);
            end
        end
        
        % Initialise the EXPORT menu item
        function InitExportMenu(this,~)
            % construct the menu item
            uimenu(this.menu, ...
               "Label","Export Figure", ...
               "Callback",@callback);
           
            function callback(~,~)
                % Construct a new figure
                fig = figure();    
                
                % Change mouse cursor to hourglass
                set(fig,"Pointer","watch");
                drawnow;
                
                % Copy the plot data to the new figure
                axnew = copyobj(this.ax,fig);

                % Allow the user to hit everything in ax1new
                objs = findobj(axnew,"-property", "HitTest");
                set(objs,"HitTest","on");
                
                % Change mouse cursor to arrow
                set(fig,"Pointer","arrow");
                drawnow;
            end
        end

        % Initialise the CLOSE menu item
        function InitCloseMenu(this,~)
            % construct the menu item
            uimenu(this.menu, ...
                   "Label","Close", ...
                   "Callback",@(~,~) this.close());
        end
        
        % Initialise the subpanel
        function InitSubpanel(this,control)
            % construct the subpanel
            [this.ax,cmenu] = bdPanel.Subpanel(this.tab,[0 0 1 1],[0 0 1 1]);
            
            % default axes for the case of no auxiliary function 
            title(this.ax,"No Auxiliary Functions");
            text(0.5,0.5,"No auxiliary plotting functions are defined for this model", ...
                "HorizontalAlignment","center", "Parent",this.ax);

            % construct the selector menu for the auxiliary functions
            naux = numel(control.sys.panels.bdAuxiliary.auxfun);
            for indx=1:naux
                UserData.auxfun = control.sys.panels.bdAuxiliary.auxfun{indx};
                UserData.label = func2str(UserData.auxfun);
                UserData.rootmenu = cmenu;
                menuitem = uimenu("Parent",cmenu, ...
                    "Label",UserData.label, ...
                    "Checked","off", ...
                    "Tag","auxmenu", ...
                    "UserData",UserData, ...
                    "Callback",@callback);                
                if indx==1
                    menuitem.Checked = "on";
                    this.submenu = menuitem;
                end
            end
            
            
            % Callback function for the subpanel selector menu
            function callback(menuitem,~)
                % check "on" the selected menu item and check "off" all others
                bdPanel.SelectorCheckItem(menuitem);
                
                % update our handle to the selected menu item
                this.submenu = menuitem;
                
                % clear the axis, reset to defaults
                cla(this.ax,"reset");
                
                % redraw the panel
                this.redraw(control);
            end
        end
        
        % Redraw the data plots
        function redraw(this,control)
            %disp("bdAuxiliary.redraw()")
            
            % if the submenu is empty (because no auxiliary functions were defined) then break
            if isempty(this.submenu)
                return
            end
            
            % make the axis current
            axes(this.ax);

            % if "hold" menu is checked then ...
            switch this.holdmenu.Checked
                case "off"
                    % Clear the plot axis
                    cla(this.ax);
            end

            % get the details of the currently selected plot function
            auxfun  = this.submenu.UserData.auxfun;
            
            % Execute the auxiliary plot function.
            % The type of the solver function determines the parameters we call it with. 
            switch control.solvertype
                case "odesolver"
                    % case of an ODE solver (eg ode45)
                    parcell = struct2cell(control.par)";
                    %this.auxdata = feval(auxfun,this.ax,control.sys.tval,control.sol,parcell{:});
                    if nargout(auxfun)==0
                        auxfun(this.ax,control.sys.tval,control.sol,parcell{:});
                    else
                        this.UserData = auxfun(this.ax,control.sys.tval,control.sol,parcell{:});
                    end

                case "ddesolver"
                    % case of a DDE solver (eg dde23)
                    lagcell = struct2cell(control.lag)";
                    parcell = struct2cell(control.par)";
                    allcell = {lagcell{:} parcell{:}};
                    %feval(auxfun,this.ax,control.sys.tval,control.sol,allcell{:});
                    if nargout(auxfun)==0
                        auxfun(this.ax,control.sys.tval,control.sol,allcell{:});
                    else
                        this.UserData = auxfun(this.ax,control.sys.tval,control.sol,allcell{:});
                    end

                case "sdesolver"
                    % case of an SDE solver
                    parcell = struct2cell(control.par)";
                    %feval(auxfun,this.ax,control.sys.tval,control.sol,parcell{:});
                    %this.auxdata = feval(auxfun,this.ax,control.sys.tval,control.sol,parcell{:});
                    if nargout(auxfun)==0
                        auxfun(this.ax,control.sys.tval,control.sol,parcell{:});
                    else
                        this.UserData = auxfun(this.ax,control.sys.tval,control.sol,parcell{:});
                    end                    
            end        
        end

    end
    
    methods (Static)
        function syspanel = syscheck(sys)
            % Assign default values to missing fields in sys.panels.bdAuxiliary

            % Default panel settings
            syspanel.title = bdAuxiliary.title;
            syspanel.hold = false;
            syspanel.auxfun = [];
            
            % Nothing more to do if sys.panels.bdAuxiliary is undefined
            if ~isfield(sys,"panels") || ~isfield(sys.panels,"bdAuxiliary")
                return;
            end
            
            % sys.panels.bdAuxiliary.title
            if isfield(sys.panels.bdAuxiliary,"title")
                syspanel.title = sys.panels.bdAuxiliary.title;
            end
            
            % sys.panels.bdAuxiliary.hold
            if isfield(sys.panels.bdAuxiliary,"hold")
                syspanel.hold = sys.panels.bdAuxiliary.hold;
            end
            
            % sys.panels.bdAuxiliary.auxfun
            if isfield(sys.panels.bdAuxiliary,"auxfun")
                syspanel.auxfun = sys.panels.bdAuxiliary.auxfun;
            end
        end
        
        function auxdefault(ax,varargin)
            text(0.5,0.5,"No auxiliary plotting functions are defined for this system", ...
                "HorizontalAlignment","center", "Parent",ax);
            title(ax,"No Auxiliary Functions");
        end
    end
    
end
