% Plot the temporal and spatial one-dimensional Gabor functions
function varargout = temp_space_rfs(varargin)
% TEMP_SPACE_RFS M-file for temp_space_rfs.fig
%      TEMP_SPACE_RFS, by itself, creates a new TEMP_SPACE_RFS or raises the existing
%      singleton*.
%
%      H = TEMP_SPACE_RFS returns the handle to a new TEMP_SPACE_RFS or the handle to
%      the existing singleton*.
%
%      TEMP_SPACE_RFS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEMP_SPACE_RFS.M with the given input arguments.
%
%      TEMP_SPACE_RFS('Property','Value',...) creates a new TEMP_SPACE_RFS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before temp_space_rfs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to temp_space_rfs_OpeningFcn via varargin.
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @temp_space_rfs_OpeningFcn, ...
                   'gui_OutputFcn',  @temp_space_rfs_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before temp_space_rfs is made visible.
function temp_space_rfs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to temp_space_rfs (see VARARGIN)

% Choose default command line output for temp_space_rfs
handles.output = hObject;


eta_t = 8; %in cycles/sec
eta_t = 2*pi*eta_t; %circular frequency in rad/sec
sigma_t = 31/1000;
delay = 86/1000;

Ct = 1/(sqrt(2*pi)*sigma_t);

t = -0.025:0.001:0.250;

g_t1 = exp(-(t-delay).^2/(2*sigma_t^2)).*cos(eta_t*(t-delay));
g_t1 = Ct*g_t1;
g_t2 = exp(-(t-delay).^2/(2*sigma_t^2)).*sin(eta_t*(t-delay));
g_t2 = Ct*g_t2;

line('Parent',handles.axes1,'XData',t,'YData',g_t1,'Color','k');
line('Parent',handles.axes1,'XData',t,'YData',g_t2,'Color','r');
set(handles.axes1,'XLim',[-0.025 0.250],'YLim',[-10 15]);
xlabel(handles.axes1,'time (s)');
ylabel(handles.axes1,'firing rate (arbitrary units)');

%circular frequency corresponding to 4.2 cycles/deg
k_x = 2*pi*4.2;

%in degrees
s_x = 0.1; 

dx = 0.004;            %units are degrees
nx = 128;

x = ((-nx/2+1)*dx:dx:nx/2*dx);  %units are degrees, 128 points

%make sure that we are in line format
x = x(:)';

%compute the spatial filters

%length of the spatial vector
n = length(x);
r_es = zeros(1,n);
Cx = 1/(sqrt(2*pi)*s_x);
r_es(1,1:n) = exp(-x.^2/(2*s_x^2)).*cos(k_x*x);
r_es = Cx*r_es;


r_os = zeros(1,n);
r_os(1,1:n) = exp(-x.^2/(2*s_x^2)).*sin(k_x*x);
r_os = Cx*r_os;

line('Parent',handles.axes2,'XData',x,'YData',r_es,'Color','k');
line('Parent',handles.axes2,'XData',x,'YData',r_os,'Color','r');
set(handles.axes2,'XLim',[-0.3 0.3]);
xlabel('space (deg)');
ylabel('firing rate (arbitrary units)');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes temp_space_rfs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = temp_space_rfs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function save_fig_Callback(hObject, eventdata, handles)
% hObject    handle to save_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

print(handles.figure1,'-depsc2','temp_space_rfs.eps');

