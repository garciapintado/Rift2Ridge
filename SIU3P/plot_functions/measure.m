function varargout = measure(varargin)
% MEASURE MATLAB code for measure.fig
%      MEASURE, by itself, creates a new MEASURE or raises the existing
%      singleton*.
%
%      H = MEASURE returns the handle to a new MEASURE or the handle to
%      the existing singleton*.
%
%      MEASURE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEASURE.M with the given input arguments.
%
%      MEASURE('Property','Value',...) creates a new MEASURE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before measure_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to measure_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help measure

% Last Modified by GUIDE v2.5 30-Aug-2017 10:34:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @measure_OpeningFcn, ...
                   'gui_OutputFcn',  @measure_OutputFcn, ...
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


% --- Executes just before measure is made visible.
function measure_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to measure (see VARARGIN)

% Choose default command line output for measure
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
if isempty(varargin)
    figure(1)
else
    figure(varargin{1})
end
ax = gca;
p = get(ax,'children');
cax = get(ax,'CLim');
hax = handles.plot;
cla(hax)
copyobj(p,hax)
caxis(hax,cax);

% UIWAIT makes measure wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = measure_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in xb.
function xb_Callback(hObject, eventdata, handles)
% hObject    handle to xb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text2,'string','Distance')
hax = handles.plot;
clear_points
x = [0 0];
y = [0 0];
for n = 1:2
    [x(n),y(n)] = ginput(1);
    hold(hax,'on')
    plot(x(n),y(n),'+k')
end
last_p = get(gca, 'children');
set(last_p(1:2),'Tag','deleteme')
dist = num2str(abs(x(1)-x(2)));
set(handles.dist_txt,'string',num2str(dist))

% --- Executes on button press in yb.
function yb_Callback(hObject, eventdata, handles)
% hObject    handle to yb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text2,'string','Distance')
hax = handles.plot;
clear_points
x = [0 0];
y = [0 0];
for n = 1:2
    [x(n),y(n)] = ginput(1);
    hold(hax,'on')
    plot(x(n),y(n),'+k')
end
last_p = get(gca, 'children');
set(last_p(1:2),'Tag','deleteme')
dist = num2str(abs(y(1)-y(2)));
set(handles.dist_txt,'string',num2str(dist))

% --- Executes on button press in bline.
function bline_Callback(hObject, eventdata, handles)
% hObject    handle to bline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to yb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text2,'string','Distance')
hax = handles.plot;
clear_points
x = [0 0];
y = [0 0];
for n = 1:2
    [x(n),y(n)] = ginput(1);
    hold(hax,'on')
    plot(x(n),y(n),'+k')
end
last_p = get(gca, 'children');
set(last_p(1:2),'Tag','deleteme')
dist = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2);
set(handles.dist_txt,'string',num2str(dist))

function dist_txt_Callback(hObject, eventdata, handles)
% hObject    handle to dist_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dist_txt as text
%        str2double(get(hObject,'String')) returns contents of dist_txt as a double


% --- Executes during object creation, after setting all properties.
function dist_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dist_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Clear previous drawn points
function clear_points()
try
    while 1
        last_p = get(gca, 'children');
        if strcmp(last_p(1).Tag,'deleteme')
            delete(last_p(1))
        else
            break
        end
    end
end


% --- Executes on button press in slope.
function slope_Callback(hObject, eventdata, handles)
% hObject    handle to slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text2,'string','Angle')
hax = handles.plot;
clear_points
x = [0 0];
y = [0 0];
for n = 1:2
    [x(n),y(n)] = ginput(1);
    hold(hax,'on')
    plot(x(n),y(n),'+k')
end
plot(x,y,'r')
last_p = get(gca, 'children');
set(last_p(1:3),'Tag','deleteme')
angle = atand(abs((y(1)-y(2))/(x(1)-x(2))));
set(handles.dist_txt,'string',num2str(angle))
