function varargout = view_geom(varargin)
% VIEW_GEOM MATLAB code for view_geom.fig
%      VIEW_GEOM, by itself, creates a new VIEW_GEOM or raises the existing
%      singleton*.
%
%      H = VIEW_GEOM returns the handle to a new VIEW_GEOM or the handle to
%      the existing singleton*.
%
%      VIEW_GEOM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_GEOM.M with the given input arguments.
%
%      VIEW_GEOM('Property','Value',...) creates a new VIEW_GEOM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view_geom_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view_geom_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view_geom

% Last Modified by GUIDE v2.5 01-Nov-2018 16:39:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @view_geom_OpeningFcn, ...
                   'gui_OutputFcn',  @view_geom_OutputFcn, ...
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


% --- Executes just before view_geom is made visible.
function view_geom_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to view_geom (see VARARGIN)

% Choose default command line output for view_geom
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% LOAD INI FILE AND PRESETS
o = check_ini();
if o
    load_ini = load('ini.mat');
    dir_p = load_ini.dir_p;
    handles.folder.String = dir_p;
    show_files(dir_p,handles);
    if isfield(load_ini,'n1') && isfield(load_ini,'n2') && ...
            isfield(load_ini,'n3') && isfield(load_ini,'p1') && ...
            isfield(load_ini,'p2')
        handles.n1.String = load_ini.n1;
        handles.n2.String = load_ini.n2;
        handles.n3.String = load_ini.n3;
        handles.p1.String = load_ini.p1;
        handles.p2.String = load_ini.p2;
    end
end

% Variables to load
VAR = {'GEOMETRY','Geo_id','km'};
setappdata(0,'VAR',VAR)
axis_p = [];
setappdata(0,'axis_p',axis_p)

% Load bar
handles.loadbar.XLim = [0 100];
handles.loadbar.YLim = [0 1];
handles.loadbar.XTick = [];
handles.loadbar.YTick = [];
axes(handles.loadbar)
patch([0 100 100 0],[0 0 1 1],'w')
hold on


% UIWAIT makes view_geom wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = view_geom_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function folder_Callback(hObject, eventdata, handles)
% hObject    handle to folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of folder as text
%        str2double(get(hObject,'String')) returns contents of folder as a double


% --- Executes during object creation, after setting all properties.
function folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
o = check_ini;
if o
    load_ini = load('ini.mat');
    dir_p = load_ini.dir_p;
else
    dir_p = uigetdir('~/');
end
show_files(dir_p,handles);
set(handles.folder,'String',dir_p)
save('ini','dir_p')


% --- Executes on selection change in data.
function data_Callback(hObject, eventdata, handles)
% hObject    handle to data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data


% --- Executes during object creation, after setting all properties.
function data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n1_Callback(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n1 as text
%        str2double(get(hObject,'String')) returns contents of n1 as a double
o = check_ini;
n1 = handles.n1.String;
if o
    save('ini','n1','-append')
else
    save('ini','n1')
end



% --- Executes during object creation, after setting all properties.
function n1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p1_Callback(hObject, eventdata, handles)
% hObject    handle to p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p1 as text
%        str2double(get(hObject,'String')) returns contents of p1 as a double
o = check_ini;
p1 = handles.p1.String;
if o
    save('ini','p1','-append')
else
    save('ini','p1')
end


% --- Executes during object creation, after setting all properties.
function p1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n2_Callback(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n2 as text
%        str2double(get(hObject,'String')) returns contents of n2 as a double
o = check_ini;
n2 = handles.n2.String;
if o
    save('ini','n2','-append')
else
    save('ini','n2')
end


% --- Executes during object creation, after setting all properties.
function n2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p2_Callback(hObject, eventdata, handles)
% hObject    handle to p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p2 as text
%        str2double(get(hObject,'String')) returns contents of p2 as a double
o = check_ini;
p2 = handles.p2.String;
if o
    save('ini','p2','-append')
else
    save('ini','p2')
end


% --- Executes during object creation, after setting all properties.
function p2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n3_Callback(hObject, eventdata, handles)
% hObject    handle to n3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n3 as text
%        str2double(get(hObject,'String')) returns contents of n3 as a double
o = check_ini;
n3 = handles.n3.String;
if o
    save('ini','n3','-append')
else
    save('ini','n3')
end


% --- Executes during object creation, after setting all properties.
function n3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% Check if there is an initial file to load
function [o] = check_ini()
direct = dir('./');
o = false;
for n = 1:length(direct)
    if strcmp(direct(n).name,'ini.mat')
        o = true;
    end
end


% Show files of the folder in a listbox
function show_files(dir_p,handles)
filesINdir = dir(dir_p);
for n = 1:length(filesINdir)-2
    handles.files.String{n} = filesINdir(n+2).name;
end


% --- Executes on selection change in files.
function files_Callback(hObject, eventdata, handles)
% hObject    handle to files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns files contents as cell array
%        contents{get(hObject,'Value')} returns selected item from files


% --- Executes during object creation, after setting all properties.
function files_CreateFcn(hObject, eventdata, handles)
% hObject    handle to files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x = str2num(handles.p1.String);
y = str2num(handles.p2.String);
XX = repmat(x,length(y),1);
YY = repmat(y',1,length(x));

X = XX(:);
Y = YY(:);

dir_p = handles.folder.String;
n1 = handles.n1.String;
n2 = handles.n2.String;
n3 = handles.n3.String;

step = 10;

VAR = getappdata(0,'VAR');

axes(handles.loadbar)
c = patch([0 100 100 0],[0 0 1 1],'w');

for n = 1:length(X)
    dirn = [dir_p,'/',n1,num2str(X(n)),n2,num2str(Y(n)),n3];
    laststep = lastest(dirn);
    steps = 1:step:laststep;
    finalstep(n) = laststep;
    for t = 1:length(steps)
        data{n,t} = load([dirn,'/_',num2str(steps(t)),'.mat'],VAR{:});
        
        % Load bar
        perc = 100*(n-1+t/length(steps))/length(X);
        delete(c);
        title([num2str(perc),' %'])
        c = patch([0 perc perc 0],[0 0 1 1],'r');
        drawnow
    end
end
xd{1} = x;
yd{1} = y;
setappdata(0,'xd',xd)
setappdata(0,'yd',yd)
setappdata(0,'data',data);
setappdata(0,'finalstep',finalstep);
setappdata(0,'step',step);
handles.timebar.Max = max(finalstep);
handles.timebar.Value = str2double(handles.time.String);
handles.timebar.SliderStep = [step/max(finalstep) step/max(finalstep)*10];


% --- Executes on slider movement.
function timebar_Callback(hObject, eventdata, handles)
% hObject    handle to timebar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
finalstep = getappdata(0,'finalstep');
step = getappdata(0,'step');
list_step = 1:step:max(finalstep);
[~,ord] = sort(abs(list_step-handles.timebar.Value));
handles.timebar.Max = max(finalstep);
handles.timebar.Value = str2double(handles.time.String);
handles.timebar.SliderStep = [step/max(finalstep) step/max(finalstep)*10];
handles.timebar.Value = list_step(ord(1));
handles.time.String = num2str(handles.timebar.Value);

xd = getappdata(0,'xd');
yd = getappdata(0,'yd');
l = length(xd{1});
h = length(yd{1});
data = getappdata(0,'data');
for n = 1:l*h
    u(n) = subplot(h,l,n);
    u(n).Position = [u(n).Position(1)*0.7 u(n).Position(2) ...
        u(n).Position(3)*0.7 u(n).Position(4)*0.7];
    plot_geom(u(n),handles,data{n,find(handles.timebar.Value==list_step)})
end
linkaxes(u)


% --- Executes during object creation, after setting all properties.
function timebar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timebar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function time_Callback(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time as text
%        str2double(get(hObject,'String')) returns contents of time as a double
handles.timebar.Value = str2double(handles.time.String);


% --- Executes during object creation, after setting all properties.
function time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Plot geometries
function plot_geom(fig,handles,data)

GEOMETRY = data.GEOMETRY;
Geo_id = data.Geo_id;
km = data.km;

% Load variables from the correct format
if size(GEOMETRY,1) == 1
    geometry_p = GEOMETRY.bnd;
    Geo_id = GEOMETRY.id;
else
    geometry_p = GEOMETRY;
end

% Generate random colors for the different interfaces
color_r = rand(max(Geo_id),3);

% Loop to plot the different interfaces
for layer_plot = 1:max(Geo_id)
    plot(geometry_p(1,Geo_id==layer_plot)/km, ...
        geometry_p(2,Geo_id==layer_plot)/km,'-', ...
        'Color','k','LineWidth',1)
    hold on
end

% Headings
title('Geometry')
xlabel('Distance [km]')
ylabel('Depth [km]')

axis_p = getappdata(0,'axis_p');
if ~isempty(axis_p)
    axis(axis_p)
else
    axis_p = [fig.XLim fig.YLim];
    setappdata(0,'axis_p',axis_p)
end

hold off
