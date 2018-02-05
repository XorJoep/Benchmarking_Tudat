function varargout = SFA_Orbit_Models(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SFA_Orbit_Models_OpeningFcn, ...
                   'gui_OutputFcn',  @SFA_Orbit_Models_OutputFcn, ...
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

% --- Executes just before SFA_Orbit_Models is made visible.
function SFA_Orbit_Models_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SFA_Orbit_Models (see VARARGIN)

% Choose default command line output for SFA_Orbit_Models
global stop;
stop = false;
handles.output = hObject;
handles.prop_buttons = [handles.rb_prop1; handles.rb_prop2; handles.rb_prop3; handles.rb_prop4];
handles.intg_buttons = [handles.rb_int1; handles.rb_int2; handles.rb_int3; handles.rb_int4];
% This sets up the initial plot - only do when we are invisible
% so window can get raised using SFA_Orbit_Models.
% Update handles structure
guidata(hObject, handles);
grid on
rotate3d on
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

function stop_button_Callback(hObject, eventdata, handles)
global stop;
stop = true;
set(hObject, 'String', 'Draw!', 'Enable', 'on', 'CallBack', @(hObject,eventdata)SFA_Orbit_Models('draw_button_Callback',hObject,eventdata,guidata(hObject)))

% --- Executes on button press in draw_button.
function draw_button_Callback(hObject, eventdata, handles)
% hObject    handle to draw_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stop
stop = false;
%% Get propagator and integrator combinations
set(hObject, 'String', 'Busy..', 'Enable', 'off')
prop_but = cell2mat(get(handles.prop_buttons, 'Value'));
intg_but = cell2mat(get(handles.intg_buttons, 'Value'));
prop_str = ["cowell","encke","gaussKeplerian","gaussModifiedEquinoctial"];
intg_str = ["Runge-Kutta78","Bulirsch-Stoer", "Runge-Kutta4","Euler"];
props = prop_str(prop_but==1);
intgs = intg_str(intg_but==1);
[A,B] = meshgrid(props,intgs);
c=cat(2,A',B');
combinations=reshape(c,[],2);
str_combs = combinations(:,1)+"_"+combinations(:,2);
handles.ISSD = false;
%% Get case
switch handles.Case_but_grp.SelectedObject.Tag
    case "rb_c0"
        str_combs = "MultiCase_keplerorbit_satellite_Out\Case_0_" + str_combs + "_E*.dat";
    case "rb_c1"
        str_combs = "MultiCase_keplerorbit_satellite_Out\Case_1_" + str_combs + "_E*.dat";
    case "rb_c2"
        str_combs = "MultiCase_capsule_reentry_Out\Case_0_" + str_combs + "_E*.dat";
    case "rb_ISSD"
        str_combs = "MultiCase_ISSD_Out\*_" + str_combs + "_E*.dat";
        handles.ISSD = true;
    otherwise
        disp("Error: case not found");
        return
end
acc = 0;
switch handles.acc_but_grp.SelectedObject.Tag
    case "rb_high"
        acc = 1;
    case "rb_med"
        acc = 2;
    case "rb_low"
        acc = 3;
    otherwise
        disp("unimplemented accuracy")
        return
end

guidata(hObject, handles);

%% Get files
dirs = char(str_combs);
for i = 1:size(dirs,1)
    files(i,:) = dir(dirs(i,:));
end
if(size(files(),2) == 0)
    disp("No files found for prop/integrator/accuracy combo");
    set(hObject, 'String', 'Draw!', 'Enable', 'on', 'CallBack', @(hObject,eventdata)SFA_Orbit_Models('draw_button_Callback',hObject,eventdata,guidata(hObject)))
    return;
end
files = struct2cell(files);
paths = files(1,:);
folder = files(2,1);

disp('loading...')
for i = 1:length(paths)
    data(i,:,:) = {importdata(strcat(folder{1},"\",paths{i}))};
end
disp('drawing...')

%% Plotting the data
cla reset
ax = handles.axes1;

t0 = handles.time_step_slid.Value;
L = 20; % length of comet tail
set(hObject, 'String', 'Stop!', 'Enable', 'on', 'CallBack', @(hObject,eventdata)SFA_Orbit_Models('stop_button_Callback',hObject,eventdata,guidata(hObject)))

timer = tic;
if(handles.ISSD)
    
    [x1,x2,x3] = sphere();
    r_sun = 695.7e6;
    surf(x1*r_sun, x2*r_sun, x3*r_sun, 'FaceAlpha', 0.1);
    
    radius = 2.5e11*handles.zoom_slid.Value;
    axis(ax,[-radius radius -radius radius -radius radius])
    disp('ISSD')
    for planet = 1:length(data)/3
        planet_data(planet, :, :) = data{(planet-1)*3+acc}(:,1:4);
        bodies(planet, 1) = line('parent',ax, 'color', 'red', 'marker','o', ...
            'xdata',planet_data(planet, 1, 2),'ydata',planet_data(planet, 1, 3),'zdata',planet_data(planet, 1, 4));
        bodies(planet, 2) = animatedline('parent',ax, 'color', 'red','linestyle','-',...
                        'MaximumNumPoints',L, 'LineWidth', 2);
        bodies(planet, 3) = animatedline('parent',ax,'linestyle','-.');
    end
    T = planet_data(1, :, 1); % zelfde voor elke planet
    
    for t = 1:length(T)
        for planet = 1:length(data)/3
            addpoints(bodies(planet,2), planet_data(planet,t,2), planet_data(planet,t,3), planet_data(planet,t,4));
            set(bodies(planet,1), 'xdata', planet_data(planet,t,2), 'ydata', planet_data(planet,t,3), 'zdata', planet_data(planet,t,4));
            if(t>L) addpoints(bodies(planet,3), planet_data(planet,t,2), planet_data(planet,t,3), planet_data(planet,t,4)); end
            if (stop) return; end
        end
        if T(t)>t0
            drawnow
            t0 = T(t) + handles.time_step_slid.Value;
            handles.time_step.String = strcat('t = ', num2str(T(t),'%.3e'), ' s');
            while toc(timer) < handles.real_time_slid.Value
            end
            timer = tic;
        end
    end
    
else

    T = data{acc}(:,1);
    X = data{acc}(:,2);
    Y = data{acc}(:,3);
    Z = data{acc}(:,4);

    [x1,x2,x3] = sphere();
    r_earth = 6.371e6;
    surf(x1*r_earth, x2*r_earth, x3*r_earth, 'FaceAlpha', 0.1);
    radius = 1e7*handles.zoom_slid.Value;
    axis(ax,[-radius radius -radius radius -radius radius])

    head = line('parent',ax, 'color', 'blue', 'marker','o', ...
        'xdata',X(1),'ydata',Y(1),'zdata',Z(1),'tag','head');
    tail = animatedline('parent',ax,'color','red','linestyle','-',...
                        'MaximumNumPoints',L, 'Tag','tail', 'LineWidth', 2);
    body = animatedline(ax);
    uistack(tail, 'top')
    for t = 1:length(T)
        addpoints(tail, X(t), Y(t), Z(t));
        if(t>L) addpoints(body, X(t-L), Y(t-L), Z(t-L)); end
        if (stop) return; end
        if T(t)>t0
            set(head,'xdata',X(t),'ydata',Y(t),'zdata',Z(t))
            drawnow
            t0 = T(t) + handles.time_step_slid.Value;
            handles.time_step.String = strcat('t = ', num2str(T(t),'%.3e'), ' s');
            while toc(timer) < handles.real_time_slid.Value
            end
            timer = tic;
        end
    end
end
set(hObject, 'String', 'Draw!', 'Enable', 'on', 'CallBack', @(hObject,eventdata)SFA_Orbit_Models('draw_button_Callback',hObject,eventdata,guidata(hObject)))

function varargout = SFA_Orbit_Models_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --- Executes on slider movement.
function time_step_slid_Callback(hObject, eventdata, handles)
% hObject    handle to time_step_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(get(hObject,'Value'),-1);
set(hObject, 'Value', val)
handles.time_step_text.String = strcat('plot time step: ', num2str(val), ' s');
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes on slider movement.
function real_time_slid_Callback(hObject, eventdata, handles)
% hObject    handle to real_time_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(get(hObject,'Value'),1);
set(hObject, 'Value', val)
handles.real_time_text.String = strcat('real time step: ', num2str(val), ' s');


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes on slider movement.
function zoom_slid_Callback(hObject, eventdata, handles)
ax = handles.axes1;
hObject.Value = round(hObject.Value,1);
if handles.ISSD
    radius = 2.5e11*hObject.Value;
else
    radius = 1e7*hObject.Value;
end
handles.zoom_text.String = strcat('zoom: ', num2str(radius,'%.1e'), 'x');
axis(ax,[-radius radius -radius radius -radius radius])
% hObject    handle to zoom_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
