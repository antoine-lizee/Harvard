function varargout = Solutions_on_the_cone(varargin)
% SOLUTIONS_ON_THE_CONE MATLAB code for Solutions_on_the_cone.fig
% --------------------------------------------------------------------------  
%   This script is the script linked to the figure of the same name.
%   Running the script will launch the GUI. Enjoy !
%
%   Copyright Antoine Lizée July 2011 
%   @Dumais lab, OEB dpt, Harvard University.
%
%---------------------------------------------------------------------------
%     SOLUTIONS_ON_THE_CONE, by itself, creates a new SOLUTIONS_ON_THE_CONE or raises the existing
%      singleton*.
%
%      H = SOLUTIONS_ON_THE_CONE returns the handle to a new SOLUTIONS_ON_THE_CONE or the handle to
%      the existing singleton*.
%
%      SOLUTIONS_ON_THE_CONE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOLUTIONS_ON_THE_CONE.M with the given input arguments.
%
%      SOLUTIONS_ON_THE_CONE('Property','Value',...) creates a new SOLUTIONS_ON_THE_CONE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Solutions_on_the_cone_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Solutions_on_the_cone_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%   
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Solutions_on_the_cone

% Last Modified by GUIDE v2.5 06-Oct-2011 16:43:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Solutions_on_the_cone_OpeningFcn, ...
                   'gui_OutputFcn',  @Solutions_on_the_cone_OutputFcn, ...
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


% --- Executes just before Solutions_on_the_cone is made visible.
function Solutions_on_the_cone_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Solutions_on_the_cone (see VARARGIN)

set(hObject, 'CurrentAxes',handles.axes1);

%Initialisation of user parameters :
set(handles.check2, 'Value', 1); %Check the 2 edges checkbox
set(handles.check_display, 'Value', 1);
set(handles.check_daughter, 'Value', 1);

%Initialisation of some figure parameters :
title(handles.axes2,'Cell bounadaries on the 2D plan and after the transformation');

% Creation of the cone
[X Y]=meshgrid(-1:0.05:1,-1:0.05:1);
Z=X*0;
surf(X,Y,Z,'Facecolor',[0.3;0.7;0.1],'EdgeColor','none');

%3D Plot edition
lighting phong % Define the way of computing the lightning effects. 'phong' is an interpolation. Try flat to understand... This function set the FaceLighting and EdgeLighting properties of surfaces and patches 
l=light('Position',[-2,0,2]);
material([0.5, 0.9, 0.3]); %Set the reflectance properties of the surface
alpha(0.7); % Set the transparency
view([-27 24]);
set(gca,'ZTick',zeros(1,0), 'ZTickMode', 'manual','YTick',zeros(1,0),...
   'YTickMode', 'manual', 'XTick',zeros(1,0), 'XtickMode','manual',...
   'DataAspectRatio',[1 1 1], 'DataAspectRatioMode', 'manual',...
   'NextPlot','replacechildren');

% Remove Ticks for axes2
set(handles.axes2,'YTick',zeros(1,0),...
    'YTickMode', 'manual', 'XTick',zeros(1,0), 'XtickMode','manual',...
    'DataAspectRatio',[1 1 1], 'DataAspectRatioMode', 'manual',...
    'NextPlot','replacechildren');

% Choose default command line output for Solutions_on_the_cone
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Solutions_on_the_cone wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Solutions_on_the_cone_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in check2.
function check2_Callback(hObject, eventdata, handles)
% hObject    handle to check2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text_warning,'String','Change gamma to refresh the display','Visible','on');
% Hint: get(hObject,'Value') returns toggle state of check2


% --- Executes on button press in check3.
function check3_Callback(hObject, eventdata, handles)
% hObject    handle to check3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text_warning,'String','Change gamma to refresh the display','Visible','on');
% Hint: get(hObject,'Value') returns toggle state of check3


% --- Executes on button press in check4.
function check4_Callback(hObject, eventdata, handles)
% hObject    handle to check4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_warning,'String','Change gamma to refresh the display','Visible','on');
% Hint: get(hObject,'Value') returns toggle state of check4


% --- Executes on button press in check_display.
function check_display_Callback(hObject, eventdata, handles)
% hObject    handle to check_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')==0
    set(handles.axes2,'Visible','off');
    set(get(handles.axes2,'children'),'Visible','off');
else
    set(handles.axes2,'Visible','on');
    set(get(handles.axes2,'children'),'Visible','on');
end

% Hint: get(hObject,'Value') returns toggle state of check_display


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gamma=get(hObject,'Value');
if gamma*10~=floor(gamma*10)
    set(hObject,'Value',floor(gamma*10)/10);
    gamma=get(hObject,'Value');
end

I=[2 3 5 8]; %Index of total number of generations
G=I(get(handles.list_gen,'Value'))-1; %number of daughter generations

set(handles.textgamma, 'String', ['Gamma = ' num2str(gamma) ' deg']);


%%Drawing in the main plot

set(gcf,'currentaxes',handles.axes1);
currentview=get(gca, 'view');
cla
hold on

% Creation of the cone
a=(2*pi-gamma*pi/180)/(2*pi);
h=sqrt(1/a^2-1);
bd=1.3/sqrt(1+h);
[X Y]=meshgrid(-bd:0.05:bd,-bd:0.05:bd);
Rho=sqrt(X.^2+Y.^2);
d=0.05;
Z=-Rho*h-h*d*exp(-Rho/d);
a=Z(2,floor(size(Z,2)/2));
bool=Z<a;
Z(bool)=a;
surf(X,Y,Z,'Facecolor',[0.3;0.7;0.1],'EdgeColor','none');

%3D Plot edition
lighting phong % Define the way of computing the lightning effects. 'phong' is an interpolation. Try flat to understand... This function set the FaceLighting and EdgeLighting properties of surfaces and patches 
l=light('Position',[-2,0,2]);
material([0.5, 0.9, 0.3]); %Set the reflectance properties of the surface
alpha(0.7); % Set the transparency
view(currentview);

% Plot the cell boundaries
j=[];
for ii=2:5
    if eval(['get(handles.check' num2str(ii) ',''Value'')==1']) %testing the checkboxes
        gamma_max=[321.1 266.8 187.4 92.9];
        if gamma<=gamma_max(ii-1)
            [ V A, ~ ]=findfixedpoint(ii,gamma,1);
            [VF VFd]=draw_daughter_on_cone(V,A,1);
       
        else
            j=[j ii];
            continue;
        end
        ha(ii).m=plot3(VF(:,1),VF(:,2),VF(:,3),'LineWidth',2); %plotting the mother cell
        if get(handles.check_daughter,'Value')==1; %plotting the daughter cell
            ha(ii).d=plot3(VFd(:,1),VFd(:,2),VFd(:,3),'r','LineWidth',1.5);
            for T=2:G
              [VF VFd]=draw_daughter_on_cone(V,A,1,T);
              ha(ii).d=plot3(VFd(:,1),VFd(:,2),VFd(:,3),'k','LineWidth',1.5);
          end
        end
        
        %drawing the auxilliar plot
        set(gcf,'currentaxes', handles.axes2);
        cla;
        hold on;
        draw_2D(V,A);
        plot(VF(:,1),VF(:,2));
        axis equal
        set(gcf,'currentaxes',handles.axes1);        
        if get(handles.check_display,'Value')==0;
            set(get(handles.axes2,'children'),'Visible','off');
        end
    end
end
% for the type 2 solutions :
ii=5;
if eval(['get(handles.check' num2str(ii) '_2,''Value'')==1']) %testing the checkboxes
    gamma_max=[0 0 0 161.4];
    if gamma<=gamma_max(ii-1)
        [ V A, ~ ]=findfixedpoint(ii+0.2,gamma,1);
        [VF VFd]=draw_daughter_on_cone_2(V,A,1);
        ha(ii).m=plot3(VF(:,1),VF(:,2),VF(:,3),'LineWidth',2); %plotting the mother cell
        if get(handles.check_daughter,'Value')==1; %plotting the daughter cell
            ha(ii).d=plot3(VFd(:,1),VFd(:,2),VFd(:,3),'r','LineWidth',1.5);
            for T=2:G
              [VF VFd]=draw_daughter_on_cone_2(V,A,1,T);
              ha(ii).d=plot3(VFd(:,1),VFd(:,2),VFd(:,3),'k','LineWidth',1.5);
          end
        end
        %drawing the auxilliar plot
        set(gcf,'currentaxes', handles.axes2);
        cla;
        hold on;
        draw_2D(V,A);
        plot(VF(:,1),VF(:,2));
        axis equal
        set(gcf,'currentaxes',handles.axes1);        
        if get(handles.check_display,'Value')==0;
            set(get(handles.axes2,'children'),'Visible','off');
        end
    else
        j=[j ii+0.2];
    end
end

if ~isempty(j) %Display warning message       
    string=['the angle gamma of the cone is too big to display the solution f' num2str(j,'or %g ') ' edges'];
    set(handles.text_warning,'String',string,'Visible','on');
else
    set(handles.text_warning,'Visible','off');
end




% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9])
end
max=get(hObject,'Max');
step=1/max;
set(hObject,'Value',0,'SliderStep', [step/10, step]);


% --- Executes on button press in check_daughter.
function check_daughter_Callback(hObject, eventdata, handles)
% hObject    handle to check_daughter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')==0
    set(handles.list_gen,'Visible','off');
    set(handles.text_gen,'ForegroundColor', [0.8 0.8 0.8]);
else
    set(handles.list_gen,'Visible','on');
    set(handles.text_gen,'ForegroundColor', 'k');
end
set(handles.text_warning,'String','Change gamma to refresh the display','Visible','on');
% Hint: get(hObject,'Value') returns toggle state of check_daughter


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textgamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textgamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in check5.
function check5_Callback(hObject, eventdata, handles)
% hObject    handle to check5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_warning,'String','Change gamma to refresh the display','Visible','on');
% Hint: get(hObject,'Value') returns toggle state of check5


% --- Executes on button press in check5_2.
function check5_2_Callback(hObject, eventdata, handles)
% hObject    handle to check5_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_warning,'String','Change gamma to refresh the display','Visible','on');
% Hint: get(hObject,'Value') returns toggle state of check5_2


% --- Executes when uipanel1 is resized.
function uipanel1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in list_gen.
function list_gen_Callback(hObject, eventdata, handles)
% hObject    handle to list_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_warning,'String','Change gamma to refresh the display','Visible','on');
% Hints: contents = cellstr(get(hObject,'String')) returns list_gen contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_gen


% --- Executes during object creation, after setting all properties.
function list_gen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
