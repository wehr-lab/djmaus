%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FILE NAME: Video_Player
%   DESCRIPTION: Video Player
%
%   Name          Date          Reason
%   Sahariyaz     28-May-2015   Basic GUI concepts 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = Video_Player(varargin)
% VIDEO_PLAYER MATLAB code for Video_Player.fig
%      VIDEO_PLAYER, by itself, creates a new VIDEO_PLAYER or raises the existing
%      singleton*.
%
%      H = VIDEO_PLAYER returns the handle to a new VIDEO_PLAYER or the handle to
%      the existing singleton*.
%
%      VIDEO_PLAYER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIDEO_PLAYER.M with the given input arguments.
%
%      VIDEO_PLAYER('Property','Value',...) creates a new VIDEO_PLAYER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Video_Player_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Video_Player_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Video_Player

% Last Modified by GUIDE v2.5 03-Jul-2018 17:00:00

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Video_Player_OpeningFcn, ...
                   'gui_OutputFcn',  @Video_Player_OutputFcn, ...
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


% --- Executes just before Video_Player is made visible.
function Video_Player_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Video_Player (see VARARGIN)

% Choose default command line output for Video_Player

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.pushbutton4,'Enable','on');

% UIWAIT makes Video_Player wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Video_Player_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ************************ BROWSE *****************************************
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%datapath=    '/Volumes/C/Users/lab/Desktop/Legless crickets/Cage.2/RT/4.3';
%datapath=     'C:\Users\lab\Desktop\Legless crickets';
%datapath=     'D:\lab\Data\CamMice\8287\CamDay3'
datapath= pwd;

cd(datapath)
[ video_file_name,video_file_path ] = uigetfile({'*.avi'},'Pick a video file');      %;*.png;*.yuv;*.bmp;*.tif'},'Pick a file');
if(video_file_path == 0)
    return;
end
input_video_file = [video_file_path,video_file_name];
set(handles.edit1,'String',input_video_file);
handles.video_file_path=video_file_path;
% Acquiring video
videoObject = VideoReader(input_video_file);
% Display first frame
frame_1 = read(videoObject,1);
axes(handles.axes1);
imshow(frame_1);
drawnow;
axis(handles.axes1,'off');
% Display Frame Number
set(handles.text3,'String','1');
set(handles.text4,'String',['  /  ',num2str(videoObject.NumberOfFrames)]);
set(handles.text2,'Visible','on');
set(handles.text3,'Visible','on');
set(handles.text4,'Visible','on');
%set(handles.pushbutton2,'Enable','on');
set(handles.pushbutton1,'Enable','off');
set(handles.pushbutton3,'Enable','on'); %pause/play
handles.frameCount=1;
handles.Running=1;

%Update handles
handles.videoObject = videoObject;
guidata(hObject,handles);


% ************************ BROWSE TEXT BOX ********************************
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




% ************************ PAUSE / PLAY ***********************************
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
videoObject = handles.videoObject;
axes(handles.axes1);
Running=handles.Running;

if(strcmp(get(handles.pushbutton3,'String'),'Pause'))
    set(handles.pushbutton3,'String','Play');
    %    uiwait();
    Running=0;
    handles.Running=0;
    guidata(hObject,handles);
else
    set(handles.pushbutton3,'String','Pause');
    %uiresume();
    %play
    Running=1;
    while Running
        frameCount=handles.frameCount;
        if videoObject.CurrentTime < videoObject.Duration
        frameCount=frameCount+1;
        set(handles.text3,'String',num2str(frameCount));
        frame = read(videoObject,frameCount);
        imshow(frame);
        handles.frameCount=frameCount;
        drawnow;
        guidata(hObject,handles);
        %handles = guidata(hObject);
        %fprintf('\n%s', get(handles.pushbutton3,'String'))
        %guidata(hObject,handles);
        %pushbutton2_Callback(hObject, eventdata, handles)
        if(strcmp(get(handles.pushbutton3,'String'),'Play'))
            Running=0;
        end
        else
            Running=0;
            set(handles.pushbutton3,'String','Play');
    end
end
end

% ************************ EXIT *******************************************
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);


% --- Executes on button press in Forward.
function Forward_Callback(hObject, eventdata, handles)
% hObject    handle to Forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

videoObject = handles.videoObject;
frameCount=handles.frameCount;
frameCount=frameCount+1;
axes(handles.axes1);

    % Display frames
    set(handles.text3,'String',num2str(frameCount));
    frame = read(videoObject,frameCount);
    imshow(frame);
    drawnow;
    handles.frameCount=frameCount;
guidata(hObject,handles);



% --- Executes on button press in Backward.
function Backward_Callback(hObject, eventdata, handles)
% hObject    handle to Backward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
videoObject = handles.videoObject;
frameCount=handles.frameCount;
frameCount=frameCount-1;
axes(handles.axes1);

    % Display frames
    set(handles.text3,'String',num2str(frameCount));
    frame = read(videoObject,frameCount);
    imshow(frame);
    drawnow;
    handles.frameCount=frameCount;
guidata(hObject,handles);


% --- Executes on button press in Check_cricket_speed.
function Check_cricket_speed_Callback(hObject, eventdata, handles)
% hObject    handle to Check_cricket_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

input_video_file=get(handles.edit1,'String');
if strfind(input_video_file, 'raw')
    datafile=strrep(input_video_file, 'raw', 'data');
elseif  strfind(input_video_file, 'MT')
    datafile=strrep(input_video_file, 'MT', 'data');
else error('data file not found')
end
out=LoadBonsaiTracks(datafile);
cricketxy=out.cricketxy;
%cricketxy=cricketxy(start_frame:stop_frame,:);
cspeed=sqrt(diff(cricketxy(:,1)).^2 + diff(cricketxy(:,2)).^2);
figure
plot(cspeed)
xlabel('frames')
ylabel('cricket speed, px/s')
title('cricket speed')



% --- Executes on button press in ReversePlay.
function ReversePlay_Callback(hObject, eventdata, handles)
% hObject    handle to ReversePlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
videoObject = handles.videoObject;
axes(handles.axes1);
Running=handles.Running;

if(strcmp(get(handles.ReversePlay,'String'),'Pause'))
    set(handles.ReversePlay,'String','Reverse Play');
    %    uiwait();
    Running=0;
    handles.Running=0;
    guidata(hObject,handles);
else
    set(handles.ReversePlay,'String','Pause');
    %uiresume();
    %play
    Running=1;
    while Running
        frameCount=handles.frameCount;
        frameCount=frameCount-1;
        if frameCount>1
            set(handles.text3,'String',num2str(frameCount));
            frame = read(videoObject,frameCount);
            imshow(frame);
            handles.frameCount=frameCount;
            drawnow;
            guidata(hObject,handles);
            %handles = guidata(hObject);
            %fprintf('\n%s', get(handles.pushbutton3,'String'))
            %guidata(hObject,handles);
            %pushbutton2_Callback(hObject, eventdata, handles)
            if(strcmp(get(handles.ReversePlay,'String'),'Reverse Play'))
                Running=0;
            end
        else
            Running=0;
            set(handles.ReversePlay,'String','Reverse Play');
        end
    end
end



function JumpToFrame_Callback(hObject, eventdata, handles)
% hObject    handle to JumpToFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of JumpToFrame as text
%        str2double(get(hObject,'String')) returns contents of JumpToFrame as a double
videoObject = handles.videoObject;
frameCount=str2num(get(hObject,'String'));
axes(handles.axes1);

    % Display frames
    set(handles.text3,'String',num2str(frameCount));
    frame = read(videoObject,frameCount);
    imshow(frame);
    drawnow;
    handles.frameCount=frameCount;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function JumpToFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to JumpToFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%ira 6.28.18 

% --- fixing missing or poorly fit pupil estimates
function fixed_pupil = pupilfix_Callback(hObject, eventdata, handles)
% hObject    handle to pupilfix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



name=['fixed_pupil_fit_',handles.videoObject.name, '.mat'];
try 
    load(name); % load previous frames
catch
    fixed_pupil1=[];
end
points=ginput;

fitted_points = fit_ellipse(points(:,1), points(:,2), handles.axes1); %fitted pupil from gived points (should be at least five)
fixed_pupil = [handles.frameCount fitted_points.X0_in fitted_points.Y0_in fitted_points.long_axis/2 fitted_points.short_axis/2 fitted_points.phi]; %get info from struct to matrix


fixed_pupil1=[fixed_pupil1; fixed_pupil]; %add new frame fit [ frame# x y long short]

save(name,'fixed_pupil1') %save (there should be a better way of accumulating new fits and saving them but i have not thought of it



% --- Executes on button press in skipnext.
function skipnext_Callback(hObject, eventdata, handles)
% hObject    handle to skipnext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
videoObject = handles.videoObject;
skipFrames=10;
axes(handles.axes1);

currentFrame=handles.frameCount;
frameCount= currentFrame+skipFrames;

% Display frames
set(handles.text3,'String',num2str(frameCount));
frame = read(videoObject,frameCount);
imshow(frame);
drawnow; hold on
handles.frameCount=frameCount;
guidata(hObject,handles);
try
name=['fixed_pupil_fit_',handles.videoObject.name, '.mat'];
load(name);
L=length(fixed_pupil1);
ellipse(fixed_pupil1(L,4), fixed_pupil1(L,5), fixed_pupil1(L,6), fixed_pupil1(L,2),fixed_pupil1(L,3), 'b')
drawnow
hold off
end


% --- Executes on button press in reapplyPupilFit.
function reapplyPupilFit_Callback(hObject, eventdata, handles)
% hObject    handle to reapplyPupilFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
name=['fixed_pupil_fit_',handles.videoObject.name, '.mat'];

try 
    load(name); % load previous frames
catch
    fixed_pupil1=[];
end
fixed_pupil=[];
last_fit_frame=fixed_pupil1(end,1);
num_frames=handles.frameCount-last_fit_frame;
for i=1:num_frames
    fixed_pupil=[fixed_pupil; last_fit_frame+1 fixed_pupil1(end, 2:end)];
    last_fit_frame=fixed_pupil(end,1);
end
fixed_pupil1=[fixed_pupil1; fixed_pupil];
save(name,'fixed_pupil1')
