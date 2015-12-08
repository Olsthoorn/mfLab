function varargout = showLayers(varargin)
%SHOWLAYERS a GUI for showing and observing 3D layers with wells and mesh
%
% Example:
%    h=showLayers(gr,'H.values',well,location);            --3D layers
%    h=showLayers(gr,'HK',well,location,'contours'); --flat layers
%
% Input:
%     H.values is a string with the name of the variable
%        It may be replaced by any 3D array or stress.
%        For example you may fill in:
%             'H.values'
%             'B(3).term{2}'
%             'HK' 'RIV'  'HRIV' 'CRIV' LCRIV', 'DRN' 'WEL'
%        The L stress prefix indicates log10,
%        The H prefix denotes heads,
%        The C prefix denotes conductance.
%        In fact, H is the 5th column of the stress and C the 6th column.
%        So HWEL would yield the flow.
%        Stresses without a prefix are equivalent to the head.
%
% Once the figure is on the screen, other 3D variables including stresses can be shown
% by changing the string in the edit box.
% Implemented for 7 layers (= number of model layers in the national hydrological instrumeent
% abbreviated to NHI and available on www.NHI.nu.
%
% See also plotXSec mf_wirefram mf_3Dblock gridObj
%
% TO 120605
%
% SHOWLAYERS MATLAB code for showLayers.fig
%      SHOWLAYERS, by itself, creates a new SHOWLAYERS or raises the existing
%      singleton*.
%
%      H = SHOWLAYERS returns the handle to a new SHOWLAYERS or the handle to
%      the existing singleton*.
%
%      SHOWLAYERS('CALLBACK',hObject,~,handles,...) calls the local
%      function named CALLBACK in SHOWLAYERS.M with the given input arguments.
%
%      SHOWLAYERS('Property','Value',...) creates a new SHOWLAYERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before showLayers_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to showLayers_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help showLayers

% Last Modified by GUIDE v2.5 03-Jun-2012 22:44:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;  % TO 120604 <-- no Singleton
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showLayers_OpeningFcn, ...
                   'gui_OutputFcn',  @showLayers_OutputFcn, ...
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

% --- Executes just before showLayers is made visible.
function showLayers_OpeningFcn(hObject,hevents, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to showLayers (see VARARGIN)

% Choose default command line output for showLayers
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% % See what variables may be shown
% load('name'); load(basename);
% varnames= who; varnames = {'H.values' 'B.term{1}','C.values',varnames{:}};
% for i=length(varnames):-1:1
%    switch varnames{i}
%        case {'GREP','BACKGROUND','AXIAL','basename','gr','LAYCBD'}
%            varnames(i)='';
%    end
% end
% display(varnames)

% This sets up the initial plot - only do when we are invisible
% so window can get raised using showLayers.
%showLayers(gr,H.values,well,location);

load('name'); load(basename);

xlabel('x [m]'); ylabel('y [m]'); zlabel('z[m +MSL]');
view(3);

k=0;

k=k+1;
if isempty(varargin) || ~isa(varargin{k},'gridObj')
        error('%s: first argument must be of class<<gridObj>>',mfilename);
else
    gr=varargin{k};
    k=k+1;
end

if length(varargin)<k || ~ischar(varargin{k})
    warning(['%s: second argument must be of type char referring to the name of a variable\n',...
        'will use default ''H(end).values'''], mfiles);
    varargin{k}='H(end).values';
end
handles.variable=varargin{k};
try
    Values = eval(handles.variable);
catch %#ok<CTCH>
    warning('mflab:showLayers:noInitialValues',...
        '%s Can''t find values for your initial variable <<%s>>, check the call to %s',...
        mfilename,handles.variable,mfilename);
    Values = input(sprintf('How to load this variable <<%s>>??',handles.variable));
%    Values = eval(handles.variable);
end
    
k=k+1;
if length(varargin)<k || ~isa(varargin{k},'wellObj') || ~isa(varargin{k},'MNW1Obj')
    well=[];
else
    well=varargin{k};
    k=k+1;
end

if length(varargin)<k || ~isa(varargin{k},'char')
    location='unkknown';
else
    location=varargin{k};
    k=k+1;
end

if length(varargin)>=k
    varargin=varargin(k:end);
else
    varargin=[];
end

%% if varargin contains the string 'con*'for contours of any kind
% TO DO: Make this surf +contour draped over it
% TO 120604
grey = get(gcf,'color');

if strmatchi('con',varargin)
    ic=strmatchi('con',varargin);
    if length(varargin) ==1, varargin=[];
    elseif length(varargin)==-ic,varargin=varargin(1:ic-1);
    else varargin = varargin([1:ic-1 ic+1:end]);
    end
    gr = gridObj(gr.xGr,gr.yGr,size(gr.zGr,3):-1:1,gr.LAYCBD,gr.MINDZ,gr.AXIAL);

    if isempty(varargin)
        handles.hLayer = gr.plotLayers(gca,1:gr.Nlay,Values,'edgecolor','none');
        handles.hMesh  = gr.plotMesh  (gca,'edgecolor',grey);
    else
        handles.hLayer = gr.plotLayers(gca,1:gr.Nlay,Values,varargin);
        handles.hMesh  = gr.plotMesh  (gca,'edgecolor',grey,varargin);
    end
else
    if isempty(varargin)
        handles.hLayer = gr.plotLayers(gca,1:gr.Nlay,Values,'edgecolor','none');
        handles.hMesh  = gr.plotMesh(gca,'edgecolor',grey);
    else
        handles.hLayer = gr.plotLayers(gca,1:gr.Nlay,Values,varargin);
        handles.hMesh  = gr.plotMesh(gca,'edgecolor',grey,varargin);
    end
end

if isempty(well)
    handles.hWell=[];
else
    handles.hWell=well.plot3D(gca);
end

handles.basename = basename;
handles.ttl      = sprintf('%%s for location "%s"',location);
handles.ht       = title(sprintf(handles.ttl,handles.variable));

guidata(hObject,handles);

%%
hObject = findobj('tag','edit1');
hObject = hObject(1); %###
set(hObject,'string',handles.variable);
edit1_Callback(hObject,hevents,handles)

% UIWAIT makes showLayers wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = showLayers_OutputFcn(hObject, ~, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, ~, handles)  %#ok
% hObject    handle to FileMenu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, ~, handles) %#ok
% hObject    handle to OpenMenuItem (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, ~, handles) %#ok
% hObject    handle to PrintMenuItem (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, ~, handles) %#ok
% hObject    handle to CloseMenuItem (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, ~, handles) %#ok
% hObject    handle to popupmenu2 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
switch get(hObject,'Value')
    case 1, set(handles.hWell,'color',[0.8 0.8 0.8]);
    case 2, set(handles.hWell,'color','b');
    case 3, set(handles.hWell,'color','k');
    case 4, set(handles.hWell,'color','r');
    case 5, set(handles.hWell,'color','g');
    case 6, set(handles.hWell,'color','m');
    case 7, set(handles.hWell,'color','c');
    case 8, set(handles.hWell,'color','y');
    case 9, set(handles.hWell,'color','b');
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, ~, handles) %#ok
% hObject    handle to popupmenu2 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, ~, handles) %#ok
% hObject    handle to checkbox15 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15
iLay=1;
if get(hObject,'value')
    set(handles.hLayer(iLay),'visible','on');
else
    set(handles.hLayer(iLay),'visible','off');
end


% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, ~, handles) %#ok
% hObject    handle to checkbox16 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16
iLay=2;
if get(hObject,'value')
    set(handles.hLayer(iLay),'visible','on');
else
    set(handles.hLayer(iLay),'visible','off');
end


% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, ~, handles)  %#ok
% hObject    handle to checkbox17 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17
iLay=3;
if get(hObject,'value')
    set(handles.hLayer(iLay),'visible','on');
else
    set(handles.hLayer(iLay),'visible','off');
end

% --- Executes on button press in checkbox18.
function checkbox18_Callback(hObject, ~, handles)  %#ok
% hObject    handle to checkbox18 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox18
iLay=4;
if get(hObject,'value')
    set(handles.hLayer(iLay),'visible','on');
else
    set(handles.hLayer(iLay),'visible','off');
end


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, ~, handles)  %#ok
% hObject    handle to checkbox19 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19
iLay=5;
if get(hObject,'value')
    set(handles.hLayer(iLay),'visible','on');
else
    set(handles.hLayer(iLay),'visible','off');
end

% --- Executes on button press in checkbox20.
function checkbox20_Callback(hObject, ~, handles) %#ok
% hObject    handle to checkbox20 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox20
iLay=6;
if get(hObject,'value')
    set(handles.hLayer(iLay),'visible','on');
else
    set(handles.hLayer(iLay),'visible','off');
end


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, ~, handles)  %#ok
% hObject    handle to checkbox21 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21
iLay=7;
if get(hObject,'value')
    set(handles.hLayer(iLay),'visible','on');
else
    set(handles.hLayer(iLay),'visible','off');
end


% --- Executes on button press in checkbox29.
function checkbox29_Callback(hObject, ~, handles)  %#ok
% hObject    handle to checkbox29 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox29
iFace=1;
if get(hObject,'Value')
    set(handles.hMesh(iFace),'visible','on');
else
    set(handles.hMesh(iFace),'visible','off');
end

% --- Executes on button press in checkbox30.
function checkbox30_Callback(hObject, ~, handles)  %#ok
% hObject    handle to checkbox30 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox30
iFace=2;
if get(hObject,'Value')
    set(handles.hMesh(iFace),'visible','on');
else
    set(handles.hMesh(iFace),'visible','off');
end


% --- Executes on button press in checkbox31.
function checkbox31_Callback(hObject, ~, handles)  %#ok
% hObject    handle to checkbox31 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox31
iFace=3;
if get(hObject,'Value')
    set(handles.hMesh(iFace),'visible','on');
else
    set(handles.hMesh(iFace),'visible','off');
end

% --- Executes on button press in checkbox32.
function checkbox32_Callback(hObject, ~, handles)  %#ok
% hObject    handle to checkbox32 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox32
iFace=4;
if get(hObject,'Value')
    set(handles.hMesh(iFace),'visible','on');
else
    set(handles.hMesh(iFace),'visible','off');
end

% --- Executes on button press in checkbox33.
function checkbox33_Callback(hObject, ~, handles)  %#ok
% hObject    handle to checkbox33 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox33
iFace=6;
if get(hObject,'Value')
    set(handles.hMesh(iFace),'visible','on');
else
    set(handles.hMesh(iFace),'visible','off');
end

% --- Executes on button press in checkbox34.
function checkbox34_Callback(hObject, ~, handles)  %#ok
% hObject    handle to checkbox34 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox34
iFace=5;
if get(hObject,'Value')
    set(handles.hMesh(iFace),'visible','on');
else
    set(handles.hMesh(iFace),'visible','off');
end

function edit1_Callback(hObject,~, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

load('name'); load(basename);

variable = get(hObject,'String');

if iscell(variable),variable = variable{1}; end

varname  = regexp(variable,'^gr\..+','match');
if isempty(varname)
    varname  = regexp(variable,'^[A-Za-z0-9]+','match');
end

%%
% Not in load, must be read separately

varname = upper(varname{1});

switch varname
    case 'H',    H=maskH(readDat([basename '.HDS']));  %#ok
    case 'C',    C=maskC(readMT3D('MT3D0001.UCN'));    %#ok
    case 'B',    B=readBud([basename,'.BGT']);  %#ok
    otherwise
        if exist('RIV','var') && ~isempty(strfind(varname,'RIV')), [Values,varname] = gr.getBCN(varname,RIV); end
        if exist('GBB','var') && ~isempty(strfind(varname,'GHB')), [Values,varname] = gr.getBCN(varname,GHB); end
        if exist('DRN','var') && ~isempty(strfind(varname,'DRN')), [Values,varname] = gr.getBCN(varname,DRN); end
        if exist('DRT','var') && ~isempty(strfind(varname,'DRT')), [Values,varname] = gr.getBCN(varname,DRT); end
        if exist('WEL','var') && ~isempty(strfind(varname,'WEL')), [Values,varname] = gr.getBCN(varname,WEL); end
        if exist('STR','var') && ~isempty(strfind(varname,'STR')), [Values,varname] = gr.getBCN(varname,STR); end
        if exist('CHD','var') && ~isempty(strfind(varname,'CHD')), [Values,varname] = gr.getBCN(varname,CHD); end
        if exist('PNTSRC','var') && ~isempty(strfind(varname,'PNTSRC')), [Values,varname] = gr.getBCN(varname,PNTSRC); end
end

%%
% Get the data Verbatim, catch if doesn't work
try
    %% If parameter does not yield values, do nothing
    if ~exist('Values','var') || isempty(Values)
        try
            eval(['Values =' variable ';']);
        catch  %#ok
        warning('mfLab:showLayers:InputWithoutResult',...
            '%s: The input in your edit line does not yield a 3D parameter to show',mfilename);
            return;
        end
    else
        % We have values and they are numeric as they come from
        % RIV, GHB or DRN
    end
    
    if all(size(Values(:,:,1))==gr.size(1:2))

        %% Set the color limits for this variable
        range = ContourRange(Values,50);
        if isempty(range)
            m = mean(Values(:));
            if m~=0,
                range=m*[0.5 2];
            else
                range=m+[-0.5 0.5];
            end
        end
        caxis([min(range),max(range)]);

        %% Put varname in the title
        set(handles.ht,'string',sprintf(handles.ttl,varname));
        
        %% Put full variable in the edit of the object
        set(hObject,'string',variable);

        % store variable name
        handles.variable=varname;
        
                %% Update the Cdata for the layers (only the layers up to the size of 
        % either values or 7 (number of NHI layers)
        for iLay=1:min(7,size(Values,3))
            if ~isnan(Values(1,1,iLay))
                set(handles.hLayer(iLay),'Cdata',Values(:,:,iLay),'FaceColor','flat');
            else
                set(handles.hLayer(iLay),'FaceColor','none');
            end
        end
        for iLay=iLay+1:min(7,length(handles.hLayer))
            set(handles.hLayer(iLay),'FaceColor','none');
        end

        % store the handles in the figure
        guidata(get(get(hObject,'parent'),'parent'),handles);
    else
        warning('showLayers:showLayers:NoSuchVarInWorkspace',...
            '%s Can''t find variable var or var is not a 3D array of appropriate size',variable);
    end
catch %#ok
    guidata(get(get(hObject,'parent'),'parent'),handles);

    warning('showLayers:showLayers:NoSuchVarInWorkspace',...
        '%s Can''t find variable var or var is not a 3D array of appropriate size',variable);
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, hevents, handles) %#ok
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
