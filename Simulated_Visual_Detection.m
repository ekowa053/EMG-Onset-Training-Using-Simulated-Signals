function varargout = Simulated_Visual_Detection(varargin)
% SIMULATED_VISUAL_DETECTION MATLAB code for Simulated_Visual_Detection.fig
%      SIMULATED_VISUAL_DETECTION, by itself, creates a new SIMULATED_VISUAL_DETECTION or raises the existing
%      singleton*.
%
%      H = SIMULATED_VISUAL_DETECTION returns the handle to a new SIMULATED_VISUAL_DETECTION or the handle to
%      the existing singleton*.
%
%      SIMULATED_VISUAL_DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMULATED_VISUAL_DETECTION.M with the given input arguments.
%
%      SIMULATED_VISUAL_DETECTION('Property','Value',...) creates a new SIMULATED_VISUAL_DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Simulated_Visual_Detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Simulated_Visual_Detection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Simulated_Visual_Detection

% Last Modified by GUIDE v2.5 15-Aug-2016 13:45:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Simulated_Visual_Detection_OpeningFcn, ...
                   'gui_OutputFcn',  @Simulated_Visual_Detection_OutputFcn, ...
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


% --- Executes just before Simulated_Visual_Detection is made visible.
function Simulated_Visual_Detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Simulated_Visual_Detection (see VARARGIN)

% Only for a readily created sequence of Signals and Onset times
% Signals - Each column is a new trial
% Onsets  - Each row is a trial, each column is the onset for a spike
a = load('Signals.mat');
b = load('Onsets.mat');
Signals  = a.Signals;
Onsets   = b.Onsets;

handles.Signals = Signals;
handles.Onsets  = Onsets;

handles.TimesRun = 1;

getDateTime  = cellstr(datestr(clock));          % Get the current date and time
currDateTime = strrep(getDateTime,':','_');      % Replace colon with underscore
currDateTime = strrep(currDateTime,'-','_');     % Replace minus sign with underscore
currDateTime = strrep(currDateTime,' ','_');     % Replace space with underscore
currDateTime = strcat(currDateTime,'.xlsx');     % Excel file extenstion
currDateTime = strcat('Visual_Onsets_',currDateTime);
currDateTime = strjoin(currDateTime);  

handles.save = currDateTime;

% Choose default command line output for Simulated_Visual_Detection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Simulated_Visual_Detection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Simulated_Visual_Detection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in VISUALpb.
function VISUALpb_Callback(hObject, eventdata, handles)
% hObject    handle to VISUALpb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Assuming a sampling rate of 1000
[filt_1,filt_2] = butter(2,[24 400]/500,'bandpass');

% Signals is a matrix of 'padadd' signals
% Onsets  is a matrix of the true onset times
Signals = handles.Signals;
Onsets  = handles.Onsets;
num_trials = min(size(Signals));

ord = randperm(num_trials);

reordered = Signals(:,ord);
handles.Reordered = reordered;

% num_events is the number of spikes per signal 
num_events = 5;
onset_v      = zeros(num_trials,num_events + 2);
onset_v(:,1) = ord';



axes(handles.axes1)
figure('units','normalized','outerposition',[0 0 1 1])
i = 1;
while i <= num_trials
    pts = reordered(:,i);               % Get random order
    pts = pts(isnan(pts)~=1);           % Remove NaNs
    pts = filtfilt(filt_1,filt_2,pts);  % Bandpass
    
    plot(pts);
    xlim([0 numel(pts)])
    ylim([-2 2])
    title(num2str(i))
    
    [x,~] = ginput(num_events);
    x = round(x);
    
%     disp(['Your points are: ', num2str(x')])
%     onset_v(i,3:end) = num2cell(x)
    onset_v(i,3:end) = x;

    
    cla
    i = i+1;
end
close


onset_v = sortrows(onset_v);

save_ind = handles.TimesRun;

times_v      = onset_v(:,3:end) - Onsets;
sav          = cell(num_trials,num_events+3);
sav(1)       = {save_ind};
sav(:,2)     = num2cell([1:num_trials]');
sav(:,4:end) = num2cell(times_v);


% Now do Bonato
% Then AGLR


assignin('base','times',sav)



xlswrite(handles.save,sav,save_ind,'A1')
handles.TimesRun = save_ind+1;

handles.times_v = times_v;

guidata(hObject,handles)


% % % uicontrol(handles.BONATOpb);
% % % BONATOpb_Callback(handles.BONATOpb,[],handles);
% % % 
% % % 
% % % 
% % % % --- Executes on button press in BONATOpb.
% % % function BONATOpb_Callback(hObject, eventdata, handles)
% % % % hObject    handle to BONATOpb (see GCBO)
% % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % handles    structure with handles and user data (see GUIDATA)
% % % Reordered = handles.Reordered;
% % % Onsets    = handles.Onsets;
% % % 
% % % 
% % % 
% % % 
% % % % Sampling Frequency
% % % SF = 1000;
% % % %1) Bandpass for TKE
% % % [f,e] = butter(3,[30 300]/(SF/2),'bandpass');
% % % %2) Lowpass for TKE
% % % [h,g] = butter(2,50/(SF/2),'low');
% % % %3)High Pass Filter to remove bias
% % % [h_1,h_2]=butter(2,20/ (SF/2), 'high');
% % % %4)Zero-lag Low Pass Butterworth Filter
% % % [l_1,l_2]=butter(2,6/ (SF/2), 'low');
% % % 
% % % 
% % % 
% % % for i =1:10
% % %     pts = Reordered(:,i);
% % %     pts = pts(isnan(pts)~=1);
% % %     assignin('base','pts',pts)
% % %     x = pts;
% % %     data = x;
% % %     disp(['Currently on: ',num2str(i)])
% % %     
% % %     
% % %     x2 = filtfilt(h_1,h_2,x);
% % %     x3 = abs(x2);
% % %     x4= filtfilt(l_1,l_2,x3);
% % %     x=x4;
% % % 
% % % 
% % % 
% % % 
% % %     % Get SNR, Get noise variance
% % %     inp = figure('units','normalized','outerposition',[0 0 1 1]);
% % %     plot(data)
% % %     [t_sn(1),y]=ginput(1);
% % %     t_sn(1)=floor(t_sn(1));
% % %     [t_sn(2),y]=ginput(1);
% % %     t_sn(2)=floor(t_sn(2));
% % %     [t_sn(3),y]=ginput(1);
% % %     t_sn(3)=floor(t_sn(3));
% % %     [t_sn(4),y]=ginput(1);
% % %     t_sn(4)=floor(t_sn(4));
% % %     n_sigma=sum(data(t_sn(1):t_sn(2)).^2)./(t_sn(2)-t_sn(1)); %Mean Variance of noise
% % %     s_sigma=sum(data(t_sn(3):t_sn(4)).^2)./(t_sn(4)-t_sn(3)); %Mean variance of signal
% % %     snr=10*log10(s_sigma/n_sigma); %Signal to Noise ratio (in decibels)?
% % % 
% % %     % Second threshold
% % % 
% % %     prob_false = 0.05;
% % %     r0 = 1;  % sucessive points
% % % 
% % %     switch r0
% % %         case 1
% % %             c=[1 -5 10 -10 5 -prob_false];
% % %         case 2
% % %             c=[-4 15 -20 10 0 -prob_false];
% % %         case 3
% % %             c=[6 -15 10 0 0 -prob_false];
% % %         case 4
% % %             c=[-4 5 0 0 0 -prob_false];
% % %         case 5
% % %             c=[1 0 0 0 0 -prob_false];
% % %         otherwise
% % %             c=[1 -5 10 -10 5 -prob_false];
% % %     end
% % %     rt=roots(c); %solving the polynomial equation
% % %     rtr=[];
% % %     for i=1:max(size(rt)) %
% % %         if isreal(rt(i))
% % %             rtr=[rtr rt(i)];
% % %         end
% % %     end
% % %     I=find(rtr>0 & rtr<1); %finding the maximum real root of the polynomial equation
% % %     p_zeta=rtr(I); % probability for the threshold
% % %     pdk=exp(log(p_zeta)./(1+10.^(snr/10))); % obsolete - probability that kth sample is above the threshold
% % %     pd=1-binocdf(r0-1,5,pdk);               % obsolete -  probability that signal samples, albeit corrupted by noise, are correctly recognized
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % %     % Detection
% % %     if mod(length(x),2)
% % %         x = x(1:end-1);
% % %     end
% % % 
% % % 
% % %     aux = x.^2; % auxiliary sequence
% % % 
% % %     if ~isrow(aux)
% % %        aux = aux' ;
% % %     end
% % % 
% % %     aux_seq = aux(1,1:2:end) + aux(1,2:2:end); % add pairs together
% % % 
% % % 
% % %     zeta=chi2inv(1-p_zeta,2)*n_sigma; % the first threshold value 
% % %     %zeta = -log(p_zeta)*2*n_sigma;   % this is equivalent to the above
% % % 
% % % 
% % %     dtct=zeros(1,max(size(aux_seq)));   % initialize vector
% % % 
% % %     % detection scanning (window of 5)
% % %     for i=1:5:max(size(aux_seq))-5
% % %         r=0;
% % %         for j=i:i+4 %detecting window's width is 5 points
% % %             if aux_seq(j)>zeta
% % %                 r=r+1;
% % %             end
% % %         end
% % %         if r>=r0;
% % %             dtct(i:i+4)=1;
% % %         end
% % %     end
% % % 
% % % 
% % %     % Find actual change points
% % %     % Convert dtct (series of 0s and 1s - indicating change) to actual pairs of
% % %     % on/off
% % %     test = find(dtct~=0);
% % %     problem = 0;
% % % 
% % %     if isempty(test)
% % %        disp('No points found') 
% % %     end
% % %     
% % % while isempty(test) || pd < 0.95
% % %     
% % %     plot(data)
% % %     [t_sn(1),y]=ginput(1);
% % %     t_sn(1)=floor(t_sn(1));
% % %     [t_sn(2),y]=ginput(1);
% % %     t_sn(2)=floor(t_sn(2));
% % %     [t_sn(3),y]=ginput(1);
% % %     t_sn(3)=floor(t_sn(3));
% % %     [t_sn(4),y]=ginput(1);
% % %     t_sn(4)=floor(t_sn(4));
% % %     n_sigma=sum(data(t_sn(1):t_sn(2)).^2)./(t_sn(2)-t_sn(1)); %Mean Variance of noise
% % %     s_sigma=sum(data(t_sn(3):t_sn(4)).^2)./(t_sn(4)-t_sn(3)); %Mean variance of signal
% % %     snr=10*log10(s_sigma/n_sigma); %Signal to Noise ratio (in decibels)?
% % % 
% % %     % Second threshold
% % % 
% % %     prob_false = 0.05;
% % %     r0 = 1;  % sucessive points
% % % 
% % %     switch r0
% % %         case 1
% % %             c=[1 -5 10 -10 5 -prob_false];
% % %         case 2
% % %             c=[-4 15 -20 10 0 -prob_false];
% % %         case 3
% % %             c=[6 -15 10 0 0 -prob_false];
% % %         case 4
% % %             c=[-4 5 0 0 0 -prob_false];
% % %         case 5
% % %             c=[1 0 0 0 0 -prob_false];
% % %         otherwise
% % %             c=[1 -5 10 -10 5 -prob_false];
% % %     end
% % %     rt=roots(c); %solving the polynomial equation
% % %     rtr=[];
% % %     for i=1:max(size(rt)) %
% % %         if isreal(rt(i))
% % %             rtr=[rtr rt(i)];
% % %         end
% % %     end
% % %     I=find(rtr>0 & rtr<1); %finding the maximum real root of the polynomial equation
% % %     p_zeta=rtr(I);         % probability for the threshold
% % %     pdk=exp(log(p_zeta)./(1+10.^(snr/10))); % obsolete -  probability that kth sample is above the threshold
% % %     pd=1-binocdf(r0-1,5,pdk);               % obsolete -  probability that signal samples, albeit corrupted by noise, are correctly recognized
% % %     
% % % 
% % % 
% % % 
% % %     assignin('base', 'this_better_work', x)
% % % 
% % %     % Detection
% % %     if mod(length(x),2)
% % %         x = x(1:end-1);
% % %     end
% % % 
% % % 
% % %     aux = x.^2; % auxiliary sequence
% % % 
% % %     if ~isrow(aux)
% % %        aux = aux' ;
% % %     end
% % % 
% % %     aux_seq = aux(1,1:2:end) + aux(1,2:2:end); % add pairs together
% % % 
% % % 
% % %     zeta=chi2inv(1-p_zeta,2)*n_sigma; % the first threshold value 
% % %     %zeta = -log(p_zeta)*2*n_sigma;   % this is equivalent to the above
% % %     disp(['Zeta: ', num2str(zeta)])
% % %     
% % %     assignin('base', 'zeta', zeta);
% % %     assignin('base', 'aux_seq', aux_seq);
% % %     assignin('base', 'data', data);
% % % 
% % % 
% % % 
% % %     dtct=zeros(1,max(size(aux_seq)));   % initialize vector
% % % 
% % %     % detection scanning (window of 5)
% % %     for i=1:5:max(size(aux_seq))-5
% % %         r=0;
% % %         for j=i:i+4 %detecting window's width is 5 points
% % %             if aux_seq(j)>zeta
% % %                 r=r+1;
% % %             end
% % %         end
% % %         if r>=r0;
% % %             dtct(i:i+4)=1;
% % %         end
% % %     end
% % % 
% % % 
% % %     % Find actual change points
% % %     % Convert dtct (series of 0s and 1s - indicating change) to actual pairs of
% % %     % on/off
% % %     test = find(dtct~=0);  
% % % end
% % %         
% % %         
% % %         
% % % 
% % %     index = 1;
% % %     count = 1;
% % %     myVec = zeros(1,length(test));
% % %     for i = 1:numel(test)-1
% % % 
% % %         if diff([test(i),test(i+1)]) ~= 1
% % %            myVec(1,count) = test(index);
% % %            myVec(1,count+1) = test(i);
% % % 
% % %            index = i+1;
% % %            count = count+2;
% % % 
% % %         end
% % %     end
% % %     myVec(count) = test(index);
% % %     myVec(count+1) = test(end);
% % %     myVec = myVec(myVec~=0);
% % %     % Remove all isolated spurious points (less than 15*2 ms)
% % %     mVec = myVec;
% % %     for i = 3:2:numel(mVec)-2
% % %         if diff([mVec(i),mVec(i+1)]) < 15 % Find differences less than 15*2
% % %             if diff([mVec(i-1),mVec(i)]) > 20 && diff([mVec(i+1),mVec(i+2)]) > 20 % If these points are greater than 20*2 away from another segment, remove
% % %                 mVec(i) = 0;
% % %                 mVec(i+1)=0;     
% % %             end
% % %         end        
% % %     end
% % %     mVec = mVec(mVec~=0);
% % % 
% % %     % Bridge points which are within 30 (2-4-6...)
% % %     % This is "bridging" points/segments - set to bridge greater than 15*2
% % % 
% % % 
% % %     printVec = mVec;
% % %     range = 1;
% % % 
% % % 
% % %     for i = 2:2:numel(printVec)-2
% % %         if diff([printVec(i),printVec(i+1)]) < round(range*1)
% % %             printVec(i) = 0;
% % %             printVec(i+1)=0;     
% % %         end        
% % %     end
% % %     printVec = printVec(printVec~=0);
% % % 
% % % 
% % % 
% % %     % Removing points (1-3-5-7...)
% % %     % This is "removal" of segments - set exactly as 15*2
% % %     for i = 1:2:numel(printVec)-1
% % %         if diff([printVec(i),printVec(i+1)]) < 15
% % %             printVec(i) = 0;
% % %             printVec(i+1)=0;     
% % %         end        
% % %     end
% % %     
% % %     
% % %     printVec = printVec(printVec~=0);
% % %     
% % %     printVec = printVec.*2;
% % %     axes(handles.axes1);
% % %     plot(data);
% % %     vline(printVec,'r')
% % %     pause(1)
% % %     cla;
% % %     close(inp);
% % % 
% % % end
% % % 
% % % 
% % % 
% % % 
% % % guidata(hObject,handles)
% % % 
% % % % uicontrol(handles.AGLRpb);
% % % % AGLRpb_Callback(handles.AGLRpb,[],handles);
% % % 
% % % 
% % % % --- Executes on button press in AGLRpb.
% % % function AGLRpb_Callback(hObject, eventdata, handles)
% % % % hObject    handle to AGLRpb (see GCBO)
% % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % handles    structure with handles and user data (see GUIDATA)
