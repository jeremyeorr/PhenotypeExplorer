%% Load data and preprocess
clear all
close all

directory = 'J:\PEOPLE\FACULTY\SANDS\O2PSG\_Scored\_MatlabAnalysis\';


addpath('E:\Work\MatlabFunctionDatabase');
addpath('E:\Work\MatlabFunctionDatabase\old matlab functions');

filename = '1309N2'; %1657/941 %checkedzeroflowactives: 1429, 941, 1264(fixed), 1722(fixed),1657,1343,1723,1469,1429,1309,815,533, 1710,1708
filenameanddir=[directory filename];

filehandle = matfile(filenameanddir);
w = whos('-file',filenameanddir);

load(filenameanddir,'PesToFlow','CCWCPAP','Pcrit','Veupnea','EdiToFlow','OSA','ArThresPes','ArThresEpi','ArthresEdi','Vactive','Vpassive','Varousal','LGediVarousal','SpO2wake') %Add more to this list

% Get channels / Rename channels
clear channelnameoptions
if 1 %load these
    channelnameoptions.Evts={'Evts','New_Evts'};
    channelnameoptions.Evts2={'Evts2','New_Evts2'};
    channelnameoptions.Epochs={'Epochs'};
    channelnameoptions.Flow={'Pnasal','Pmask','PNasal'}; %'Vflow'
    %channelnameoptions.Pmask={'Pmask','PMask'};
    %channelnameoptions.Pnasal={'Pnasal','Pmask','PNasal'}; %'Vflow'
    channelnameoptions.Position={'Position','Pos','pos','position'};
    channelnameoptions.SaO2={'SaO2','SpO2','Sat','Sao2','Spo2','O2sat','o2sat'};
    %channelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal','CO2'};
    %channelnameoptions.PO2={'PO2','pO2','O2_Ana','O2_anal'};
    %channelnameoptions.RC={'Thorax','RC','Chest','CHEST'};
    %channelnameoptions.ABD={'Abdomen','ABD','Abdom','ABDM','ABDO'};
    %channelnameoptions.ThRIP={'ThNoxRIP'};
    %channelnameoptions.AbRIP={'AbNoxRIP'};
    channelnameoptions.Pepi={'Pepi','Epi','PEpi','EPI','P6'};
    %channelnameoptions.Pes={'Pes'};
    %channelnameoptions.Edi={'Edi'};
    %channelnameoptions.GGpmax={'GGpmax','GGPmax','EMGggpmax'};
    channelnameoptions.EKG={'EKG','ECG'};
    channelnameoptions.EEG_C3_A2={'EEG_C3_A2','C3_A2'};
    channelnameoptions.EEG_C4_A1={'EEG_C4_A1','C4_A1'};
    channelnameoptions.EEG_O2_A1={'EEG_O2_A1','O2_A1'};
    channelnameoptions.EEG_F3_A2={'EEG_F3_A2','F3_A2'};
end

channelnamestemp=fieldnames(channelnameoptions);
for i=1:length(channelnamestemp)
    temp=eval(['channelnameoptions.' char(channelnamestemp(i))])
    foundamatch=0;
    for n=1:length(temp)
        %Does it exist?
        for j=1:length(w)
            if strcmp(w(j).name,char(temp(n)))
                eval([channelnamestemp{i} '=filehandle.' char(temp(n))]);
                foundamatch=1;
                break
            end
        end
        if foundamatch
            break
        end
    end
end

% Options

plotfigs=1;
saveontherun=1;
saveattheend=1;
usePnasalasflow=0;
analyseSaO2=0;
analyseO2=0;

%
if Flow.interval<0.008
    dt=0.008;
    %Flow.values=resample(Flow.values,round(dt/Flow.interval),1);
    Flow.values=downsample(Flow.values,round(dt/Flow.interval));
    Flow.interval=0.008;
    Flow.length=length(Flow.values);
end

%Setup Time array
dt=Flow.interval;
Fs=1/dt;

a=Flow.start;
Time=(a:dt:(a+dt*(Flow.length-1)))';

%% Epochs / Staging channel
EpochsXHz=interp1((Epochs.times+15),double(Epochs.codes(:,1)'),Time,'nearest','extrap');
%W=0,1=1,2=2,3=3,R=5,?=8

% Events info
%Evts list is just called Evts
if exist('New_Evts','var'), Evts=New_Evts; end
if exist('New_Evts2','var'), Evts2=New_Evts2; end

backupEvts = Evts;

%add extra row of zeros with time = inf just in case there is an event without an end in the first of two concatenated sets of events.
Evts.codes(size(Evts.codes,1)+1,:)=[0 0 0 0];
Evts.times(size(Evts.times,1)+1,:)=Inf;
Evts.text=[Evts.text;['                                            ']];

% new code to copy into other versions of code:
while 1
    if Evts.codes(1,1)==0
        Evts.codes(1,:)=[];
        Evts.times(1)=[];
        Evts.length=Evts.length-1;
        Evts.text(1,:)=[];
    else
        break
    end
end


if exist('Evts2','var')
while 1
    if Evts2.codes(1,1)==0
        Evts2.codes(1,:)=[];
        Evts2.times(1)=[];
        Evts2.length=Evts.length-1;
        Evts2.text(1,:)=[];
    else
        break
    end
end
end
if exist('Evts2','var')&&exist('Evts','var')
    %combine events lists
    Evts_backup1=Evts;
    Evts.times=[Evts.times;Evts2.times];
    Evts.codes=[Evts.codes;Evts2.codes];
    Evts.text=[Evts.text;Evts2.text];
    Evts.length=Evts.length+Evts2.length;
end


%delete second row of zeros if there are two concurrent rows of zeros
for i=length(Evts.codes):-1:2
    if Evts.codes(i)==0&&Evts.codes(i-1)==0
        i
        Evts.codes(i,:)=[];
        Evts.times(i)=[];
        Evts.text(i,:)=[];
        Evts.length=Evts.length-1;
        disp('delete_row_of_zeros');
    end
end

%add a row of zeros if there is a missing row of zeros. Assume end of event is the start of the next event (or end of file/Inf). 
for i=length(Evts.codes):-1:1
    if i==length(Evts.codes)
        if Evts.codes(i)>0
            Evts.codes=[Evts.codes;[0 0 0 0]];
            Evts.times=[Evts.times;Inf];
            Evts.text=[Evts.text;['                                            ']];
            Evts.length=Evts.length+1;
            disp('added_row_of_zeros to the end');
        end
    end
    if i>1
        if Evts.codes(i)>0&&Evts.codes(i-1)>0
            i
            Evts.codes=[Evts.codes(1:(i-1),:);[0 0 0 0];Evts.codes(i:end,:)];
            Evts.times=[Evts.times(1:i);Evts.times(i:end)];
            Evts.text=[Evts.text(1:(i-1),:);['                                            '];Evts.text(i:end,:)];
            Evts.length=Evts.length+1;
            disp('added_row_of_zeros');
        end
    end
end

% Make list of event types
clear EventTypeList_text EventTypeList_code tempi
EventTypeList_text=[];
for i=1:2:size(Evts.codes,1)
    codematch=0;
    for j=1:size(EventTypeList_text,1)
        if strcmp(char(EventTypeList_text(j,:)),Evts.text(i,:))
            'found code match'
            codematch=1;
            break
        end
    end
    if codematch==0
        tempi=size(EventTypeList_text,1);
        EventTypeList_code(tempi+1)=Evts.codes(i);
        EventTypeList_text(tempi+1,:)=Evts.text(i,:);
    end
end
EventTypeList_code=EventTypeList_code'
char(EventTypeList_text)

Evts_backup=Evts;

for i=1:size(Evts.codes,1)
    if strncmp(Evts.text(i,:),'AR',2)
        Evts.codes(i,1)=1;
    end
end


%remove zero rows:
Evts.codes(2:2:end,:)=[];

Evts.endtimes=Evts.times(2:2:end);
Evts.starttimes=Evts.times(1:2:end);
Evts.durations=Evts.endtimes-Evts.starttimes;
% Evts.timetoprevarousal=Evts.starttimes-[-Inf;Evts.endtimes(1:end-1)];
% Evts.timetonextarousal=[Evts.starttimes(2:end);Inf]-Evts.starttimes;

% Find sleep state and CPAP level for each event (arousal) %%faster code 2015-04-29
for i=1:length(Evts.starttimes)
    I=round((Evts.starttimes(i)-Time(1))/dt+1);
    Evts.epochs(i,:)=EpochsXHz(I);
%     try
%         Evts.CPAP(i,:)=Pmask.values(I);
%     catch me %for rerunning this code...
%         Evts.CPAP(i,:)=Pmask(I);
%     end
end

% Make arousals events in continuous time
EventsArXHz=zeros(length(Time),1); %Arousals.
EventsRespXHz=zeros(length(Time),1); %Arousals.
for x=1:length(Evts.codes)
    if Evts.codes(x)==1
        lefti=round(Evts.starttimes(x)-a)*Fs+1;
        righti=lefti+round(Evts.durations(x))*Fs;
        if righti==Inf
            righti=length(EventsArXHz);
        end
        %EventsArXHz(Time>=Evts.starttimes(x)&Time<Evts.endtimes(x))=Evts.codes(x); %%slower code because of search
        EventsArXHz(lefti:righti)=Evts.codes(x);
    end
    if Evts.codes(x)>1&&Evts.codes(x)<7
        lefti=round(Evts.starttimes(x)-a)*Fs+1;
        righti=lefti+round(Evts.durations(x))*Fs;
        if righti==Inf
            righti=length(EventsArXHz);
        end
        %EventsArXHz(Time>=Evts.starttimes(x)&Time<Evts.endtimes(x))=Evts.codes(x); %%slower code because of search
        EventsRespXHz(lefti:righti)=Evts.codes(x);
    end
end

clear lefti righti

clear EventTypeList_text EventTypeList_code Evts_backup Evts_backup1 Evts2 I backupEvts codematch foundamatch i j n temp tempi x

% Filter gently for use in analysis -- currently not analyzed
if 0
    filter_HFcutoff_butter0 = 4;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    PCO2.values = filtfilt(B_butter0,A_butter0,PCO2.values);
    
    if analyseO2
        pO2.values = filtfilt(B_butter0,A_butter0,pO2.values);
    else
        pO2.values = [];
    end
    
    backupflow=Flow.values;
    filter_HFcutoff_butter0 = 10;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    Flow.values = filtfilt(B_butter0,A_butter0,Flow.values);
    %Flow = filtfilt(B_butter0,A_butter0,Flow);
end

Evts.state_startevent = EpochsXHz(round((Evts.starttimes-a)/dt)+1);

%% Median event duration all example
median(Evts.durations(Evts.codes>1&Evts.codes<7))

%% Make all channels the same length as Flow, resample if needed

clear channels_all
for i=4:length(channelnamestemp)
    channels_all{i-3}=channelnamestemp{i};
end

clear channels
for i=1:length(channels_all)
    if exist(channels_all{i},'var')
        channels{i}=channels_all{i};
    end
end

clear channels_all
%Add code here which resamples if channels are not the same sampling rate
for i=1:length(channels)
    if ~isempty(channels{i})
        dtnew=eval([channels{i} '.interval']);
        if dtnew~=dt
            disp(['resampling: ' channels{i} '.interval = ' num2str(eval([channels{i} '.interval']))]);
            if mod(dt/dtnew,1)~=0 %unfinished coding
                eval([channels{i} '.values = resample(' channels{i} '.values,round(1/dt),round(1/dtnew));' ]);
            else
                %eval([channels{i} '.values = resample(' channels{i} '.values,round(dt/' channels{i} '.interval),1);' ]);
                eval([channels{i} '.values = downsample(' channels{i} '.values,round(dt/dtnew));' ]);
            end
            eval([channels{i} '.interval = dt;']);
        end
    end
end

%Make channels the same length (sometimes these are off by up to 20 samples)
for i=1:length(channels)
    if ~isempty(channels{i})
        while length(eval([channels{i} '.values']))~=length(Time)
            if length(eval([channels{i} '.values']))<length(Time)
                eval([channels{i} '.values(end+1)=' channels{i} '.values(end);']); %add a sample to the end if the channel is too short by 1
            elseif length(eval([channels{i} '.values']))>length(Time)
                eval([channels{i} '.values(end)=[];']); %delete a sample from the end if the channel is too long by 1
            end
        end
    end
end

%Rename channel data: for example, Flow = Flow.values; clear Flow structure variable
for i=1:length(channels)
    if ~isempty(channels{i})
        eval(['temp=' channels{i} '.values;']);
        eval([channels{i} '=temp;']);
        clear temp
    end
end

%channels{length(channels)+1}='Time';
channels{length(channels)+1}='EventsArXHz';
channels{length(channels)+1}='EpochsXHz';

%% Plot preliminary figure to check
figure(1);
X=5;

if exist('Pepi','var')
    X=X+1;
end

dsf=5;

ax1(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(EpochsXHz,dsf));
ax1(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(EventsArXHz,dsf));
ax1(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(EventsRespXHz,dsf));
ax1(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(SaO2,dsf));
ax1(5)=subplot(X,1,5); plot(downsample(Time,dsf),downsample(Flow,dsf));
i=X-1;

if exist('Pepi','var')
    i=i+1;
    ax1(i)=subplot(X,1,i); plot(downsample(Time,dsf),downsample(Pepi,dsf));
end

%ax1(3)=subplot(X,1,2); plot(Time,Vol2,Time(index),Vol2(index),'.',Time(index2),Vol2(index2),'.');
linkaxes(ax1,'x');


%% Beta Power analysis (~10 s)

EEGsignals = {'EEG_O2_A1','EEG_F3_A2','EEG_C4_A1','EEG_C3_A2'};
for i=length(EEGsignals):-1:1
    if ~exist(EEGsignals{i})
        EEGsignals(i)=[];
    end
end

if exist('EKG','var')
    % Load or create new "clean" EEG signals
    lefttemplatetime = 0.05; righttemplatetime = 0.05;
    contaminationmagnitudethreshold = 4;
    leftdi=round(1/dt*lefttemplatetime); rightdi=round(1/dt*righttemplatetime);
    clear temp contaminationmagnitude contaminationmagnitudepost
    processEEG=1;
    loadprocessedEEGifexists = 1; saveprocessedEEG = 1; foundamatch=0;
    if loadprocessedEEGifexists
        for n=1:length(EEGsignals)
            foundamatch=0;
            temp{n}=[EEGsignals{n} '_clean'];
            for j=1:length(w)
                foundamatch=strcmp(w(j).name,temp{n});
                if foundamatch
                    eval([temp{n} '=filehandle.' temp{n} ';']);
                    disp(['found:' temp{n}]);
                    break
                end
            end
        end
        if foundamatch
            ECG_peak_i = EKGpeakdetection(EKG,Time,dt);
            for i=1:length(EEGsignals)
                [~,template] = crosscontaminated(eval([EEGsignals{i}]),ECG_peak_i,leftdi,rightdi,1);
                contaminationmagnitude(i)=max(abs(template));
            end
            for i=1:length(EEGsignals)
                [~,template] = crosscontaminated(eval([EEGsignals{i} '_clean.values']),ECG_peak_i,leftdi,rightdi,1);
                contaminationmagnitudepost(i)=max(abs(template));
            end
        end
    end
    
    if (~foundamatch||loadprocessedEEGifexists==0)&&processEEG
        for i=1:length(EEGsignals)
            eval([EEGsignals{i} '_clean=' EEGsignals{i} ';']);
        end
        % EKG artifact removal
        ECG_peak_i = EKGpeakdetection(EKG,Time,dt);
        
        for i=1:length(EEGsignals)
            [~,template] = crosscontaminated(eval([EEGsignals{i}]),ECG_peak_i,leftdi,rightdi,1);
            contaminationmagnitude(i)=max(abs(template));
        end
        disp(['Decontaminating ' num2str(sum(contaminationmagnitude>contaminationmagnitudethreshold)) ' EEG signals']);
        for i=1:length(EEGsignals)
            if contaminationmagnitude(i)>contaminationmagnitudethreshold
                eval([EEGsignals{i} '_clean=crosscontaminated(' EEGsignals{i} ',ECG_peak_i,leftdi,rightdi,0);']);
                %e.g. EEG_F3_A2 = crosscontaminated(EEG_F3_A2_old,ECG_peak_i,leftdi,rightdi,0);
            end
        end
        %end
        % after
        for i=1:length(EEGsignals)
            [~,template] = crosscontaminated(eval([EEGsignals{i} '_clean']),ECG_peak_i,leftdi,rightdi,1);
            contaminationmagnitudepost(i)=max(abs(template));
        end
        
        if saveprocessedEEG
            for i=1:length(EEGsignals)
                varlist{i} = [EEGsignals{i} '_clean'];
                data{i}=eval(varlist{i});
            end
            ConvertToSpikeStructandSave(varlist,data,dt,directory,filename);
            clear data
        end
    else
        for i=1:length(EEGsignals)
            eval([EEGsignals{i} '_clean=' EEGsignals{i} ';']);
        end
    end
else
    for i=1:length(EEGsignals)
        eval([EEGsignals{i} '_clean=' EEGsignals{i} ';']);
    end
    contaminationmagnitude=NaN;
    contaminationmagnitudepost=NaN;
end

% BetaPower
clear WakeSleepInfo
WakeSleepInfo.contaminationmagnitude=contaminationmagnitude;
WakeSleepInfo.contaminationmagnitudepost=contaminationmagnitudepost;
for i=1:length(EEGsignals)
    WakeSleepInfo.EEGoptions{i} = [EEGsignals{i}];
end

for i=1:length(WakeSleepInfo.EEGoptions)
    if ~exist(WakeSleepInfo.EEGoptions{i},'var')
        eval([WakeSleepInfo.EEGoptions{i} '=[];']);
    end
end
%WakeSleepInfo.alpha = [8 12];
WakeSleepInfo.beta = [16 32]; %Best upper: 1723:62, 1343:32(vs 24,48), 533:62(not substantially better than 32); 815:62; 1429:32(vs 24,48); 1469:48(vs 32,62); 1710:32 (vs 62); 941:32(vs 24,62; 48 better slightly)
WakeSleepInfo.fft_length = 128; %256 %Best: 1710,1723:256, 1343:512
WakeSleepInfo.Fpower=0; %scale PSD by frequency^Fpower; Optimal for 1429=1(others are ok); 941,1723=2, 1469:2(othersok), 1343=4, 533=7; 815=3; 1309:2(3 if beta(2)=32); 1708:5(othersok); 1710:2 (3 ok, 1 worse)
WakeSleepInfo.Foverlap = 0.75; %0.75 is better than 0.5 in 1723,1343,1309
WakeSleepInfo.medianfiltertime = 6; %10 is better than 15 and 8 in 1723,1343 not 1309
WakeSleepInfo.nearbyduration=300;
WakeSleepInfo.useonlySpO2on=0; %if zero, also uses nearness to sleep (nearbyduration 180 s) to estimate whether EEG is on
WakeSleepInfo.scoredarousalsinwake=1; %Lauren's special scoring that ignores some clinical rules regarding arousals.
plotfigs=0;
if WakeSleepInfo.fft_length*(1-WakeSleepInfo.Foverlap)<32 %maintain maximum of 0.25 s step size, saving time...
    temp = 1-32/WakeSleepInfo.fft_length;
    if temp<0, temp = 0; end
    WakeSleepInfo.Foverlap=temp;
end

[WakeSleep,WakeSleepInfo]=BetaPower(WakeSleepInfo,EpochsXHz,EventsArXHz,Time,SaO2,Flow,plotfigs);
[~,bestEEG]=max(WakeSleepInfo.AUC_M);
bestEEGstr = [WakeSleepInfo.EEGoptions{bestEEG} '_clean'];

%plot
plotWakeSleep=1;
if plotWakeSleep
    figure(30); set(gcf,'color',[1 1 1]);
    ax30(1)=subplot(4,1,1); plot(Time,EventsArXHz,Time,EpochsXHz,Time,WakeSleep); ylim([-0.1 5.1]);
    box('off');
    ax30(2)=subplot(4,1,2);
    plot(Time,EventsArXHz,Time,WakeSleep); box('off'); hold('on');
    %plot(Palphabeta_time,0.5*Palphabeta_transformed_medfilt+0.5,'k'); hold('off'); box('off'); ylim([-0.1 1.1]); ylabel('Prwake');
    %
    %WakeSleep = interp1(Palphabeta_time,Palphabeta_transformed_medfilt,Time,'linear');
    %
    ax30(3)=subplot(4,1,3); plot(Time,eval(bestEEGstr),'k');  ylim([-50 50]); box('off');
    
    ax30(4)=subplot(4,1,4);
    plot(Time,Flow,'k'); ylim([-1 1]); box('off'); hold('on');
    
    tempT=Time; tempT(WakeSleep<-0.5|WakeSleep>=0)=NaN;
    tempV=Flow; tempV(WakeSleep<-0.5|WakeSleep>=0)=NaN;
    plot(tempT,tempV,'color',[0.33 0 0]); ylim([-1 1]); box('off'); hold('on');
    clear tempT tempV
    
    tempT=Time; tempT(WakeSleep<0|WakeSleep>=0.5)=NaN;
    tempV=Flow; tempV(WakeSleep<0|WakeSleep>=0.5)=NaN;
    plot(tempT,tempV,'color',[0.67 0 0]); ylim([-1 1]); box('off'); hold('on');
    clear tempT tempV
    
    tempT=Time; tempT(WakeSleep<0.5)=NaN;
    tempV=Flow; tempV(WakeSleep<0.5)=NaN;
    plot(tempT,tempV,'color',[1 0 0]); ylim([-1 1]); box('off'); hold('off');
    clear tempT tempV
    
    linkaxes(ax30,'x');
end

channels{length(channels)+1}='WakeSleep';


%% Find Arousal Threshold (Pes)

%Use to cursor to exclude (right click)
%or select Pepi breath to analyze x=time of nadir, y=baseline (start
%inspiration).

global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
signallist={'[EpochsXHz EventsArXHz]','EventsRespXHz','Pmask','Flow','Pepi'};
PlotAndSelectData(signallist,Time,clearvalues);
It = Time(find(diff(EventsArXHz)==1));

%% Analyze ARthres Pepi
% no manual analysis needed
ArThresPes.x_select = xvalues(:,2);
ArThresPes.y_select = yvalues(:,1);
%%
ArThresPes.data=[];
Pesnadirabs=[];
ArThresPes.CPAP=[]; %note this is the Pmask measured at the Pes nadir rather than during zero flow.
ArThresPes.data_t=[];
ArThresPes.data_i=[];

for i=1:length(ArThresPes.x_select)
    [Pesnadir_uncorrected,index]=min(Pepi(Time>ArThresPes.x_select(i)-1&Time<ArThresPes.x_select(i)+1));
    ArThresPes.data_t(i)=ArThresPes.x_select(i)-1+(index-1)*dt;
    ArThresPes.data_i(i)=find(Time>ArThresPes.data_t(i),1,'first');
    %plot
    lowerX=ArThresPes.x_select(i)-15;
    upperX=ArThresPes.x_select(i)+15;
    ArThresPes.data(i)=Pesnadir_uncorrected-ArThresPes.y_select(i);
    %ArThresPes.CPAP(i)=Pmask(find(Time>ArThresPes.data_t(i),1));
    %Pesnadirabs(i)=ArThresPes.data(i)+Pmask(ArThresPes.data_i(i));
    %plot
    %ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pepi,dsf),downsample(Time,dsf),downsample(Time,dsf)*0+estimatedzeroPepi,'k:',ArThresEpi.data_t(i),Pepi(ArThresEpi.data_i(i)),'r.');
    ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pepi,dsf),ArThresPes.data_t(i),Pepi(ArThresPes.data_i(i)),'r.');
    xlim([lowerX upperX]);
    set(ax2(3),'ylim',[min(Pepi(Time>lowerX&Time<upperX)) prctile(Pepi(Time>lowerX&Time<upperX),95)]);
    hold('off');
end

%% Summary analysis
ArThresPes.state=NaN*ArThresPes.data_i;
for i=1:length(ArThresPes.data_i)
    try
        ArThresPes.state(i)=max([EpochsXHz(ArThresPes.data_i(i)-15/dt) EpochsXHz(ArThresPes.data_i(i)-0/dt)]);
    catch me
    end
end
criteria = (ArThresPes.state<5);
ArThresPes.mean = nanmean(ArThresPes.data(criteria));
ArThresPes.median = nanmedian(ArThresPes.data(criteria));
ArThresPes.N = sum(~isnan(ArThresPes.data(criteria)));
ArThresPes.SEM = nanstd(ArThresPes.data(criteria))/ArThresPes.N^0.5;

ArThresPes.mean_abs = nanmean(Pesnadirabs(criteria));
ArThresPes.median_abs = nanmedian(Pesnadirabs(criteria));
ArThresPes.SEM_abs = nanstd(Pesnadirabs(criteria))/ArThresPes.N^0.5;

ArThresPes.cov = nanstd(ArThresPes.data(criteria))/-ArThresPes.mean;

ArThresPes.onCPAP = nanmean(ArThresPes.data(criteria&(ArThresPes.CPAP<-1|ArThresPes.CPAP>1)));
ArThresPes.offCPAP = nanmean(ArThresPes.data(criteria&ArThresPes.CPAP>=-1&ArThresPes.CPAP<=1));
ArThresPes.offCPAPN = sum(criteria&~isnan(ArThresPes.data)&ArThresPes.CPAP>=-1&ArThresPes.CPAP<=1);

%% --save
save(filenameanddir,'ArThresPes','-append');