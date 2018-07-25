%% Load data and preprocess
clear all
close all
%F:\Work\Projects SS\Overweight and obese nonapneics\Draft and Data\1676N1BariatricTrait
%filenameanddir='J:\PEOPLE\POST DOC\SANDS\Phenotyping and Edi\1722';
if 0
    directory = 'E:\Work\Projects SS\Phenotyping and Edi\';
else
    directory='J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\SpikeFiles\';
end
addpath('E:\Work\MatlabFunctionDatabase');
addpath('E:\Work\MatlabFunctionDatabase\old matlab functions');
%directory = 'J:\PEOPLE\FACULTY\SANDS\Nox\';
filename = '1771'; %1657/941 %checkedzeroflowactives: 1429, 941, 1264(fixed), 1722(fixed),1657,1343,1723,1469,1429,1309,815,533, 1710,1708
filenameanddir=[directory filename];

%note 1727 Pes swings are large during wake, larger vs Pmask by 1/0.67-1/0.77.
%PmusVsEdi: 1657 (¬0.55), 1429 and 941 (<=0.5).

%1469 - could not get Ccw
filehandle = matfile(filenameanddir);
w = whos('-file',filenameanddir);

load(filenameanddir,'LG','PesToFlow','CCWCPAP','Pcrit','Veupnea','EdiToFlow','OSA','ArThresPes','ArThresEpi','ArthresEdi','Vactive','Vpassive','Varousal','LGediVarousal','SpO2wake') %Add more to this list
% load(filenameanddir,'OSA') %Add more to this list

% Get channels / Rename channels
clear channelnameoptions
if 1 %load these
    channelnameoptions.Evts={'Evts','New_Evts'};
    channelnameoptions.Evts2={'Evts2','New_Evts2'};
    channelnameoptions.Epochs={'Epochs'};
    channelnameoptions.Flow={'Flow','Vflow','VFlow'}; %'Vflow'
    channelnameoptions.Pmask={'Pmask','PMask'};
    channelnameoptions.Position={'Position','Pos','pos','position'};
    channelnameoptions.SaO2={'SaO2','SpO2','Sat','Sao2','Spo2','O2sat','o2sat'};
    channelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal','CO2'};
    channelnameoptions.PO2={'PO2','pO2','O2_Ana','O2_anal'};
    channelnameoptions.RC={'Thorax','RC','Chest','CHEST'};
    channelnameoptions.ABD={'Abdomen','ABD','Abdom','ABDM','ABDO'};
    channelnameoptions.ThRIP={'ThNoxRIP'};
    channelnameoptions.AbRIP={'AbNoxRIP'};
    channelnameoptions.Pepi={'Pepi','Epi','PEpi','EPI','P6'};
    channelnameoptions.Pes={'Pes_clean','Pes'};
    channelnameoptions.Edi={'Edi'};
    channelnameoptions.GGpmax={'GGpmax','GGPmax','EMGggpmax'};
    channelnameoptions.EKG={'EKG','ECG'};
    channelnameoptions.EEG_C3_A2={'EEG_C3_A2','C3_A2'};
    channelnameoptions.EEG_C4_A1={'EEG_C4_A1','C4_A1'};
    channelnameoptions.EEG_O2_A1={'EEG_O2_A1','O2_A1'};
    channelnameoptions.EEG_F3_A2={'EEG_F3_A2','F3_A2'};
    channelnameoptions.WakeSleep={'WakeSleep'};
    channelnameoptions.FlowEdi={'FlowEdi'};
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

% if ThNoxRIP.interval~=0.04
%     ThNoxRIP.interval=0.04;
%     AbNoxRIP.interval=0.04;
% end

%Setup Time array
dt=Flow.interval;
Fs=1/dt;

a=Flow.start;
Time=(a:dt:(a+dt*(Flow.length-1)))';

%% Epochs / Staging channel
EpochsXHz=interp1((Epochs.times+15),double(Epochs.codes(:,1)'),Time,'nearest','extrap');

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
    Evts.times=[Evts.times;Evts2.times];
    Evts.codes=[Evts.codes;Evts2.codes];
    Evts.text=[Evts.text;Evts2.text];
    Evts.length=Evts.length+Evts2.length;
end


%delete second row of zeros if there are two concurrent rows of zeros
for i=length(Evts.codes):-1:2
    if Evts.codes(i)==0&&Evts.codes(i-1)==0
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
    try
        Evts.CPAP(i,:)=Pmask.values(I);
    catch me %for rerunning this code...
        Evts.CPAP(i,:)=Pmask(I);
    end
end

% Make arousals events in continuous time
EventsArXHz=zeros(length(Time),1); %Arousals.
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

% Plot preliminary figure to check
figure(1);
X=4;

if exist('PCO2','var')
    X=X+1;
end

if exist('Pepi','var')
    X=X+1;
elseif exist('Pes','var')
    X=X+1;
end
if exist('GGpmax','var')
    X=X+1;
end
if exist('Edi','var')
    X=X+1;
end
dsf=20;
ax1(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(Flow,dsf));
ax1(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Pmask,dsf));
ax1(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(EpochsXHz,dsf));
ax1(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(EventsArXHz,dsf));
i=4;
if exist('PCO2','var')
    i=i+1;
    ax1(i)=subplot(X,1,i); plot(downsample(Time,dsf),downsample(PCO2,dsf));
end
if exist('Pes','var')
    i=i+1;
    ax1(i)=subplot(X,1,i); plot(downsample(Time,dsf),downsample(Pes,dsf));
elseif exist('Pepi','var')
    i=i+1;
    ax1(i)=subplot(X,1,i); plot(downsample(Time,dsf),downsample(Pepi,dsf));
end
if exist('GGpmax','var')
    i=i+1;
    ax1(i)=subplot(X,1,i); plot(downsample(Time,dsf),downsample(GGpmax,dsf));
end
if exist('Edi','var')
    i=i+1;
    ax1(i)=subplot(X,1,i); plot(downsample(Time,dsf),downsample(Edi,dsf));
end
%ax1(3)=subplot(X,1,2); plot(Time,Vol2,Time(index),Vol2(index),'.',Time(index2),Vol2(index2),'.');
linkaxes(ax1,'x');

% Plot preliminary figure RIP
if 0
    figure(2);
    X=4;
    
    dsf=1;
    ax1(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(Flow,dsf));
    ax1(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Pmask,dsf));
    ax1(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(ThNoxRIP,dsf));
    ax1(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(AbNoxRIP,dsf));
    linkaxes(ax1,'x');
end

%% Position channel fix using text
% Fix Pos data if necessary
textfilename=[directory filename '_pos.txt'];
if exist(textfilename,'file')==2
    'found position file'
    %Position.values_original=Position.values;
    [col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
    for i=1:(length(col1)-1)
        Position(Time>=col1(i)&Time<col1(i+1))=col2(i);
    end
    Position(Time>=col1(end))=col2(end);
end

%% Remove artifact in all signals using text files
if 0 %causes trouble!
    for j=1:length(channels)
        textfilename=[directory filename '_' channels{j} '_art.txt'];
        if exist(textfilename,'file')==2
            ['found ' channels{j} ' artifact file and removed artifact']
            %Position.values_original=Position.values;
            [col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
            for i=1:length(col1)
                %dttemp=eval([channels{j} '.interval']);
                lefti=round((col1(i)-Time(1))/dt+1);
                if lefti<1, lefti=1; end
                righti=round((col2(i)-Time(1))/dt+1);
                if righti>length(eval(channels{j})), righti=length(eval(channels{j})); end
                eval([channels{j} '(lefti:righti)=NaN;']);
            end
            Percent_removed = sum(isnan(eval(channels{j})))/length(eval(channels{j}))*100
            figure(5); plot(Time,eval([channels{j}]));
        end
    end
end

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
%% Save channel WakeSleep
if ~isstruct(WakeSleep)
    temp = WakeSleep;
end
clear WakeSleep;
WakeSleep=WakeSleepInfo
WakeSleep.values=temp;
WakeSleep.interval=dt;
if saveontherun
    save(filenameanddir,'WakeSleep','-append');
end
clear WakeSleep
WakeSleep=temp;


%% Pes decontaminate
Pes_clean = PesRemoveEKG(Pes,EKG,dt);
if ~isstruct(Pes_clean)
    temp = Pes_clean;
end
clear Pes_clean;
Pes_clean.values=temp;
Pes_clean.interval=dt;
if saveontherun
    save(filenameanddir,'Pes_clean','-append');
end
clear Pes_clean
Pes_clean=temp;

%%

Pes=Pes_clean;
clear Pes_clean
%% Find Pcrit runs, using Xminute screens scrolling by Yminutes

global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Flow','Pes'};
PlotAndSelectData(signallist,Time,clearvalues);
It = Time(find(diff(EventsArXHz)==1));


%% Test scroll scrolling by Yminutes
% 
% global xvalues yvalues range ax2
% clearvalues=1;
% range=5*60;
% signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Flow','Pes','Edi'};
% PlotAndSelectData(signallist,Time,clearvalues);
% It = Time(find(diff(EventsArXHz)==1));
% 
% EpochsXHz
% EventsArXHz
% WakeSleep


%%
Pcrit.runtimes=xvalues;

if 0&saveontherun
    save(filenameanddir,'Pcrit','-append');
end

%% Analyze Pcrit

%left click ok | right click if zero flow | ...
%right click left of screen to exclude | left click left of screen to go back

i=1;
exclude=[];
Pcrit.VE=[];
Pcrit.Peakflow=[];
Pcrit.CPAP=[];
deltaX=3;
while i<=size(Pcrit.runtimes,1)
    try
        lowerX=Pcrit.runtimes(i,1)-deltaX;
        upperX=Pcrit.runtimes(i,2)+deltaX;
        I=find(Time>lowerX&Time<upperX);
        figure(33);
        
        ax33(1)=subplot(2,1,1); plot(Time(I),cumsum(Flow(I)-median(Flow(I)))*dt);
        ax33(2)=subplot(2,1,2); plot(Time(I),Flow(I)-median(Flow(I)));
        xlabel(num2str(i))
        exclude(i)=0; %default
        
        [x1,~,button(i)] = ginput(1);
        
        if x1<Time(I(1))&&button(i)==3 %exclude for other reason
            exclude(i)=1;
            i=i+1;
            continue
        elseif button(i)==3 %zero flow
            Pcrit.Peakflow_list{i}=0;
            Pcrit.Peakflow(i)=0;
            Pcrit.VE_list{i}=0;
            Pcrit.VE(i)=0;
            Pcrit.CPAP_list{i}=Pmask(find(Time>Pcrit.runtimes(i,1),1));
            Pcrit.CPAP(i)=Pcrit.CPAP_list{i};
            exclude(i)=0;
            i=i+1;
            continue
        elseif x1<Time(I(1))&&button(i)==1&&i>1 %go back
            i=i-2;
            continue
        end
        
        [~,~,BB_i_start,~,~,BB_t,VI,VE,Ttot,leak,IEratio,VT,Vpeak,Vpeakmean]=Vflowanalysis3(Time(I),Flow(I),0);
        xlabel(num2str(i))
        Pcrit.Peakflow_list{i}=Vpeak;
        Pcrit.Peakflow(i)=median(Vpeak);
        Pcrit.VE_list{i}=VI*60;
        Pcrit.VE(i)=median(VI)*60;
        Pcrit.CPAP_list{i}=Pmask(I(1)-1+BB_i_start);
        Pcrit.CPAP(i)=median(Pcrit.CPAP_list{i});
        Pcrit.State(i)=EpochsXHz(I(1)-1+BB_i_start(1));
        
        pause
    catch me
        exclude(i)=1;
        ['error:i=' num2str(i)]
    end
    i=i+1;
    
end

Pcrit.VE(button==3)=0;
Pcrit.Peakflow(button==3)=0;
Pcrit.VE(exclude==1|isnan(Pcrit.CPAP))=[];
Pcrit.Peakflow(exclude==1|isnan(Pcrit.CPAP))=[];
Pcrit.CPAP(exclude==1|isnan(Pcrit.CPAP))=[];


%% --save
save(filenameanddir,'Pcrit','-append');

%% Pcrit analysis and plots
plotwithnumbers=1;
fh=figure;
global fixedslope upperlimit;
fixedslope=NaN;
upperlimit.on=0;
Pcrit.minPmask=99;
%Pcrit.forcedslopes = [3.6519 0.9120];
Pcrit.forcedslopes = [NaN NaN];

if 1
    Pcrit.include_data=1:length(Pcrit.CPAP); %1:length(Pcrit.CPAP) [20:45];
    removelist1 = []; %[10 18 4 5 12 26]
    removelist1 = [removelist1 find(Pcrit.CPAP>Pcrit.minPmask)];
    Pcrit.include_data(removelist1)=[];
end

Pmaskdata=Pcrit.CPAP(Pcrit.include_data);
VdotMaxdata=Pcrit.Peakflow(Pcrit.include_data)*60;
VdotEdata=Pcrit.VE(Pcrit.include_data);

set(fh,'color',[1 1 1]);
ax1(1)=subplot(2,1,1);
plotpoints=1;

if isfield(Pcrit,'forcedslopes')&&~isnan(Pcrit.forcedslopes(1))
    fixedslope = Pcrit.forcedslopes(1);
end

[Pcrit.VdotcritSEM,Pcrit.Vdotcritvalue,Pcrit.PcritSEM,Pcrit.Pcritvalue,Pcrit.PVSlope_Vmax]=plotregressionwithSEM_linearsegment(Pmaskdata,VdotMaxdata,0,plotpoints)
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('Vpeak','fontname','arial narrow','fontsize',12);
ax1(2)=subplot(2,1,2);

if isfield(Pcrit,'forcedslopes')&&~isnan(Pcrit.forcedslopes(1))
    fixedslope = Pcrit.forcedslopes(2);
end

[Pcrit.VcritSEM,Pcrit.Vcritvalue,~,~,Pcrit.PVSlope_VE]=plotregressionwithSEM_linearsegment([Pmaskdata],[VdotEdata],0,plotpoints);
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('VE','fontname','arial narrow','fontsize',12);
linkaxes(ax1,'x');

if plotwithnumbers
    for i=1:length(Pcrit.include_data)
        ax2(1)=subplot(2,1,2);
        h=text(Pmaskdata(i)+0.1,VdotEdata(i)+0.1,num2str(Pcrit.include_data(i)),'fontname','arial narrow','fontsize',9);
    end
end
xlim([min(Pmaskdata) max(Pmaskdata)])
%% Pcrit analysis and plots - median binned data
Pcrit.usealldataformedianbinned=0;
binwidth=1;
CPAPbins1=(floor(min(Pcrit.CPAP))-0):binwidth:(ceil(max(Pcrit.CPAP))-binwidth*0.5);
Ndata=[]; VdotMaxdata=[]; VdotEdata=[]; Pmaskdata=[];
VdotEdataU=[];
VdotEdataL=[];
VdotMaxdataU=[];
VdotMaxdataL=[];

if 1
    Pcrit.include_data2=1:length(Pcrit.CPAP);
    %removelist2 = [];
    removelist2 = removelist1;
    removelist2 = [removelist2 find(Pcrit.CPAP>Pcrit.minPmask)];
    Pcrit.include_data2(removelist2)=[];
end

if Pcrit.usealldataformedianbinned
    Pmaskdata_=Pcrit.CPAP;
    VdotMaxdata_=Pcrit.Peakflow*60;
    VdotEdata_=Pcrit.VE;
else
    Pmaskdata_=Pcrit.CPAP(Pcrit.include_data2);
    VdotMaxdata_=Pcrit.Peakflow(Pcrit.include_data2)*60;
    VdotEdata_=Pcrit.VE(Pcrit.include_data2);
end

for i=1:length(CPAPbins1)
    tempI=find(Pmaskdata_>CPAPbins1(i)&Pmaskdata_<=(CPAPbins1(i)+binwidth));
    Ndata(i)=length(tempI);
    VdotMaxdata(i)=median(VdotMaxdata_(tempI));
    VdotMaxdataU(i)=prctile(VdotMaxdata_(tempI),75);
    VdotMaxdataL(i)=prctile(VdotMaxdata_(tempI),25);
    VdotEdata(i)=median(VdotEdata_(tempI));
    VdotEdataU(i)=prctile(VdotEdata_(tempI),75);
    VdotEdataL(i)=prctile(VdotEdata_(tempI),25);
    Pmaskdata(i)=mean(Pmaskdata_(tempI));
end

plotwithnumbers=1;
fh=figure;
fixedslope=NaN;

I2 = find(Pmaskdata>Pcrit.minPmask|isnan(Pmaskdata));
Ndata(I2)=[];
VdotMaxdata(I2)=[];
VdotEdata(I2)=[];
Pmaskdata(I2)=[];
VdotEdataU(I2)=[];
VdotEdataL(I2)=[];
VdotMaxdataU(I2)=[];
VdotMaxdataL(I2)=[];

set(fh,'color',[1 1 1]);
ax1(1)=subplot(2,1,1);
plotpoints=1;


if isfield(Pcrit,'forcedslopes')&&~isnan(Pcrit.forcedslopes(1))
    fixedslope = Pcrit.forcedslopes(1);
end

[Pcrit.VdotcritSEM_medianmethod,Pcrit.Vdotcritvalue_medianmethod,Pcrit.PcritSEM_medianmethod,Pcrit.Pcritvalue_medianmethod,Pcrit.PVSlope_Vmax_medianmethod]=plotregressionwithSEM_linearsegment(Pmaskdata,VdotMaxdata,0,plotpoints)
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('Vpeak','fontname','arial narrow','fontsize',12);
for i=1:length(Ndata)
    plot([Pmaskdata(i) Pmaskdata(i)],[VdotMaxdataL(i) VdotMaxdataU(i)],'k');
end

if isfield(Pcrit,'forcedslopes')&&~isnan(Pcrit.forcedslopes(1))
    fixedslope = Pcrit.forcedslopes(2);
end

ax1(2)=subplot(2,1,2);
[Pcrit.VcritSEM_medianmethod,Pcrit.Vcritvalue_medianmethod,~,~,Pcrit.PVSlope_VE_medianmethod]=plotregressionwithSEM_linearsegment([Pmaskdata],[VdotEdata],0,plotpoints);
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('VE','fontname','arial narrow','fontsize',12);
linkaxes(ax1,'x');
for i=1:length(Ndata)
    plot([Pmaskdata(i) Pmaskdata(i)],[VdotEdataL(i) VdotEdataU(i)],'k');
end

%% --save
save(filenameanddir,'Pcrit','-append');


%% Misc select respiratory data
global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Flow','Pes'};
PlotAndSelectData(signallist,Time,clearvalues);

%% Analyze misc respiratory data
i=1
I=find(Time>xvalues(i,1)&Time<xvalues(i,2));
Time1=Time(I);
Flow1=Flow(I);
[~,~,BB_i_start,BB_i_mid,BB_i_end,BB_t,VI,VE,Ttot,leak,IEratio,VT,Vpeak,Vpeakmean]=Vflowanalysis3(Time1,Flow1,0)


%% Find individual Pepi/GG segments, using Xminute screens scrolling by Yminutes

figure(2)
Xminute=10;
Yminute=5;
dsf=5;
X=4;
ax2(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(EpochsXHz,dsf),'b',downsample(Time,dsf),downsample(EventsArXHz,dsf),'r');
if exist('WakeSleep'), hold('on'); plot(downsample(Time,dsf),downsample(WakeSleep,dsf),'color',[1 0.3 0.3]); hold('off'); end
ax2(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Pmask,dsf));
ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pepi,dsf));
ax2(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(GGpmax,dsf));
linkaxes(ax2,'x')
upperX=Xminute*60
lowerX=0;
xlim([lowerX upperX])
ylim([prctile(GGpmax,1) prctile(GGpmax,95)])
rangedata=[];
while 1
    clear x1 x2 y1 y2 button
    [x1,y1,button]=ginput(1);
    if button==3
        break
    end
    if x1>upperX
        upperX=upperX+Yminute*60;
        lowerX=lowerX+Yminute*60;
    elseif x1<lowerX
        upperX=upperX-Yminute*60;
        lowerX=lowerX-Yminute*60;
    else
        [x2,y2,button]=ginput(1);
        rangedata=[rangedata;x1,x2]
    end
    xlim([lowerX upperX])
    set(ax2(3),'ylim',[prctile(Pepi(Time>lowerX&Time<upperX),5)-5 prctile(Pepi(Time>lowerX&Time<upperX),95)+5]);
    set(ax2(4),'ylim',[prctile(GGpmax(Time>lowerX&Time<upperX),1)-1 prctile(GGpmax(Time>lowerX&Time<upperX),99)+1]);
end

GGvsPepiruntimes=rangedata;
%% --save
save(filenameanddir,'GGvsPepiruntimes','-append');


%% GG analysis

MR = struct('Cpap', [], 'Pepi', [], 'GGpmax', [], 'Time', []);

%Pcritruns=1:length(Pcritruntimes);
MRruns=1:length(GGvsPepiruntimes);
for m = MRruns
    a = find(Time>=GGvsPepiruntimes(m,1)&Time<=GGvsPepiruntimes(m,2));
    figure(4);
    screen=get(0,'ScreenSize');
    set(4,'Position',[screen(3)*.25 39 screen(3)*.5 screen(4)-120]);
    
    dsf=1;
    X=4;
    timetemp=Time(a);
    Pepitemp=Pepi(a);
    Pmasktemp=Pmask(a);
    GGpmaxtemp=GGpmax(a);
    ax4(1)=subplot(X,1,1); plot(downsample(Time(a),dsf),downsample(EpochsXHz(a),dsf),'b',downsample(Time(a),dsf),downsample(EventsArXHz(a),dsf),'r');
    ax4(2)=subplot(X,1,2); plot(downsample(Time(a),dsf),downsample(Pmask(a),dsf));
    ax4(3)=subplot(X,1,3); plot(downsample(Time(a),dsf),downsample(Pepitemp,dsf));
    ax4(4)=subplot(X,1,4); plot(downsample(Time(a),dsf),downsample(GGpmax(a),dsf));
    
    [x1,estimatedzeroPepi,button]=ginput(1);
    
    if button==3
        continue
    end
    
    ax4(3)=subplot(X,1,3); plot(downsample(Time(a),dsf),downsample(Pepitemp,dsf),Time(a),Time(a)*0+estimatedzeroPepi,'k:');
    
    %Filter
    filter_HFcutoff_butter1 = 1;
    filter_order1 = 2;
    [B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_HFcutoff_butter1]/(1/dt/2),'low');
    Pepifilt = filtfilt(B_butterHcut,A_butterHcut,Pepitemp);
    ax4(3)=subplot(X,1,3); plot(downsample(Time(a),dsf),downsample(Pepitemp,dsf),Time(a),Pepifilt,'r',Time(a),Time(a)*0+estimatedzeroPepi,'k:');
    
    startinsp_i=[];
    endinsp_i=[];
    for i=3:length(Pepifilt)
        if (Pepifilt(i)<estimatedzeroPepi)&&(Pepifilt(i-1)>estimatedzeroPepi)
            startinsp_i=[startinsp_i;i];
        elseif (Pepifilt(i)>estimatedzeroPepi)&&(Pepifilt(i-1)<estimatedzeroPepi)
            endinsp_i=[endinsp_i;i];
        end
    end
    if startinsp_i(1)>endinsp_i(1)
        endinsp_i(1)=[];
    end
    
    clear Pepinadir_iB Pepidelta_B GGpeak_iB GGpeak_B CPAP_B
    %breath-breath values
    for i=1:(length(endinsp_i))
        [Pepidelta_B(i),Pepinadir_iB(i)] = min(Pepitemp(startinsp_i(i):endinsp_i(i)));
        Pepidelta_B(i)=Pepidelta_B(i)-estimatedzeroPepi;
        Pepinadir_iB(i)=Pepinadir_iB(i)+startinsp_i(i)-1;
        [GGpeak_B(i),GGpeak_iB(i)] = max(GGpmaxtemp(startinsp_i(i):endinsp_i(i)));
        GGpeak_iB(i)=GGpeak_iB(i)+startinsp_i(i)-1;
        CPAP_B(i) = Pmasktemp(startinsp_i(i));
    end
    Pepiabs_B=Pepidelta_B+CPAP_B;
    GGpeak_Blessmin = GGpeak_B - min(GGpmaxtemp);
    
    ax4(3)=subplot(X,1,3);
    plot(downsample(Time(a),dsf),downsample(Pepitemp,dsf),Time(a),Pepifilt,'r',Time(a),Time(a)*0+estimatedzeroPepi,'k:',timetemp(startinsp_i),Pepifilt(startinsp_i),'.r');
    hold('on')
    plot(timetemp(Pepinadir_iB),Pepitemp(Pepinadir_iB),'.b');
    hold('off')
    
    ax4(4)=subplot(X,1,4); plot(downsample(Time(a),dsf),downsample(GGpmax(a),dsf),timetemp(GGpeak_iB),GGpmaxtemp(GGpeak_iB),'.b');
    
    figure(5)
    plot(Pepiabs_B,GGpeak_B,'.k','markersize',14)
    
    MR(m).Cpap = CPAP_B;
    MR(m).Pepi = Pepiabs_B;
    MR(m).GGpmax = GGpeak_B;
    MR(m).GGpeak_Blessmin = GGpeak_Blessmin;
    MR(m).Time = timetemp(Pepinadir_iB);
    
    pause
    close all
    clear a
end
%% --save
save(filenameanddir,'MR','-append');

%% Analyse MR - GG summary
Pepiabs_B=[];
GGpeak_B=[];
Pepiabs_B_all=[];
GGpeak_B_all=[];
include_runs=1:length(GGvsPepiruntimes) %manually select which runs to include. Eg. [1]-will graph run 1 [1 2]-will graph runs 1 and 2
clear GGvsEpislope GGvsEpiCI

try
    close(5)
catch me
end

for i=include_runs
    
    tempsize=size(MR(i).Pepi);
    N=tempsize(1)*tempsize(2);
    
    if 1
        %Pepiabs_B=reshape(MR(i).Pepi,1,N);
        Pepiabs_B=reshape(MR(i).Pepi,1,N)-reshape(MR(i).Cpap,1,N);
        % GGpeak_B=reshape(MR(i).GGpmax,1,N);
        GGpeak_B=reshape(MR(i).GGpeak_Blessmin,1,N);
    else
        Pepiabs_B=[Pepiabs_B reshape(MR(i).Pepi,1,N)-reshape(MR(i).Cpap,1,N)];
        %Pepiabs_B=[Pepiabs_B reshape(MR(i).Pepi,1,N)];
        %GGpeak_B=[GGpeak_B reshape(MR(i).GGpmax,1,N)];
        GGpeak_B=[GGpeak_B reshape(MR(i).GGpeak_Blessmin,1,N)];
    end
    
    GGvsEpislope(i)=NaN;
    GGvsEpiCI(i)=NaN;
    
    if isempty(GGpeak_B)||length(GGpeak_B)<3
        continue
    end
    
    figure(5)
    
    [p] = polyfit(Pepiabs_B,GGpeak_B,1);
    cf = fit(Pepiabs_B',GGpeak_B','poly1');
    cf_confint = confint(cf)
    GGvsEpiCI(i) = (cf_confint(2,1) - cf_confint(1,1))/2;
    
    GGvsEpislope(i)=p(1);
    
    xplot=min(Pepiabs_B):0.1:max(Pepiabs_B);
    yplot=polyval(p,xplot);
    color1=rand(1,3);
    
    plot(Pepiabs_B,GGpeak_B,'.','markersize',14,'color',color1)
    hold('on')
    plot(xplot,yplot,'color',color1)
    
    % average slope
    GGvsEpislope_mean=nanmedian(GGvsEpislope);
    
    % combined data
    Pepiabs_B_all=[Pepiabs_B Pepiabs_B_all];
    GGpeak_B_all=[GGpeak_B GGpeak_B_all];
    
end

%overall fit
[p2] = polyfit(Pepiabs_B_all,GGpeak_B_all,1);

GGvsEpislope_alldatacombined = p2(1);

cf2 = fit(Pepiabs_B_all',GGpeak_B_all','poly1');
cf_confint2 = confint(cf2)
GGvsEpiCI2 = (cf_confint2(2,1) - cf_confint2(1,1))/2;

xplot=min(Pepiabs_B_all):0.1:max(Pepiabs_B_all);
yplot=polyval(p2,xplot);

plot(xplot,yplot,'k','linewidth',2)
%% --save
if saveontherun
    save(filenameanddir,'MR','GGvsEpislope_mean','GGvsEpislope_alldatacombined','-append');
end

%% Find Arousal Threshold (Pepi)

%Use to cursor to exclude (right click)
%or select Pepi breath to analyze x=time of nadir, y=baseline (start
%inspiration).

%Find arousals.
ArThresEpi.x_select=[];
ArThresEpi.y_select=[];

figure(2)
Xminute=3;
Yminute=1;
dsf=5;
X=4;
ax2(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(EpochsXHz,dsf),'b',downsample(Time,dsf),downsample(EventsArXHz,dsf),'r');
%ax2(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Pmask,dsf));
if exist('WakeSleep'), hold('on'); plot(downsample(Time,dsf),downsample(WakeSleep,dsf),'color',[1 0.3 0.3]); hold('off'); end
ax2(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Flow,dsf));
ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pepi,dsf));
if 0
    if exist('GGpmax','var')
        ax2(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(GGpmax,dsf));
    end
else
    if exist('Pes','var')
        ax2(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(Pes,dsf));
    end
end
linkaxes(ax2,'x')

%arousals are at temp==1:
temp=diff(EventsArXHz);
I=find(temp==1);
rerun=0;
for i=1:length(I)
    if rerun
        if ArThresEpi.button(i)==3
            continue
        end
    end
    lowerX=Time(I(i))-Yminute*60;
    upperX=Time(I(i))+Yminute*60;
    xlim([lowerX upperX])
    set(ax2(3),'ylim',[min(Pepi(Time>lowerX&Time<upperX)) prctile(Pepi(Time>lowerX&Time<upperX),99.5)]);
    [ArThresEpi.x_select(i),ArThresEpi.y_select(i),ArThresEpi.button(i)]=ginput(1);
end
% --save
save(filenameanddir,'ArThresEpi','-append');

%% Analyze ARthres
% no manual analysis needed

ArThresEpi.data=[];
Pepinadirabs=[];
ArThresEpi.CPAP=[]; %note this is the Pmask measured at the Pepi nadir rather than during zero flow.
ArThresEpi.data_t=[];
ArThresEpi.data_i=[];

for i=1:length(I)
    ArThresEpi.data(i)=NaN;
    Pepinadirabs(i)=NaN;
    ArThresEpi.CPAP(i)=NaN;
    ArThresEpi.data_t(i)=NaN;
    ArThresEpi.data_i(i)=NaN;
    if ArThresEpi.button(i)==3
        continue
    end
    [Pepinadir_uncorrected,index]=min(Pepi(Time>ArThresEpi.x_select(i)-1&Time<ArThresEpi.x_select(i)+1));
    ArThresEpi.data_t(i)=ArThresEpi.x_select(i)-1+(index-1)*dt;
    ArThresEpi.data_i(i)=find(Time>ArThresEpi.data_t(i),1,'first');
    %plot
    lowerX=ArThresEpi.x_select(i)-15;
    upperX=ArThresEpi.x_select(i)+15;
    ArThresEpi.data(i)=Pepinadir_uncorrected-ArThresEpi.y_select(i);
    ArThresEpi.CPAP(i)=Pmask(find(Time>ArThresEpi.data_t(i),1));
    Pepinadirabs(i)=ArThresEpi.data(i)+Pmask(ArThresEpi.data_i(i));
    %plot
    %ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pepi,dsf),downsample(Time,dsf),downsample(Time,dsf)*0+estimatedzeroPepi,'k:',ArThresEpi.data_t(i),Pepi(ArThresEpi.data_i(i)),'r.');
    ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pepi,dsf),ArThresEpi.data_t(i),Pepi(ArThresEpi.data_i(i)),'r.');
    xlim([lowerX upperX]);
    set(ax2(3),'ylim',[min(Pepi(Time>lowerX&Time<upperX)) prctile(Pepi(Time>lowerX&Time<upperX),95)]);
    hold('off');
end

%% Summary analysis
ArThresEpi.mean = nanmean(ArThresEpi.data);
ArThresEpi.median = nanmedian(ArThresEpi.data);
ArThresEpi.N = sum(~isnan(ArThresEpi.data));
ArThresEpi.SEM = nanstd(ArThresEpi.data)/ArThresEpi.N^0.5;

ArThresEpi.mean_abs = nanmean(Pepinadirabs);
ArThresEpi.median_abs = nanmedian(Pepinadirabs);
ArThresEpi.SEM_abs = nanstd(Pepinadirabs)/ArThresEpi.N^0.5;

ArThresEpi.cov = nanstd(ArThresEpi.data)/-ArThresEpi.mean;

%% --save
save(filenameanddir,'ArThresEpi','-append');

%% Find Arousal Threshold (Pes)

%Use to cursor to exclude (right click)
%or select Pepi breath to analyze x=time of nadir, y=baseline (start
%inspiration).

global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Flow','Pes'};
PlotAndSelectData(signallist,Time,clearvalues);
It = Time(find(diff(EventsArXHz)==1));


%% Analyze ARthres Pes
% no manual analysis needed
ArThresPes.x_select = xvalues(:,1);
ArThresPes.y_select = yvalues(:,1);
%%
ArThresPes.data=[];
Pesnadirabs=[];
ArThresPes.CPAP=[]; %note this is the Pmask measured at the Pes nadir rather than during zero flow.
ArThresPes.data_t=[];
ArThresPes.data_i=[];

for i=1:length(ArThresPes.x_select)
    %     ArThresPes.data(i)=NaN;
    %     Pesnadirabs(i)=NaN;
    %     ArThresPes.CPAP(i)=NaN;
    %     ArThresPes.data_t(i)=NaN;
    %     ArThresPes.data_i(i)=NaN;
    %     if ArThresPes.button(i)==3
    %         continue
    %     end
    [Pesnadir_uncorrected,index]=min(Pes(Time>ArThresPes.x_select(i)-1&Time<ArThresPes.x_select(i)+1));
    ArThresPes.data_t(i)=ArThresPes.x_select(i)-1+(index-1)*dt;
    ArThresPes.data_i(i)=find(Time>ArThresPes.data_t(i),1,'first');
    %plot
    lowerX=ArThresPes.x_select(i)-15;
    upperX=ArThresPes.x_select(i)+15;
    ArThresPes.data(i)=Pesnadir_uncorrected-ArThresPes.y_select(i);
    ArThresPes.CPAP(i)=Pmask(find(Time>ArThresPes.data_t(i),1));
    Pesnadirabs(i)=ArThresPes.data(i)+Pmask(ArThresPes.data_i(i));
    %plot
    %ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pepi,dsf),downsample(Time,dsf),downsample(Time,dsf)*0+estimatedzeroPepi,'k:',ArThresEpi.data_t(i),Pepi(ArThresEpi.data_i(i)),'r.');
    ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pes,dsf),ArThresPes.data_t(i),Pes(ArThresPes.data_i(i)),'r.');
    xlim([lowerX upperX]);
    set(ax2(3),'ylim',[min(Pes(Time>lowerX&Time<upperX)) prctile(Pes(Time>lowerX&Time<upperX),95)]);
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

%% Find Veupneas
global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Flow','Pes'};
PlotAndSelectData(signallist,Time,clearvalues);


%%
Veupnea.ranges=xvalues;

if 0&saveontherun
    save(filenameanddir,'Veupnea','-append');
end
%% Calculate Veupnea
for i=size(Veupnea.ranges,1):-1:1
    if Veupnea.ranges(i,2)<Veupnea.ranges(i,1)
        Veupnea.ranges(i,:)=[];
    end
end
clear Veupnea_n Veupnea_Ttot_n T0 Veupnea_Vpeak VIlist
VIlist=[];
minimum_figs=0;
for i=1:size(Veupnea.ranges,1)
    try
        irange=(Time>Veupnea.ranges(i,1)-2)&(Time<=Veupnea.ranges(i,2)+2);
        %[Ix,VI,VTi,Ti,Te]=Vflowanalysis1(Flow(irange),Time(irange),dt,minimum_figs);
        [~,~,~,~,~,~,VI,~,Ttot,~,~,~,Vpeak,~]=Vflowanalysis3(Time(irange),Flow(irange),0);
        T0(i)=sum(Ttot);
        VIlist = [VIlist VI];
        Veupnea_n(i)=sum(VI.*Ttot)/T0(i);
        Veupnea_Ttot_n(i)=mean(Ttot);
        Veupnea_Vpeak(i)=mean(Vpeak);
    catch me
        T0(i)=0;
        Veupnea_n(i)=NaN;
        Veupnea_Ttot_n(i)=NaN;
        Veupnea_Vpeak(i)=NaN;
    end
end

Veupnea.mean=nansum(Veupnea_n.*T0)/sum(T0)*60;
Veupnea.duration=sum(T0)/60;
Veupnea.cov=nanstd(Veupnea_n)/nanmean(Veupnea_n);
Veupnea.meanVpeak=nanmean(Veupnea_n.*T0)/sum(T0)*60;
Veupnea.cov_breathtobreath=nanstd(VIlist*60)/Veupnea.mean;

%% Save Veupnea
save(filenameanddir,'Veupnea','-append');

%% Vpassive

Vpassive.Feupnea = Pcrit.Vcritvalue/Veupnea.mean;
Vpassive.Direct_Feupnea_median = median(Pcrit.VE(Pcrit.CPAP>-1&Pcrit.CPAP<1))/Veupnea.mean;
Vpassive.Direct_Feupnea_N = length(Pcrit.CPAP>-1&Pcrit.CPAP<1);

%%
% save(filenameanddir,'Vpassive','-append');
save(filenameanddir,'Vpassive','Veupnea','-append');

%% Test
% Select left and right of ranges. Aim for times with mean PETO2=100mmHg
global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Pes','Flow'};
PlotAndSelectData(signallist,Time,clearvalues);


%% SpO2 stable wake
% Select left and right of ranges. Aim for times with mean PETO2=100mmHg
global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Flow','PO2','SaO2'};
PlotAndSelectData(signallist,Time,clearvalues);


%%
SpO2wake.ranges=xvalues;

%make larger array of appended SpO2 data
SpO2appendeddata=[];
for i=1:size(SpO2wake.ranges,1)
    irange=Time>SpO2wake.ranges(i,1)&Time<=SpO2wake.ranges(i,2);
    SpO2appendeddata = [SpO2appendeddata;round(SaO2(irange))];
end
SpO2appendeddata(SpO2appendeddata>100)=100;

SpO2wake.median=median(SpO2appendeddata);
SpO2wake.mean=mean(SpO2appendeddata);

clear PpO2
deltaSpO2=-3:0.1:-(100-max(SpO2appendeddata));
for i=1:length(deltaSpO2)
    PpO2(i) = mean(S_O2toP_O2(1,SpO2appendeddata+deltaSpO2(i),23400,150));
end
deltaSpO2(isnan(PpO2)|PpO2>150)=[];
PpO2(isnan(PpO2)|PpO2>150)=[];

SpO2wake.deltaSpO2forPO2of100mmHg=interp1(PpO2,deltaSpO2,100,'spline');
figure(); plot(deltaSpO2,PpO2,SpO2wake.deltaSpO2forPO2of100mmHg,100,'r.'); ylabel('PO2'); xlabel('delta SpO2');
SpO2_test = [100 99 98 97];
PO2_test = S_O2toP_O2(1,SpO2_test+SpO2wake.deltaSpO2forPO2of100mmHg,23400,150);

if saveontherun
    save(filenameanddir,'SpO2wake','-append');
end

%% CCWCPAP: ChestWallCompliance CPAP using Pcrit drops

%errors: deltavol: [NaN NaN etc] and deltaCPAP==deltaPes
CCWCPAP.runtimes = Pcrit.runtimes;
%CCWCPAP.runtimes(Pcrit.Peakflow==0,:)=[]; %deletes Pcrit drops with zero flow

%instructions:
%1. using two clicks, select left and right of breaths ON CPAP (select just OUTSIDE the desired breaths)
%---aim for last 3 full breaths on CPAP and skip first 2 full breaths off CPAP
%1b then same on lower CPAP
%2. zoom in Y axis on Pes using two clicks upper and lower Pes limit (first zoom click right click to indicate poor deltavol
%3. select Pes at zero flow on high CPAP and low CPAP (first select click use right click to indicate poor deltavol

%3rd click right click for exclude delta vol
%5th click right click for exclude delta pes

%right click for do not analyze
%clear button buttonPes deltavol deltaCPAP deltaPes
for i=[1:size(CCWCPAP.runtimes,1)]
    try
        close 3
    catch me
    end
    try
        deltaX=45;
        lowerX=CCWCPAP.runtimes(i,1)-deltaX;
        upperX=CCWCPAP.runtimes(i,2)+20;
        I=find(Time>lowerX&Time<upperX);
        deltaX=6;
        if Pcrit.Peakflow(i)==0
            deltaX=0;
        end
        lowerX=CCWCPAP.runtimes(i,1)-deltaX;
        upperX=CCWCPAP.runtimes(i,2)+deltaX;
        
        
        vol1_=cumsum(Flow(I)-mean(Flow(I)))*dt;
        figure(3);
        ax3(1)=subplot(4,1,1); plot(Time(I),Flow(I),'k',[CCWCPAP.runtimes(i,1) CCWCPAP.runtimes(i,1)],[-1 1],'r',[CCWCPAP.runtimes(i,2) CCWCPAP.runtimes(i,2)],[-1 1],'r');
        hold('on');
        ax3(2)=subplot(4,1,2); plot(Time(I),vol1_,'k',[CCWCPAP.runtimes(i,1) CCWCPAP.runtimes(i,1)],[-1 1],'r',[CCWCPAP.runtimes(i,2) CCWCPAP.runtimes(i,2)],[-1 1],'r');
        hold('on');
        ax3(3)=subplot(4,1,3); plot(Time(I),Pmask(I),'k',[CCWCPAP.runtimes(i,1) CCWCPAP.runtimes(i,1)],[-1 1],'r',[CCWCPAP.runtimes(i,2) CCWCPAP.runtimes(i,2)],[-1 1],'r');
        hold('on');
        ax3(4)=subplot(4,1,4); plot(Time(I),Pes(I),'k',[CCWCPAP.runtimes(i,1) CCWCPAP.runtimes(i,1)],[-1 1],'r',[CCWCPAP.runtimes(i,2) CCWCPAP.runtimes(i,2)],[-1 1],'r');
        hold('on');
        linkaxes(ax3,'x');
        %first get flow indices at baseline
        [xtemp,temp]=ginput(2);
        I2=find(Time>xtemp(1)&Time<xtemp(2));
        
        %then get the flow indices during CPAP drop
        
        %find the max Pes values and auto zoom y axis
        
        vol1=cumsum(Flow(I2)-mean(Flow(I2)))*dt;
        
        
        filter_HFcutoff_butter1 = 2;
        filter_order1 = 2;
        [B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_HFcutoff_butter1]/(Fs/2),'low');
        vol1 = filtfilt(B_butterHcut,A_butterHcut,vol1);
        
        thresh=std(vol1)/10;
        
        if Pcrit.Peakflow(i)~=0
            
            [min_list, max_list] = peakdet(-vol1,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
            
        else
            tempmini = [1:Fs:(length(vol1)-round(Fs/2)) length(vol1)]';
            min_list = [tempmini vol1(tempmini)];
            tempmaxi = tempmini(1:end-1)+round(diff(tempmini)/2);
            max_list = [tempmaxi vol1(tempmaxi)];
        end
        %[indVI, indVE] = CheckVind(max_list, min_list);
        %TidalVol = V(indVI)-V(indVE);
        ax3(2)=subplot(4,1,2); plot(Time(I),vol1_,'k',[CCWCPAP.runtimes(i,1) CCWCPAP.runtimes(i,1)],[-1 1],'r',[CCWCPAP.runtimes(i,2) CCWCPAP.runtimes(i,2)],[-1 1],'r');
        
        plot(Time(I2(1)-1+max_list(:,1)),vol1_(I2(1)-I(1)+max_list(:,1)),'k.');
        plot(Time(I2(1)-1+min_list(:,1)),vol1_(I2(1)-I(1)+min_list(:,1)),'r.');
        %[ymin,temptemp]=get(ax2(4),'ylim')
        linkaxes(ax3,'x');
        
        lefti=find(Time(I2(1)-1+min_list(:,1))>xtemp(1),1,'first');
        righti=find(Time(I2(1)-1+min_list(:,1))<xtemp(2),1,'last')-1;
        if isempty(lefti)
            lefti=1;
        end
        if isempty(righti)
            righti=lefti;
        end
        if righti<lefti
            righti=lefti;
        end
        ax3(2)=subplot(4,1,2); plot([Time(I2(1)-1+min_list(lefti,1)) Time(I2(1)-1+min_list(lefti,1))],[-1 1],'r:',[Time(I2(1)-1+min_list(righti+1,1)) Time(I2(1)-1+min_list(righti+1,1))],[-1 1],'r:');
        
        xpoly = Time(I2(1)-1+min_list(lefti:(righti+1),1));
        ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
        ypoly = vol1(min_list(lefti:(righti+1),1));
        ppoly = polyfit(xpoly,ypoly,1);
        Yvoldrift = polyval(ppoly,Time(I2));
        
        if Pcrit.Peakflow(i)~=0
            %repeat find end-exp and end-insp now that we have done some leak correction:
            vol3=cumsum(Flow(I2)-mean(Flow(I2))-ppoly(1))*dt;
            vol3 = filtfilt(B_butterHcut,A_butterHcut,vol3);
            thresh=std(vol3)/10;
            [min_list, max_list] = peakdet(-vol3,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
            
            plot(Time(I2(1)-1+max_list(:,1)),vol1_(I2(1)-I(1)+max_list(:,1)),'ko');
            plot(Time(I2(1)-1+min_list(:,1)),vol1_(I2(1)-I(1)+min_list(:,1)),'ro');
            
            lefti=find(Time(I2(1)-1+min_list(:,1))>xtemp(1),1,'first');
            righti=find(Time(I2(1)-1+min_list(:,1))<xtemp(2),1,'last')-1;
            if isempty(lefti)
                lefti=1;
            end
            if isempty(righti)||righti==0
                righti=size(min_list,1)-1;
            end
            plot([Time(I2(1)-1+min_list(lefti,1)) Time(I2(1)-1+min_list(lefti,1))],[-1 1],'r:',[Time(I2(1)-1+min_list(righti+1,1)) Time(I2(1)-1+min_list(righti+1,1))],[-1 1],'r:');
            
            xpoly = Time(I2(1)-1+min_list(lefti:(righti+1),1));
            ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
            ypoly = vol1(min_list(lefti:(righti+1),1));
            ppoly = polyfit(xpoly,ypoly,1);
            Yvoldrift = polyval(ppoly,Time(I2));
            
            %     else
        end
        ypoly_ = vol1_(I2(1)-I(1)+min_list(lefti:(righti+1),1));
        ppoly_ = polyfit(xpoly,ypoly_,1);
        Yvoldrift_ = polyval(ppoly_,Time(I2));
        
        
        Yvol1_exp = ypoly(1:end-1); %end exp
        Yvol1_insp = vol1(max_list(lefti:(righti),1)); %end insp
        Yvoldrift_exp = polyval(ppoly,xpoly); Yvoldrift_exp(end)=[];
        Yvoldrift_insp = polyval(ppoly,Time(I2(1)-1+max_list(lefti:(righti),1)));
        
        plot(Time(I2),Yvoldrift_);
        %    hold('off'); hold('on');
        
        ax3(3)=subplot(4,1,3); plot(Time(I),Pmask(I)); hold('on');
        plot(Time(I2(1)-1+max_list(:,1)),Pmask(I2(1)-1+max_list(:,1)),'k.');
        plot(Time(I2(1)-1+min_list(:,1)),Pmask(I2(1)-1+min_list(:,1)),'r.');
        ax3(4)=subplot(4,1,4); plot(Time(I),Pes(I)); hold('on');
        plot(Time(I2(1)-1+max_list(:,1)),Pes(I2(1)-1+max_list(:,1)),'k.');
        plot(Time(I2(1)-1+min_list(:,1)),Pes(I2(1)-1+min_list(:,1)),'r.');
        linkaxes(ax3,'x');
        
        endexpvol=Yvol1_exp-Yvoldrift_exp;
        endinspvol=Yvol1_insp-Yvoldrift_insp;
        volbreath=endinspvol-endexpvol;
        
        [CCWCPAP.runtimes2(i,1),~,button(i)] = ginput(1);
        %     if button(i)==3
        %         continue
        %     end
        [CCWCPAP.runtimes2(i,2 ),~] = ginput(1);
        lowerX=CCWCPAP.runtimes2(i,1);
        upperX=CCWCPAP.runtimes2(i,2);
        I4=find(Time>lowerX&Time<upperX);
        
        vol4=cumsum(Flow(I4)-mean(Flow(I4)))*dt;
        vol4 = filtfilt(B_butterHcut,A_butterHcut,vol4);
        thresh=std(vol4)/10;
        [min_list2, max_list2] = peakdet(-vol4,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
        ax3(4)=subplot(4,1,2);
        plot(Time(I4(1)-1+max_list2(:,1)),vol1_(I4(1)-I(1)+max_list2(:,1)),'k.');
        plot(Time(I4(1)-1+min_list2(:,1)),vol1_(I4(1)-I(1)+min_list2(:,1)),'r.');
        %[ymin,temptemp]=get(ax2(4),'ylim')
        linkaxes(ax3,'x');
        
        lefti=1;%find(Time(I4(1)-1+min_list2(:,1))<CCWCPAP.runtimes2(i,1),1,'last');
        righti=length(min_list2)-1;%find(Time(I4(1)-1+min_list2(:,1))>CCWCPAP.runtimes2(i,2),1,'first')-1;
        ax3(2)=subplot(4,1,2); plot([Time(I4(1)-1+min_list2(lefti,1)) Time(I4(1)-1+min_list2(lefti,1))],[-1 1],'r:',[Time(I4(1)-1+min_list2(righti+1,1)) Time(I4(1)-1+min_list2(righti+1,1))],[-1 1],'r:');
        
        xpoly2 = Time(I4(1)-1+min_list2(lefti:(righti+1),1));
        ttot = (xpoly2(end)-xpoly2(1))/(size(xpoly2,1)-1);
        ypoly2 = vol4(min_list2(lefti:(righti+1),1));
        ppoly2 = polyfit(xpoly2,ypoly2,1);
        Yvoldrift2 = polyval(ppoly2,Time(I4));
        
        %repeat find end-exp and end-insp now that we have done some leak correction:
        vol5=cumsum(Flow(I4)-mean(Flow(I4))-ppoly2(1))*dt;
        vol5 = filtfilt(B_butterHcut,A_butterHcut,vol5);
        thresh=std(vol4)/10;
        [min_list2, max_list2] = peakdet(-vol4,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
        
        plot(Time(I4(1)-1+max_list2(:,1)),vol1_(I4(1)-I(1)+max_list2(:,1)),'ko');
        plot(Time(I4(1)-1+min_list2(:,1)),vol1_(I4(1)-I(1)+min_list2(:,1)),'ro');
        
        lefti=1;%find(Time(I4(1)-1+min_list2(:,1))<CCWCPAP.runtimes2(i,1),1,'last');
        righti=length(min_list2)-1;%find(Time(I4(1)-1+min_list2(:,1))>CCWCPAP.runtimes2(i,2),1,'first')-1;
        plot([Time(I4(1)-1+min_list2(lefti,1)) Time(I4(1)-1+min_list2(lefti,1))],[-1 1],'r:',[Time(I4(1)-1+min_list2(righti+1,1)) Time(I4(1)-1+min_list2(righti+1,1))],[-1 1],'r:');
        
        xpoly2 = Time(I4(1)-1+min_list2(lefti:(righti+1),1));
        ttot = (xpoly2(end)-xpoly2(1))/(size(xpoly2,1)-1);
        ypoly2 = vol4(min_list2(lefti:(righti+1),1));
        ppoly2 = polyfit(xpoly2,ypoly2,1);
        Yvoldrift2 = polyval(ppoly2,Time(I4));
        
        ypoly2_ = vol1_(I4(1)-I(1)+min_list2(lefti:(righti+1),1));
        ppoly2_ = polyfit(xpoly2,ypoly2_,1);
        Yvoldrift2_ = polyval(ppoly2_,Time(I4));
        
        Yvol1_exp = ypoly2(1:end-1); %end exp
        Yvol1_insp = vol5(max_list2(lefti:(righti),1)); %end insp
        Yvoldrift_exp = polyval(ppoly2,xpoly2); Yvoldrift_exp(end)=[];
        Yvoldrift_insp = polyval(ppoly2,Time(I4(1)-1+max_list2(lefti:(righti),1)));
        
        
        
        subplot(4,1,2);
        plot(Time(I4),Yvoldrift2_);
        %    hold('off'); hold('on');
        
        ax3(3)=subplot(4,1,3); %plot(Time(I),Pmask(I)); hold('on');
        plot(Time(I2(1)-1+max_list(:,1)),Pmask(I2(1)-1+max_list(:,1)),'k.');
        plot(Time(I2(1)-1+min_list(:,1)),Pmask(I2(1)-1+min_list(:,1)),'r.');
        plot(Time(I4(1)-1+max_list2(:,1)),Pmask(I4(1)-1+max_list2(:,1)),'k.');
        plot(Time(I4(1)-1+min_list2(:,1)),Pmask(I4(1)-1+min_list2(:,1)),'r.');
        ax3(4)=subplot(4,1,4); plot(Time(I),Pes(I)); hold('on');
        plot(Time(I2(1)-1+max_list(:,1)),Pes(I2(1)-1+max_list(:,1)),'k.');
        plot(Time(I2(1)-1+min_list(:,1)),Pes(I2(1)-1+min_list(:,1)),'r.');
        plot(Time(I4(1)-1+max_list2(:,1)),Pes(I4(1)-1+max_list2(:,1)),'k.');
        plot(Time(I4(1)-1+min_list2(:,1)),Pes(I4(1)-1+min_list2(:,1)),'r.');
        linkaxes(ax3,'x');
        
        %assuming linear relationship between leak and CPAP level
        CCWCPAP.CPAPright(i) = mean(Pmask(I2(1)-1+max_list(:,1)));
        CCWCPAP.CPAPleft(i) = mean(Pmask(I4(1)-1+max_list2(:,1)));
        leak(1)=ppoly_(1);
        leak(2)=ppoly2_(1);
        ppolyleak = polyfit([CCWCPAP.CPAPright(i) CCWCPAP.CPAPleft(i)],leak,1);
        Leak_t = polyval(ppolyleak,Pmask(I))+mean(Flow(I));
        
        ax3(1)=subplot(4,1,1); %
        plot(Time(I),Flow(I)-Leak_t,'r:'); %hold('on');
        vol6 = cumsum(Flow(I)-Leak_t)*dt;
        vol6 = vol6 - mean(vol6(I4(1)-I(1)+min_list2(:,1)));
        deltavol(i) = -mean(vol6(I2(1)-I(1)+min_list(:,1)));
        deltaCPAP(i) = CCWCPAP.CPAPleft(i)-CCWCPAP.CPAPright(i);
        
        ax3(2)=subplot(4,1,2); %
        plot(Time(I),vol6,'r:');
        
        
        [~,Y(1),buttontemp]=ginput(1); %can still right click after the fact, but not needed (click 3)
        if buttontemp==3
            button(i)=3;
        end
        [~,Y(2)]=ginput(1); %click 4
        ax3(4)=subplot(4,1,4);
        ylim([min(Y) max(Y)]);
        [~,PesResting(1),buttonPes(i)]=ginput(1); %click 5 determines whether to include Pes data
        if buttonPes(i)==1
            [~,PesResting(2)]=ginput(1); %get a click 6 if including Pes
            deltaPes(i) = PesResting(1)-PesResting(2);
        else
            deltaPes(i) = NaN;
        end
    catch me
        button(i)=NaN;
        CCWCPAP.deltavol(i)=NaN;
        CCWCPAP.deltaCPAP(i)=NaN;
        CCWCPAP.deltaPes(i)=NaN;
        CCWCPAP.CPAPleft(i)=NaN;
        CCWCPAP.CPAPright(i)=NaN;
    end
end

CCWCPAP.button = button;
CCWCPAP.buttonPes = buttonPes;
CCWCPAP.deltavol = deltavol;
CCWCPAP.deltaCPAP = deltaCPAP;
CCWCPAP.deltaPes = -deltaPes;

CCWCPAP.deltavol(CCWCPAP.deltavol>=0)=NaN;
CCWCPAP.deltavol(CCWCPAP.button==3)=NaN;
CCWCPAP.deltavol(CCWCPAP.button==0)=NaN;

CCWCPAP.deltaCPAP(CCWCPAP.deltaCPAP==0) = NaN;
CCWCPAP.deltaPes(CCWCPAP.deltaPes==0)=NaN;
CCWCPAP.deltaPes(CCWCPAP.deltaPes>CCWCPAP.deltaCPAP)=CCWCPAP.deltaPes(CCWCPAP.deltaPes>CCWCPAP.deltaCPAP);
CCWCPAP.Ecw_Ftotal = nanmedian(CCWCPAP.deltaPes./CCWCPAP.deltaCPAP);

CCWCPAP.Crs = nanmedian(CCWCPAP.deltavol./CCWCPAP.deltaCPAP);

CCWCPAP.CrsStd = nanstd(CCWCPAP.deltavol./CCWCPAP.deltaCPAP);

CCWCPAP.Ecw_FtotalStd = nanstd(CCWCPAP.deltaPes./CCWCPAP.deltaCPAP);

Ers = 1/CCWCPAP.Crs;
Ecw = CCWCPAP.Ecw_Ftotal*Ers;
CCWCPAP.Ccw = 1/Ecw;
CCWCPAP.Clung = 1/(Ers-Ecw);

figure(1001);
subplot(2,1,1);plot(deltaCPAP,CCWCPAP.deltavol./CCWCPAP.deltaCPAP,'.'); ylabel('dVol/dCPAP');
subplot(2,1,2);plot(deltaCPAP,CCWCPAP.deltaPes./CCWCPAP.deltaCPAP,'.'); ylabel('dPes/dCPAP');

%% --save
if saveontherun
    save(filenameanddir,'CCWCPAP','-append');
end

%% Find Varousal

%%new

global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
Ptrend = filtfilt(B_butterHcut,A_butterHcut,Pes);
signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Flow','[Pes Ptrend-ArThresPes.median*1/3 Ptrend+ArThresPes.median*2/3]'};
PlotAndSelectData(signallist,Time,clearvalues);

%% --save
Varousal.ranges=xvalues;

if saveontherun
    save(filenameanddir,'Varousal','-append');
end

%% Analyze Varousal
%load(filenameanddir,'VARx1');

%idea for getting Varousal from spontaneous breathing (and improving the gold standard)
%find 8 breath sequences during sleep, and if relatively stable [criteria], use survival analysis to determine the minimum tolerable ventilation.



for i=1:size(Varousal.ranges,1)
    %try
    
    deltaX=6;
    lowerX=Varousal.ranges(i,1)-deltaX;
    upperX=Varousal.ranges(i,2)+deltaX;
    figure(3)
    I=find(Time>lowerX&Time<upperX);
    subplot(2,1,1); plot(Time(I),Flow(I),'k',[Varousal.ranges(i,1) Varousal.ranges(i,1)],[-1 1],'r',[Varousal.ranges(i,2) Varousal.ranges(i,2)],[-1 1],'r');
    vol1=cumsum(Flow(I)-mean(Flow(I)))*dt;
    title(i);
    pause(0.33);
    thresh=std(vol1)/10;
    [min_list, max_list] = peakdet(-vol1,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
    %[indVI, indVE] = CheckVind(max_list, min_list);
    %TidalVol = V(indVI)-V(indVE);
    subplot(2,1,2); plot(Time(I),vol1,'k',[Varousal.ranges(i,1) Varousal.ranges(i,1)],[-1 1],'r',[Varousal.ranges(i,2) Varousal.ranges(i,2)],[-1 1],'r');
    hold('on')
    subplot(2,1,2); plot(Time(I(1)-1+max_list(:,1)),vol1(max_list(:,1)),'k.');
    subplot(2,1,2); plot(Time(I(1)-1+min_list(:,1)),vol1(min_list(:,1)),'r.');
    %[ymin,temptemp]=get(ax2(4),'ylim')
    
    lefti=find(Time(I(1)-1+min_list(:,1))<Varousal.ranges(i,1),1,'last');
    righti=find(Time(I(1)-1+min_list(:,1))>Varousal.ranges(i,2),1,'first')-1;
    if isempty(lefti)
        lefti=1;
    end
    if isempty(righti)
        righti=lefti;
    end
    if righti<lefti
        righti=lefti;
    end
    subplot(2,1,2); plot([Time(I(1)-1+min_list(lefti,1)) Time(I(1)-1+min_list(lefti,1))],[-1 1],'r:',[Time(I(1)-1+min_list(righti+1,1)) Time(I(1)-1+min_list(righti+1,1))],[-1 1],'r:');
    
    xpoly = Time(I(1)-1+min_list(lefti:(righti+1),1));
    ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
    ypoly = vol1(min_list(lefti:(righti+1),1));
    ppoly = polyfit(xpoly,ypoly,1);
    Yvoldrift = polyval(ppoly,Time(I));
    
    %repeat find end-exp and end-insp now that we have done some leak correction:
    vol3=cumsum(Flow(I)-mean(Flow(I))-ppoly(1))*dt;
    thresh=std(vol3)/10;
    [min_list, max_list] = peakdet(-vol3,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
    
    subplot(2,1,2); plot(Time(I(1)-1+max_list(:,1)),vol1(max_list(:,1)),'ko');
    subplot(2,1,2); plot(Time(I(1)-1+min_list(:,1)),vol1(min_list(:,1)),'ro');
    
    lefti=find(Time(I(1)-1+min_list(:,1))<Varousal.ranges(i,1),1,'last');
    righti=find(Time(I(1)-1+min_list(:,1))>Varousal.ranges(i,2),1,'first')-1;
    if isempty(lefti)
        lefti=1;
    end
    if isempty(righti)||righti==0
        righti=size(min_list,1)-1;
    end
    subplot(2,1,2); plot([Time(I(1)-1+min_list(lefti,1)) Time(I(1)-1+min_list(lefti,1))],[-1 1],'r:',[Time(I(1)-1+min_list(righti+1,1)) Time(I(1)-1+min_list(righti+1,1))],[-1 1],'r:');
    
    xpoly = Time(I(1)-1+min_list(lefti:(righti+1),1));
    ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
    ypoly = vol1(min_list(lefti:(righti+1),1));
    ppoly = polyfit(xpoly,ypoly,1);
    Yvoldrift = polyval(ppoly,Time(I));
    
    Yvol1_exp = ypoly(1:end-1); %end exp
    Yvol1_insp = vol1(max_list(lefti:(righti),1)); %end insp
    Yvoldrift_exp = polyval(ppoly,xpoly); Yvoldrift_exp(end)=[];
    Yvoldrift_insp = polyval(ppoly,Time(I(1)-1+max_list(lefti:(righti),1)));
    
    subplot(2,1,2); plot(Time(I),Yvoldrift);
    hold('off');
    
    endexpvol=Yvol1_exp-Yvoldrift_exp;
    endinspvol=Yvol1_insp-Yvoldrift_insp;
    volbreath=endinspvol-endexpvol;
    Varousal.data(i)=mean(volbreath)/ttot*60;
    Varousal.CPAP(i)=mean(Pmask(I(1)-1+min_list(lefti:(righti+1),1)));
    %     catch me
    %     end
end

Varousal.state=EpochsXHz(round((Varousal.ranges(:,1)-Time(1))/dt+1));
criteria = (Varousal.state<5);

Varousal.mean = mean(Varousal.data(criteria));
Varousal.cov = std(Varousal.data(criteria))/mean(Varousal.data(criteria));
Varousal.median = median(Varousal.data(criteria));
Varousal.Feupnea = Varousal.median/Veupnea.mean;

Varousal.Feupnea_onCPAP = median(Varousal.data(criteria'&(Varousal.CPAP>1|Varousal.CPAP<-1)))/Veupnea.mean;
Varousal.Feupnea_onCPAP_N = sum(criteria'&(Varousal.CPAP>1|Varousal.CPAP<-1));
Varousal.Feupnea_offCPAP = median(Varousal.data(criteria'&(Varousal.CPAP<=1&Varousal.CPAP>=-1)))/Veupnea.mean;
Varousal.Feupnea_offCPAP_N = sum(criteria'&(Varousal.CPAP<=1&Varousal.CPAP>=-1));

%% --save
if saveontherun
    save(filenameanddir,'Varousal','-append');
end

if saveontherun
    save(filenameanddir,'Pcrit','Veupnea','ArThresPes','Varousal','CCWCPAP','SpO2wake','-append');
end

%% pressure-VE graph on top of passive pressure-flow plots

figure(4);
set(4,'color',[1 1 1]);
plotpoints=0;
ax2(1)=subplot(1,1,1);
[~,Vcritvalue,~,Pcritvalue,PVSlope]=pcrit_model_run(Pmaskdata,VdotEdata,0,plotpoints) %Pmaskdata VdotMaxdata VdotEdata
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('VE','fontname','arial narrow','fontsize',12);

plot(Varousal.CPAP,Varousal.data,'r.')

global fixedslope
fixedslope=NaN;

%% Find Vactive data

global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
if 1
    try
        newsignal = 0*Flow;
        for i=1:size(Vactive.ranges,1)
            newsignal(Time>Vactive.ranges(i,1)&Time<Vactive.ranges(i,2))=1;
        end
        signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','[Flow newsignal]','Pes','Edi'};
    catch me
    end
else
    signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Flow','Pes','Edi'};
end

PlotAndSelectData(signallist,Time,clearvalues);

%%

Vactive.ranges=xvalues;

if 0&saveontherun
    save(filenameanddir,'Vactive','-append');
end

%% Analyze Vactive

%left click ok | right click if zero flow | ...
%right click left of screen to exclude | left click left of screen to go back

i=1;
exclude=[];
Vactive.Peakflow_list=[];
Vactive.Vactive_list=[];
Vactive.CPAP_list=[];
Vactive.State=[];
%Vactive.Position=[];

deltaX=10;
while i<=size(Vactive.ranges,1)
    lowerX=Vactive.ranges(i,1)-deltaX;
    upperX=Vactive.ranges(i,2)+deltaX;
    I=find(Time>lowerX&Time<upperX);
    
    [~,~,BB_i_start,~,BB_i_end,BB_t,VI,VE,Ttot,leak,IEratio,VT,Vpeak,Vpeakmean]=Vflowanalysis3(Time(I),Flow(I),0);
    I2=Time((I(1)-1+BB_i_end))>Vactive.ranges(i,1)&Time((I(1)-1+BB_i_start))<Vactive.ranges(i,2)
    hold('on')
    plot([Vactive.ranges(i,1) Vactive.ranges(i,1)],[-1 1],'r')
    plot([Vactive.ranges(i,2) Vactive.ranges(i,2)],[-1 1],'r')
    hold('off');
    Vactive.Peakflow_list(i)=median(Vpeak(I2));
    Vactive.Vactive_list(i)=median(VI(I2))*60;
    Vactive.CPAP_list(i)=median(Pmask(I(1)-1+BB_i_start(I2==1)));
    Vactive.State(i)=EpochsXHz(I(1));
    %Vactive.State(i)=Position(I(1));
    exclude(i)=0;
    [x1,~,button(i)] = ginput(1);
    if x1<Time(I(1))&&button(i)==3 %exclude for other reason
        exclude(i)=1;
    elseif button(i)==3 %zero flow
        Vactive.Peakflow_list(i)=0;
        Vactive.Vactive_list(i)=0;
        Vactive.CPAP_list(i)=Pmask(find(Time>Vactive.ranges(i,1),1));
        exclude(i)=0;
    elseif x1<Time(I(1))&&button(i)==1&&i>1 %go back
        i=i-2;
    end
    i=i+1;
end
%%
exclude(isnan(Vactive.Vactive_list))=1;
exclude(Vactive.State==5)=1;

Vactive.Vactive_list(exclude==1)=[];
Vactive.CPAP_list(exclude==1)=[];
Vactive.State(exclude==1)=[];
Vactive.Peakflow_list(exclude==1)=[];


%% Save Vactive
if 0&saveontherun
    save(filenameanddir,'Vactive','-append');
end

%% Draw active and passive pressure-VE plots
global fixedslope
maxPdata =  99;
figure(4);
set(4,'color',[1 1 1]);
plotpoints=2;
ax2(1)=subplot(1,1,1);


if isfield(Pcrit,'forcedslopes')&&~isnan(Pcrit.forcedslopes(1))
    fixedslope = Pcrit.forcedslopes(2);
else
    fixedslope=NaN;
end

if 1
    Pmaskdata=Pcrit.CPAP(Pcrit.include_data);
    VdotEdata=Pcrit.VE(Pcrit.include_data);
end
[~,Vcritvalue,~,Pcritvalue,PVSlope]=pcrit_model_run(Pmaskdata,VdotEdata,0,plotpoints) %Pmaskdata VdotMaxdata VdotEdata
xlabel('Pmask','fontname','arial narrow','fontsize',12);
ylabel('VE','fontname','arial narrow','fontsize',12);
%plot(Varousal_CPAP,Varousal,'r.')

plot(Vactive.CPAP_list,Vactive.Vactive_list,'k.','markersize',18)



Pdata=Vactive.CPAP_list;
Vdata=Vactive.Vactive_list;
Vdata(Pdata>maxPdata)=[];
Pdata(Pdata>maxPdata)=[];

fixedslope=Pcrit.PVSlope_VE;
%fixedslope=NaN;
[VcritvalueSEMA,VcritvalueA,PcritvalueSEMA,PcritvalueA,PVSlopeA]=pcrit_model_run(Pdata,Vdata,0,1);

Vactive.Vcrit_Feupnea = VcritvalueA/Veupnea.mean;
Vactive.Vcrit_Feupnea_SEM = VcritvalueSEMA/Veupnea.mean;
Vactive.Active_Pcrit = PcritvalueA;
Vactive.Active_Pcrit_SEM = PcritvalueSEMA;
Vactive.Direct_Feupnea_mean = mean(Vactive.Vactive_list(Vactive.CPAP_list>-1&Vactive.CPAP_list<1))/Veupnea.mean;
Vactive.Direct_Feupnea_median = median(Vactive.Vactive_list(Vactive.CPAP_list>-1&Vactive.CPAP_list<1))/Veupnea.mean;
Vactive.Direct_Feupnea_N = length(Vactive.Vactive_list(Vactive.CPAP_list>-1&Vactive.CPAP_list<1));
Vactive.Direct_Feupnea_SEM = std(Vactive.Vactive_list(Vactive.CPAP_list>-1&Vactive.CPAP_list<1))/(Vactive.Direct_Feupnea_N^0.5)/Veupnea.mean;



%% Save Vactive/Vpassive
if saveontherun
    save(filenameanddir,'Pcrit','Vactive','Vpassive','-append');
    save(filenameanddir,'Vactive','-append');
end

%% Find loop gains

% Select first 3 breaths of dial-up for loop gain

global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Flow','Pes','Edi'};
PlotAndSelectData(signallist,Time,clearvalues);


%% --save
LG.ranges=xvalues;
if saveontherun
    save(filenameanddir,'LG','-append');
end
%% Analyse loop gains

for i=1:size(LG.ranges,1)
    figure(3)
    hold('off');
    LG.data(i,1:5)=NaN; %LGdata contains disturbance, response, veupnea_local, mean(pmask disturbance), mean(pmask response)
    %     if tryLG==1
    %         continue
    %     end
    leftX=30;
    rightX=120;
    deltaX=5;
    lowerX=LG.ranges(i,1)-leftX;
    upperX=LG.ranges(i,2)+rightX;
    
    I=find(Time>lowerX&Time<upperX);
    I2=find((Time>LG.ranges(i,1)-deltaX)&(Time<LG.ranges(i,2)+deltaX));
    I3=find((Time>lowerX)&(Time<LG.ranges(i,1)));
    
    ax3(1)=subplot(4,1,1); plot(Time(I),Flow(I),'k',[LG.ranges(i,1) LG.ranges(i,1)],[-1 1],'r',[LG.ranges(i,2) LG.ranges(i,2)],[-1 1],'r');
    
    vol1=cumsum(Flow(I)-mean(Flow(I)))*dt;
    
    vol2=cumsum(Flow(I2)-mean(Flow(I2)))*dt;
    vol3=cumsum(Flow(I3)-mean(Flow(I3)))*dt;
    
    thresh2=std(vol2)/5;
    
    [min_list, max_list] = peakdet(-vol2,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
    
    ax3(3)=subplot(4,1,3); plot(Time(I),Pmask(I),'k');
    %ax3(4)=subplot(4,1,4); plot(Time(I),Pepi(I),'k');
    ax3(2)=subplot(4,1,2); plot(Time(I2),vol2,'k',[LG.ranges(i,1) LG.ranges(i,1)],[0 1],'r',[LG.ranges(i,2) LG.ranges(i,2)],[0 1],'r');
    hold('on')
    ax3(2)=subplot(4,1,2); plot(Time(I2(1)-1+max_list(:,1)),vol2(max_list(:,1)),'k.');
    ax3(2)=subplot(4,1,2); plot(Time(I2(1)-1+min_list(:,1)),vol2(min_list(:,1)),'r.');
    %[ymin,temptemp]=get(ax2(4),'ylim')
    
    vol3=cumsum(Flow(I3)-mean(Flow(I3)))*dt;
    
    thresh2=std(vol3)/5;
    [min_list2, max_list2] = peakdet(-vol3,thresh2); %in this configuration, max_list(i,1) > minlist(i,1).
    
    ax3(2)=subplot(4,1,2); plot(Time(I3),vol3,'k',[LG.ranges(i,1) LG.ranges(i,1)],[0 1],'r',[LG.ranges(i,2) LG.ranges(i,2)],[0 1],'r');
    ax3(2)=subplot(4,1,2); plot(Time(I3(1)-1+max_list2(:,1)),vol3(max_list2(:,1)),'ko');
    ax3(2)=subplot(4,1,2); plot(Time(I3(1)-1+min_list2(:,1)),vol3(min_list2(:,1)),'ro');
    linkaxes(ax3,'x');
    set(gca,'xlim',[Time(I(1)) Time(I(end))]);
    %Use data?
    [~,~,button] = ginput(1);
    if button==3
        continue
    end
    
    %%RESPONSE%%
    lefti=find(Time(I2(1)-1+min_list(:,1))<LG.ranges(i,1),1,'last'); %find range for response analysis
    righti=find(Time(I2(1)-1+min_list(:,1))>LG.ranges(i,2),1,'first')-1; %find range for response analysis
    if isempty(righti)
        righti=length(min_list)-1;
    end
    
    ax3(2)=subplot(4,1,2); plot([Time(I2(1)-1+min_list(lefti,1)) Time(I2(1)-1+min_list(lefti,1))],[-1 1],'r:',[Time(I2(1)-1+min_list(righti+1,1)) Time(I2(1)-1+min_list(righti+1,1))],[-1 1],'r:');
    
    %calculate drift over the full response range
    xpoly = Time(I2(1)-1+min_list(lefti:(righti+1),1));
    ypoly = vol2(min_list(lefti:(righti+1),1));
    ppoly = polyfit(xpoly,ypoly,1);
    Yvoldrift = polyval(ppoly,Time(I));
    ax3(2)=subplot(4,1,2); plot(xpoly,polyval(ppoly,xpoly));
    
    %calculate tidal vol of breath 1 of response
    Yvol1_exp = ypoly(1); %end exp
    Yvol1_insp = vol2(max_list(lefti,1)); %end insp
    Yvoldrift_exp = polyval(ppoly,xpoly(1));
    Yvoldrift_insp = polyval(ppoly,xpoly(2));
    ttot = (xpoly(2)-xpoly(1));
    endexpvol=Yvol1_exp-Yvoldrift_exp;
    endinspvol=Yvol1_insp-Yvoldrift_insp;
    volbreath=(Yvol1_insp-Yvoldrift_insp)-(Yvol1_exp-Yvoldrift_exp); %endinsp vol minus start insp vol
    VEbreathR=volbreath/ttot*60;
    
    CPAPR=mean(Pmask(I2(1)-1+min_list(lefti:(righti+1),1)));
    
    %%DISTURBANCE%%
    rightDi=find(Time(I3(1)-1+min_list2(:,1))<LG.ranges(i,1),1,'last')-1; %find range for response analysis
    if rightDi>=5 %assume zero flow for now, will need to work on this.
        
        leftDi=rightDi-4; %find range for response analysis
        
        ax3(2)=subplot(4,1,2); plot([Time(I3(1)-1+min_list2(leftDi,1)) Time(I3(1)-1+min_list2(leftDi,1))],[-1 1],'g:',[Time(I3(1)-1+min_list2(rightDi+1,1)) Time(I3(1)-1+min_list2(rightDi+1,1))],[-1 1],'g:');
        
        xpoly = Time(I3(1)-1+min_list2(leftDi:(rightDi),1));
        ypoly = vol3(min_list2(leftDi:(rightDi),1));
        ppoly = polyfit(xpoly,ypoly,1);
        Yvoldrift = polyval(ppoly,Time(I));
        ax3(2)=subplot(4,1,2); plot(xpoly,polyval(ppoly,xpoly));
        
        Yvol1_exp = ypoly(1:end-1); %end exp
        Yvol1_insp = vol3(max_list2(leftDi:(rightDi-1),1)); %end insp
        Yvoldrift_exp = polyval(ppoly,xpoly); Yvoldrift_exp(end)=[];
        Yvoldrift_insp = polyval(ppoly,Time(I3(1)-1+max_list2(leftDi:(rightDi-1),1)));
        ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
        endexpvol=Yvol1_exp-Yvoldrift_exp;
        endinspvol=Yvol1_insp-Yvoldrift_insp;
        volbreath=(Yvol1_insp-Yvoldrift_insp)-(Yvol1_exp-Yvoldrift_exp); %endinsp vol minus start insp vol
        VEbreathD=volbreath/ttot*60;
        
        CPAPD=mean(Pmask(I3(1)-1+min_list2(leftDi:(rightDi+1),1)));
        
        
    else
        VEbreathD=0;
        CPAPD=Pmask(I3(1)-5*Fs);
    end
    %%LOCAL EUPNEA%%
    %select eupnea first left then right
    if 1
        
        %%%%up to here%%%%
        
        disp('select Veupnea local');
        [LGeupneaX(i,1),~,button] = ginput(1);
        if button==3
            VEeupnea_local=NaN;
        else
            [LGeupneaX(i,2),~,button] = ginput(1);
            
            I4=find((Time>LGeupneaX(i,1)-deltaX)&(Time<LGeupneaX(i,2)+deltaX));
            vol4=cumsum(Flow(I4)-mean(Flow(I4)))*dt;
            thresh4=std(vol4)/5;
            [min_list4, max_list4] = peakdet(-vol4,thresh4); %in this configuration, max_list(i,1) > minlist(i,1).
            
            %ax3(1)=subplot(4,1,1); plot(Time(I4),Flow(I4),'k',[LG.ranges(i,1) LG.ranges(i,1)],[-1 1],'r',[LG.ranges(i,2) LG.ranges(i,2)],[-1 1],'r');
            %set(ax3(1),'xlim',[Time(I2(1)) Time(I4(end))])
            
            ax3(2)=subplot(4,1,2); plot(Time(I4),vol4,'k',[LGeupneaX(i,1) LGeupneaX(i,1)],[0 1],'r',[LGeupneaX(i,2) LGeupneaX(i,2)],[0 1],'r');
            ax3(2)=subplot(4,1,2); plot(Time(I4(1)-1+max_list4(:,1)),vol4(max_list4(:,1)),'ko');
            ax3(2)=subplot(4,1,2); plot(Time(I4(1)-1+min_list4(:,1)),vol4(min_list4(:,1)),'ro');
            linkaxes(ax3,'x');
            set(gca,'xlim',[Time(I(1)) Time(I(end))]);
            
            leftEi=find(Time(I4(1)-1+min_list4(:,1))<LGeupneaX(i,1),1,'last'); %find range for response analysis
            if isempty(leftEi)
                leftEi=1;
            end
            
            rightEi=find(Time(I4(1)-1+min_list4(:,1))>LGeupneaX(i,2),1,'first')-1; %find range for response analysis
            
            if isempty(rightEi) %assuming we are at the right hand limit.
                rightEi=size(min_list4,1)-1;
            end
            ax3(2)=subplot(4,1,2); plot([Time(I4(1)-1+min_list4(leftEi,1)) Time(I4(1)-1+min_list4(leftEi,1))],[-1 1],'b:',[Time(I4(1)-1+min_list4(rightEi+1,1)) Time(I4(1)-1+min_list4(rightEi+1,1))],[-1 1],'b:');
            
            xpoly = Time(I4(1)-1+min_list4(leftEi:(rightEi+1),1));
            ypoly = vol4(min_list4(leftEi:(rightEi+1),1));
            ppoly = polyfit(xpoly,ypoly,1);
            Yvoldrift = polyval(ppoly,Time(I4));
            ax3(2)=subplot(4,1,2); plot(xpoly,polyval(ppoly,xpoly));
            linkaxes(ax3,'x');
            set(gca,'xlim',[Time(I(1)) Time(I(end))]);
            Yvol1_exp = ypoly(1:end-1); %end exp
            Yvol1_insp = vol4(max_list4(leftEi:(rightEi),1)); %end insp
            Yvoldrift_exp = polyval(ppoly,xpoly); Yvoldrift_exp(end)=[];
            Yvoldrift_insp = polyval(ppoly,Time(I4(1)-1+max_list4(leftEi:(rightEi),1)));
            ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
            endexpvol=Yvol1_exp-Yvoldrift_exp;
            endinspvol=Yvol1_insp-Yvoldrift_insp;
            volbreath=(Yvol1_insp-Yvoldrift_insp)-(Yvol1_exp-Yvoldrift_exp); %endinsp vol minus start insp vol
            VEbreathEupnea=volbreath/ttot*60;
            VEeupnea_local=mean(VEbreathEupnea)
        end
    else
        VEeupnea_local=NaN;
    end
    
    
    LG.data(i,:)=[mean(VEbreathD) mean(VEbreathR) VEeupnea_local CPAPD CPAPR];
    hold('off');
    %     pause
end

%% --save
if saveontherun
    save(filenameanddir,'LG','-append');
end
%% LG summary data
mindisturbance=Veupnea.mean/20;
LG.data2 = LG.data(sum(~isnan(LG.data(:,:))')'~=0,:);
used_local_Veup = isnan(LG.data2(:,3));

ignorelocaleupnea=0
if ~ignorelocaleupnea
    LG.data2(isnan(LG.data2(:,3)),3)=Veupnea.mean;
else
    LG.data2(:,3)=Veupnea.mean;
end
%Disturbance in col 6
LG.data2(:,6)=LG.data2(:,3)-LG.data2(:,1);

%remove line if Vdist is <0.5L/min
LG.data2((LG.data2(:,6)<mindisturbance),:)=[];

%Response in col 7
LG.data2(:,7)=LG.data2(:,2)-LG.data2(:,3);

%remove line if Vresponse is <0
LG.data2((LG.data2(:,7)<0),:)=[];

LG.median = median(LG.data2(:,7)./LG.data2(:,6));
LG.mean = mean(LG.data2(:,7)./LG.data2(:,6));
LG.SEM = std(LG.data2(:,7)./LG.data2(:,6))/sqrt(size(LG.data2,1));
LG.N = size(LG.data2,1);

ArThres.Direct_Feupnea = median(LG.data2(:,2)./LG.data2(:,3));
ArThres.Indirect_Feupnea = 1+(1-Varousal.Feupnea)*(LG.median);

ArThres.LGmeasuredatArThres=0;

if ArThres.LGmeasuredatArThres
    ArThres.Feupnea = ArThres.Direct_Feupnea;
else
    ArThres.Feupnea = ArThres.Indirect_Feupnea;
end
%% --save
if saveontherun&&~saveattheend
    save(filenameanddir,'LG','ArThres','-append');
end
%% save all
if saveattheend
    save(filenameanddir,'LG','ArThres','ArThresPes','Pcrit','Vactive','Varousal','Veupnea','-append');
    save(filenameanddir,'ArThresPes','Pcrit','Vactive','Varousal','Veupnea','CCWCPAP','SpO2wake','-append');
end


%% Find Edi -> Vflow
% for calculating respiratory impedance and Vdrive
% Select left and right of ranges
clear EdiToFlow
global xvalues yvalues range ax2
clearvalues=1;
range=5*60;
signallist={'[EpochsXHz EventsArXHz WakeSleep]','Pmask','Flow','Pes','Edi'};
PlotAndSelectData(signallist,Time,clearvalues);

%%
EdiToFlow.ranges=xvalues;
if 0
    EdiToFlow.ranges = [rangedata;EdiToFlow.ranges];
end
%%
if saveontherun
    save(filenameanddir,'EdiToFlow','-append');
end

%% Edi to Flow: Find best power
EdiToFlow.Powers = 1.2:-0.1:0.5;
[EdiToFlow.meanRsquared] = FindPowerRespSysMech(EdiToFlow,EdiToFlow.Powers,Edi,Flow,Time,dt);
[EdiToFlow.meanRsquaredmax,tempi]=max(EdiToFlow.meanRsquared);
EdiToFlow.bestPower = EdiToFlow.Powers(tempi);
figure(345)
plot(EdiToFlow.Powers,EdiToFlow.meanRsquared,'.'); ylabel('Rsq'); xlabel('Edi exponent');
hold('on');
plot(EdiToFlow.bestPower,EdiToFlow.meanRsquaredmax,'ro');
hold('off');
if tempi>1&tempi<length(EdiToFlow.Powers)
    range = (tempi-1:tempi+1);
    [x_peak,y_peak] = PeakFitQuadratic(EdiToFlow.Powers(range),EdiToFlow.meanRsquared(range));
    EdiToFlow.bestPower = x_peak;
    hold('on');
    plot(x_peak,y_peak,'r.');
    hold('off');
end

%% Estimate respiratory impedance (Edi)

% if ~isfield('EdiToFlow','power_maxcorrPmus')
%     EdiToFlow.usepower = 0;
% end

parameters_x=[];
for x=1:size(EdiToFlow.ranges,1)
    
    %get data
    I=find(Time>EdiToFlow.ranges(x,1)&Time<EdiToFlow.ranges(x,2));
    
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    tempEdi = Edi(I);
    tempEdi(tempEdi<0) = 0;
    tempEdi = tempEdi.^EdiToFlow.bestPower;
    
    xdata_filtered = filtfilt(B_butter0,A_butter0,tempEdi);
    
    downsamplefactor=1;
    xdata.data = downsample(xdata_filtered,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    ydata = downsample(Flow(I),downsamplefactor);
    
    lsqoptions=optimset('display','off');
    parameters=[6 10 0 0 0]; %guess: R E VL(1)-FRC alpha Vleak delX pwr %alpha = Younes exponent for length-tension, 0 is no effect
    lower=[0.1 0.1 -1  -0.1 -60]; %R E VL(1)-FRC alpha Vleak delX pwr
    upper=[60 60 1 0.5 +60]; %R C VL(1)-FRC alpha Vleak delX  pwr
    
    [modelYguess] = EdiToVflowModel(parameters,xdata);
    
    for i=1:5
        [parameters,resnormC1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@EdiToVflowModel,parameters,xdata,ydata,lower,upper,lsqoptions);
    end
    
    [modelY] = EdiToVflowModel(parameters,xdata);
    
    rsquared_x(x) = sum((modelY-ydata).^2)/sum((ydata-mean(ydata)).^2);
    
    figure(222)
    ax(1)=subplot(3,2,1);
    plot(xdata.time,ydata,'k');
    hold('on');
    plot(xdata.time,modelY,'r');
    hold('off');
    ax(2)=subplot(3,2,3);
    plot(xdata.time,cumsum(ydata-parameters(4))*dt,'k');
    hold('on');
    plot(xdata.time,cumsum(modelY-parameters(4))*dt,'r');
    hold('off');
    ax(3)=subplot(3,2,5);
    plot(xdata.time,xdata.data,'k');
    linkaxes(ax,'x');
    
    subplot(2,2,2);
    plot(ydata,modelY,'k');
    
    subplot(2,2,4);
    plot(cumsum(ydata-parameters(4))*dt,cumsum(modelY-parameters(4))*dt,'k');
    
    parameters_x(x,:)=real(parameters);
    
end

tauRC = parameters_x(:,1)./parameters_x(:,2);
tauRC_median = median(tauRC);
ZRC = abs(parameters_x(:,1)+1i*parameters_x(:,2));

parameters_median = median(parameters_x);
parameters_x = [parameters_x;parameters_median];
close all

%% Test: Use parameters from X on segment Y
EdiToFlow.Parameters = parameters_x;
EdiToFlow.Rsquared=[];

x=size(parameters_x,1)
for i=1:size(parameters_x,1)-1
    y=i;
    
    I=find(Time>EdiToFlow.ranges(y,1)&Time<EdiToFlow.ranges(y,2));
    
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    
    tempEdi = Edi(I);
    tempEdi(tempEdi<0) = 0;
    tempEdi = tempEdi.^EdiToFlow.bestPower;
    
    xdata_filtered = filtfilt(B_butter0,A_butter0,tempEdi);
    
    downsamplefactor=1;
    xdata.data = downsample(xdata_filtered,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    
    ydata = downsample(Flow(I),downsamplefactor);
    
    modelY = EdiToVflowModel(parameters_x(x,:),xdata);
    
    SSE = sum(((modelY-mean(modelY))-(ydata-mean(ydata))).^2);
    SStot = sum((modelY-mean(modelY)).^2);
    
    
    figure()
    ax(1)=subplot(3,2,1)
    plot(xdata.time,ydata-mean(ydata),'k');
    hold('on');
    plot(xdata.time,modelY-mean(modelY),'r');
    hold('off');
    ax(2)=subplot(3,2,3)
    plot(xdata.time,cumsum(ydata-mean(ydata))*dt,'k');
    hold('on')
    plot(xdata.time,cumsum(modelY-mean(modelY))*dt,'r');
    hold('off');
    ax(3)=subplot(3,2,5)
    plot(xdata.time,xdata.data,'k')
    linkaxes(ax,'x')
    
    subplot(2,2,2)
    plot(ydata,modelY,'k');
    
    subplot(2,2,4);
    plot(cumsum(ydata-parameters(4))*dt,cumsum(modelY-parameters(4))*dt,'k')
    
    EdiToFlow.Rsquared(i) = 1-SSE/SStot;
end

close all
%% Calculate Vdrive (FlowEdi) for whole night

Edi(isnan(Edi))=0;

EdiToFlow.useparameterset=size(EdiToFlow.Parameters,1);
parametersOSA=EdiToFlow.Parameters(EdiToFlow.useparameterset,:)
parametersOSA(5)=0;
EdiToFlow.Frescale=1;
parametersOSA(1:2)=parametersOSA(1:2)*EdiToFlow.Frescale;

filter_HFcutoff_butter0 = 8;
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');

tempEdi = Edi;
tempEdi(tempEdi<0) = 0;
tempEdi = tempEdi.^EdiToFlow.bestPower;


Edi1_filtered = filtfilt(B_butter0,A_butter0,tempEdi);

xdata.data = Edi1_filtered;
xdata.time = Time;
xdata.dt = dt;
clear Edi1_filtered tempEdi
FlowEdi = EdiToVflowModel(parametersOSA,xdata);
clear xdata

channels{length(channels)+1}='FlowEdi';

%% Save channel FlowEdi
if ~isstruct(FlowEdi)
    temp = FlowEdi;
end
clear FlowEdi;
FlowEdi=EdiToFlow;
FlowEdi.values=temp;
FlowEdi.interval=dt;
if saveontherun
    save(filenameanddir,'FlowEdi','EdiToFlow','-append');
end
clear FlowEdi
FlowEdi=temp;
%% SKIP: Estimate respiratory impedance (Edi): Power1

parameters_x=[];
for x=1:size(EdiToFlow.ranges,1)
    
    %get data
    I=find(Time>EdiToFlow.ranges(x,1)&Time<EdiToFlow.ranges(x,2));
    
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    tempEdi = Edi(I);
    xdata_filtered = filtfilt(B_butter0,A_butter0,tempEdi);
    
    downsamplefactor=1;
    xdata.data = downsample(xdata_filtered,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    ydata = downsample(Flow(I),downsamplefactor);
    
    lsqoptions=optimset('display','off');
    parameters=[6 10 0 0 0]; %guess: R E VL(1)-FRC alpha Vleak delX pwr %alpha = Younes exponent for length-tension, 0 is no effect
    lower=[0.5 0.1 -1  -0.1 -60]; %R E VL(1)-FRC alpha Vleak delX pwr
    upper=[60 60 1 0.5 +60]; %R C VL(1)-FRC alpha Vleak delX  pwr
    
    [modelYguess] = EdiToVflowModel(parameters,xdata);
    
    for i=1:5
        [parameters,resnormC1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@EdiToVflowModel,parameters,xdata,ydata,lower,upper,lsqoptions);
    end
    
    [modelY] = EdiToVflowModel(parameters,xdata);
    
    rsquared_x(x) = sum((modelY-ydata).^2)/sum((ydata-mean(ydata)).^2);
    
    figure(222)
    ax(1)=subplot(3,2,1);
    plot(xdata.time,ydata,'k');
    hold('on');
    plot(xdata.time,modelY,'r');
    hold('off');
    ax(2)=subplot(3,2,3);
    plot(xdata.time,cumsum(ydata-parameters(4))*dt,'k');
    hold('on');
    plot(xdata.time,cumsum(modelY-parameters(4))*dt,'r');
    hold('off');
    ax(3)=subplot(3,2,5);
    plot(xdata.time,xdata.data,'k');
    linkaxes(ax,'x');
    
    subplot(2,2,2);
    plot(ydata,modelY,'k');
    
    subplot(2,2,4);
    plot(cumsum(ydata-parameters(5))*dt,cumsum(modelY-parameters(4))*dt,'k');
    
    parameters_x(x,:)=real(parameters);
    
end

tauRC = parameters_x(:,1)./parameters_x(:,2);
tauRC_median = median(tauRC);
ZRC = abs(parameters_x(:,1)+1i*parameters_x(:,2));

parameters_median = median(parameters_x);
parameters_x = [parameters_x;parameters_median];
close all

%% SKIP: Test: Use parameters from X on segment Y
EdiToFlow.ParametersPower1 = parameters_x;
EdiToFlow.RsquaredPower1=[];

x=size(parameters_x,1)
for i=1:size(parameters_x,1)-1
    y=i;
    
    I=find(Time>EdiToFlow.ranges(y,1)&Time<EdiToFlow.ranges(y,2));
    
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    
    
    tempEdi = Edi(I);
    
    xdata_filtered = filtfilt(B_butter0,A_butter0,tempEdi);
    
    downsamplefactor=1;
    xdata.data = downsample(xdata_filtered,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    
    ydata = downsample(Flow(I),downsamplefactor);
    
    modelY = EdiToVflowModel(parameters_x(x,:),xdata);
    
    SSE = sum(((modelY-mean(modelY))-(ydata-mean(ydata))).^2);
    SStot = sum((modelY-mean(modelY)).^2);
    
    
    figure()
    ax(1)=subplot(3,2,1);
    plot(xdata.time,ydata-mean(ydata),'k');
    hold('on');
    plot(xdata.time,modelY-mean(modelY),'r');
    hold('off');
    ax(2)=subplot(3,2,3);
    plot(xdata.time,cumsum(ydata-mean(ydata))*dt,'k');
    hold('on');
    plot(xdata.time,cumsum(modelY-mean(modelY))*dt,'r');
    hold('off');
    ax(3)=subplot(3,2,5);
    plot(xdata.time,xdata.data,'k');
    linkaxes(ax,'x');
    
    subplot(2,2,2);
    plot(ydata,modelY,'k');
    
    subplot(2,2,4);
    plot(cumsum(ydata-parameters(4))*dt,cumsum(modelY-parameters(4))*dt,'k');
    
    EdiToFlow.RsquaredPower1(i) = 1-SSE/SStot;
end

close all
%% SKIP: Calculate Vdrive (FlowEdi) for whole night

EdiToFlow.useparameterset=size(parameters_x,1);
parametersOSA=parameters_x(EdiToFlow.useparameterset,:)
parametersOSA(5)=0;

parametersOSA(1:2)=parametersOSA(1:2)*EdiToFlow.Frescale;

filter_HFcutoff_butter0 = 8;
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');

tempEdi = Edi;
tempEdi(tempEdi<0) = 0;
Edi1_filtered = filtfilt(B_butter0,A_butter0,tempEdi);

xdata.data = Edi1_filtered;
xdata.time = Time;
xdata.dt = dt;
clear Edi1_filtered tempEdi
FlowEdiPower1 = EdiToVflowModel(parametersOSA,xdata);
clear xdata

channels{length(channels)+1}='FlowEdiPower1';

% %% save
% save(filenameanddir,'EdiToFlow','-append');
%  close all

%% SKIP: Estimate respiratory impedance (EdiSqrt)

EdiToFlow.usepower = 1;
% if ~isfield('EdiToFlow','power_maxcorrPmus')
%     EdiToFlow.usepower = 0;
% end

parameters_x=[];
for x=1:size(EdiToFlow.ranges,1)
    
    %get data
    I=find(Time>EdiToFlow.ranges(x,1)&Time<EdiToFlow.ranges(x,2));
    
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    tempEdi = Edi(I);
    tempEdi(tempEdi<0) = 0;
    if EdiToFlow.usepower
        if 1
            tempEdi = tempEdi.^0.5;
        else
            tempEdi = tempEdi.^EdiToFlow.power_maxcorrPmus;
        end
        
    end
    
    xdata_filtered = filtfilt(B_butter0,A_butter0,tempEdi);
    
    
    
    downsamplefactor=1;
    xdata.data = downsample(xdata_filtered,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    ydata = downsample(Flow(I),downsamplefactor);
    
    lsqoptions=optimset('display','off');
    parameters=[6 10 0 0 0]; %guess: R E VL(1)-FRC alpha Vleak delX pwr %alpha = Younes exponent for length-tension, 0 is no effect
    lower=[0.5 0.1 -1  -0.1 -60]; %R E VL(1)-FRC alpha Vleak delX pwr
    upper=[60 60 1 0.5 +60]; %R C VL(1)-FRC alpha Vleak delX  pwr
    
    [modelYguess] = EdiToVflowModel(parameters,xdata);
    
    for i=1:5
        [parameters,resnormC1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@EdiToVflowModel,parameters,xdata,ydata,lower,upper,lsqoptions);
    end
    
    [modelY] = EdiToVflowModel(parameters,xdata);
    
    rsquared_x(x) = sum((modelY-ydata).^2)/sum((ydata-mean(ydata)).^2);
    
    figure()
    ax(1)=subplot(3,2,1)
    plot(xdata.time,ydata,'k')
    hold('on')
    plot(xdata.time,modelY,'r')
    hold('off');
    ax(2)=subplot(3,2,3)
    plot(xdata.time,cumsum(ydata-parameters(4))*dt,'k')
    hold('on')
    plot(xdata.time,cumsum(modelY-parameters(4))*dt,'r')
    hold('off');
    ax(3)=subplot(3,2,5)
    plot(xdata.time,xdata.data,'k')
    linkaxes(ax,'x')
    
    subplot(2,2,2)
    plot(ydata,modelY,'k');
    
    subplot(2,2,4);
    plot(cumsum(ydata-parameters(5))*dt,cumsum(modelY-parameters(4))*dt,'k')
    
    parameters_x(x,:)=real(parameters);
    
end

tauRC = parameters_x(:,1)./parameters_x(:,2);
tauRC_median = median(tauRC);
ZRC = abs(parameters_x(:,1)+1i*parameters_x(:,2));

parameters_median = median(parameters_x);
parameters_x = [parameters_x;parameters_median];
close all

%% SKIP: Test: Use parameters from X on segment Y
EdiToFlow.ParametersSqrt = parameters_x;
EdiToFlow.RsquaredSqrt=[];

x=size(parameters_x,1)
for i=1:size(parameters_x,1)-1
    y=i;
    
    I=find(Time>EdiToFlow.ranges(y,1)&Time<EdiToFlow.ranges(y,2));
    
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    
    
    tempEdi = Edi(I);
    tempEdi(tempEdi<0) = 0;
    if EdiToFlow.usepower
        if 1
            tempEdi = tempEdi.^0.5;
        else
            tempEdi = tempEdi.^EdiToFlow.power_maxcorrPmus;
        end
    end
    xdata_filtered = filtfilt(B_butter0,A_butter0,tempEdi);
    
    downsamplefactor=1;
    xdata.data = downsample(xdata_filtered,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    
    ydata = downsample(Flow(I),downsamplefactor);
    
    modelY = EdiToVflowModel(parameters_x(x,:),xdata);
    
    SSE = sum(((modelY-mean(modelY))-(ydata-mean(ydata))).^2);
    SStot = sum((modelY-mean(modelY)).^2);
    
    
    figure()
    ax(1)=subplot(3,2,1)
    plot(xdata.time,ydata-mean(ydata),'k');
    hold('on');
    plot(xdata.time,modelY-mean(modelY),'r');
    hold('off');
    ax(2)=subplot(3,2,3)
    plot(xdata.time,cumsum(ydata-mean(ydata))*dt,'k');
    hold('on')
    plot(xdata.time,cumsum(modelY-mean(modelY))*dt,'r');
    hold('off');
    ax(3)=subplot(3,2,5)
    plot(xdata.time,xdata.data,'k')
    linkaxes(ax,'x')
    
    subplot(2,2,2)
    plot(ydata,modelY,'k');
    
    subplot(2,2,4);
    plot(cumsum(ydata-parameters(4))*dt,cumsum(modelY-parameters(4))*dt,'k')
    
    EdiToFlow.RsquaredSqrt(i) = 1-SSE/SStot;
end

close all
%% SKIP: Calculate Vdrive (FlowEdiSqrt) for whole night

EdiToFlow.useparametersetSqrt=size(parameters_x,1);
parametersOSA=parameters_x(EdiToFlow.useparametersetSqrt,:)
parametersOSA(5)=0;

parametersOSA(1:2)=parametersOSA(1:2)*EdiToFlow.Frescale;

filter_HFcutoff_butter0 = 8;
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');

tempEdi = Edi;
tempEdi(tempEdi<0) = 0;

tempEdi = tempEdi.^0.5;


Edi1_filtered = filtfilt(B_butter0,A_butter0,tempEdi);

xdata.data = Edi1_filtered;
xdata.time = Time;
xdata.dt = dt;
clear Edi1_filtered tempEdi
FlowEdiSqrt = EdiToVflowModel(parametersOSA,xdata);
clear xdata

channels{length(channels)+1}='FlowEdiSqrt';
% %%
% save(filenameanddir,'EdiToFlow','-append');
%  close all

%% SKIP: Estimate respiratory impedance (EdiPower75)

% if ~isfield('EdiToFlow','power_maxcorrPmus')
%     EdiToFlow.usepower = 0;
% end

parameters_x=[];
for x=1:size(EdiToFlow.ranges,1)
    
    %get data
    I=find(Time>EdiToFlow.ranges(x,1)&Time<EdiToFlow.ranges(x,2));
    
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    tempEdi = Edi(I);
    tempEdi(tempEdi<0) = 0;
    xdata_filtered = filtfilt(B_butter0,A_butter0,tempEdi.^0.75);
    
    downsamplefactor=1;
    xdata.data = downsample(xdata_filtered,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    ydata = downsample(Flow(I),downsamplefactor);
    
    lsqoptions=optimset('display','off');
    parameters=[6 10 0 0 0]; %guess: R E VL(1)-FRC alpha Vleak delX pwr %alpha = Younes exponent for length-tension, 0 is no effect
    lower=[0.5 0.1 -1  -0.1 -60]; %R E VL(1)-FRC alpha Vleak delX pwr
    upper=[60 60 1 0.5 +60]; %R C VL(1)-FRC alpha Vleak delX  pwr
    
    [modelYguess] = EdiToVflowModel(parameters,xdata);
    
    for i=1:5
        [parameters,resnormC1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@EdiToVflowModel,parameters,xdata,ydata,lower,upper,lsqoptions);
    end
    
    [modelY] = EdiToVflowModel(parameters,xdata);
    
    rsquared_x(x) = sum((modelY-ydata).^2)/sum((ydata-mean(ydata)).^2);
    
    figure()
    ax(1)=subplot(3,2,1)
    plot(xdata.time,ydata,'k')
    hold('on')
    plot(xdata.time,modelY,'r')
    hold('off');
    ax(2)=subplot(3,2,3)
    plot(xdata.time,cumsum(ydata-parameters(4))*dt,'k')
    hold('on')
    plot(xdata.time,cumsum(modelY-parameters(4))*dt,'r')
    hold('off');
    ax(3)=subplot(3,2,5)
    plot(xdata.time,xdata.data,'k')
    linkaxes(ax,'x')
    
    subplot(2,2,2)
    plot(ydata,modelY,'k');
    
    subplot(2,2,4);
    plot(cumsum(ydata-parameters(5))*dt,cumsum(modelY-parameters(4))*dt,'k')
    
    parameters_x(x,:)=real(parameters);
    
end

tauRC = parameters_x(:,1)./parameters_x(:,2);
tauRC_median = median(tauRC);
ZRC = abs(parameters_x(:,1)+1i*parameters_x(:,2));

parameters_median = median(parameters_x);
parameters_x = [parameters_x;parameters_median];
close all

%% SKIP: Test: Use parameters from X on segment Y
EdiToFlow.ParametersPower75 = parameters_x;
EdiToFlow.RsquaredPower75=[];

x=size(parameters_x,1)
for i=1:size(parameters_x,1)-1
    y=i;
    
    I=find(Time>EdiToFlow.ranges(y,1)&Time<EdiToFlow.ranges(y,2));
    
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    
    tempEdi = Edi(I);
    tempEdi(tempEdi<0) = 0;
    xdata_filtered = filtfilt(B_butter0,A_butter0,tempEdi.^0.75);
    
    downsamplefactor=1;
    xdata.data = downsample(xdata_filtered,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    
    ydata = downsample(Flow(I),downsamplefactor);
    
    modelY = EdiToVflowModel(parameters_x(x,:),xdata);
    
    SSE = sum(((modelY-mean(modelY))-(ydata-mean(ydata))).^2);
    SStot = sum((modelY-mean(modelY)).^2);
    
    
    figure()
    ax(1)=subplot(3,2,1)
    plot(xdata.time,ydata-mean(ydata),'k');
    hold('on');
    plot(xdata.time,modelY-mean(modelY),'r');
    hold('off');
    ax(2)=subplot(3,2,3)
    plot(xdata.time,cumsum(ydata-mean(ydata))*dt,'k');
    hold('on')
    plot(xdata.time,cumsum(modelY-mean(modelY))*dt,'r');
    hold('off');
    ax(3)=subplot(3,2,5)
    plot(xdata.time,xdata.data,'k')
    linkaxes(ax,'x')
    
    subplot(2,2,2)
    plot(ydata,modelY,'k');
    
    subplot(2,2,4);
    plot(cumsum(ydata-parameters(4))*dt,cumsum(modelY-parameters(4))*dt,'k')
    
    EdiToFlow.RsquaredPower75(i) = 1-SSE/SStot;
end

close all

%% SKIP: Calculate Vdrive (FlowEdiPower75) for whole night

EdiToFlow.useparametersetPower75=size(parameters_x,1);
parametersOSA=parameters_x(EdiToFlow.useparametersetPower75,:)
parametersOSA(5)=0;

parametersOSA(1:2)=parametersOSA(1:2)*EdiToFlow.Frescale;

filter_HFcutoff_butter0 = 8;
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');

tempEdi = Edi;
tempEdi(tempEdi<0) = 0;

tempEdi = tempEdi.^0.75;

Edi1_filtered = filtfilt(B_butter0,A_butter0,tempEdi);

xdata.data = Edi1_filtered;
xdata.time = Time;
xdata.dt = dt;
clear Edi1_filtered tempEdi
FlowEdiPower75 = EdiToVflowModel(parametersOSA,xdata);
clear xdata

channels{length(channels)+1}='FlowEdiPower75';

%% SKIP: Compared Edi models

EdiToFlow.compareEdimodels_mean = [mean(EdiToFlow.Rsquared) mean(EdiToFlow.RsquaredPower1) mean(EdiToFlow.RsquaredSqrt) mean(EdiToFlow.RsquaredPower75)];
EdiToFlow.compareEdimodels_median = [median(EdiToFlow.Rsquared) median(EdiToFlow.RsquaredPower1) median(EdiToFlow.RsquaredSqrt) median(EdiToFlow.RsquaredPower75)];

%% SKIP
if 0
    filenameanddir(1)='F';
end
save(filenameanddir,'EdiToFlow','-append');
close all
%% SKIP: ObstructedXHz and CentralapneasXHz for whole night
ObstructedXHz = 0*Time;
CentralapneasXHz = 0*Time;
channels{length(channels)+1}='ObstructedXHz';
channels{length(channels)+1}='CentralapneasXHz';
for i=1:length(Evts.starttimes)
    if sum(Evts.codes(i)==[2 4 5 7 9])==0
        continue
    end
    lefti=round((Evts.starttimes(i)-Time(1))/dt+1);
    righti=round((Evts.endtimes(i)-Time(1))/dt+1);
    if lefti<i&&righti>=1
        lefti=1;
    end
    if lefti<=length(Time)&&righti>length(Time)
        righti=length(Time);
    end
    ObstructedXHz(lefti:righti)=1;
end
for i=1:length(Evts.starttimes)
    if Evts.codes(i)~=3
        continue
    end
    lefti=round((Evts.starttimes(i)-Time(1))/dt+1);
    righti=round((Evts.endtimes(i)-Time(1))/dt+1);
    if lefti<i&&righti>=1
        lefti=1;
    end
    if lefti<=length(Time)&&righti>length(Time)
        righti=length(Time);
    end
    CentralapneasXHz(lefti:righti)=1;
end

figure(100)
plot(Time,ObstructedXHz)

close all

%% Pes To Flow: Filtered Pes (gentle low pass)
filter_HFcutoff_butter0 = 8;
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');

Pes(isnan(Pes))=nanmean(Pes);
Pes1_filtered = filtfilt(B_butter0,A_butter0,Pes);

%%

%% check Editoflow ranges for good Pes
close all

button=[];
for x=1:size(EdiToFlow.ranges,1)
    figure(200);
    X=3;
    I=find(Time>EdiToFlow.ranges(x,1)&Time<EdiToFlow.ranges(x,2));
    ax200(1)=subplot(X,1,1);plot(Time(I),Flow(I));
    ax200(2)=subplot(X,1,2);plot(Time(I),Edi(I));
    ax200(3)=subplot(X,1,3);plot(Time(I),Pes1_filtered(I));
    [~,~,button(x)]=ginput(1);
end

PesToFlow.ranges = EdiToFlow.ranges(button==1,:);

%% Estimate respiratory impedance (Pes)

parameters_x=[];
for x=1:size(PesToFlow.ranges,1)
    try
        %get data
        I=find(Time>PesToFlow.ranges(x,1)&Time<PesToFlow.ranges(x,2));
        
        xdata_filtered = Pes(I);
        
        minimum_figs=1;
        [I_edi,~,~,~,~]=Vflowanalysis2(FlowEdi(I),Time(I),dt,minimum_figs);
        baselinePes = median(xdata_filtered(I_edi.starti));
        
        downsamplefactor=1;
        xdata.data = downsample(-xdata_filtered+baselinePes,downsamplefactor); %downsample(Pes(I),downsamplefactor);
        xdata.time = downsample(Time(I),downsamplefactor);
        xdata.dt = dt;
        ydata = downsample(Flow(I),downsamplefactor);
        
        set_alpha=0;
        lsqoptions=optimset('display','off');
        parameters=[6 6 0 0 50]; %guess: R E initialVL Vleak baselinePesErr %alpha = Younes exponent for length-tension, 0 is no effect
        
        lower=[0.1 0.1 -0.01 0 0]; %R E initialVL Vleak baselinePesErr
        upper=[60 60 0.01 0.5 100]; %R E initialVL Vleak baselinePesErr
        
        [modelYguess] = PesToVflowModel(parameters,xdata);
        
        for i=1:5
            [parameters,resnormC1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@PesToVflowModel,parameters,xdata,ydata,lower,upper,lsqoptions);
        end
        
        [modelY] = PesToVflowModel(parameters,xdata);
        
        rsquared_x(x) = sum((modelY-ydata).^2)/sum((ydata-mean(ydata)).^2);
        
        figure()
        ax(1)=subplot(3,2,1)
        plot(xdata.time,ydata,'k')
        hold('on')
        plot(xdata.time,modelY,'r')
        hold('off');
        ax(2)=subplot(3,2,3)
        plot(xdata.time,cumsum(ydata-parameters(5))*dt,'k')
        hold('on')
        plot(xdata.time,cumsum(modelY-parameters(5))*dt,'r')
        hold('off');
        ax(3)=subplot(3,2,5)
        plot(xdata.time,xdata.data,'k')
        linkaxes(ax,'x')
        
        subplot(2,2,2)
        plot(ydata,modelY,'k');
        
        subplot(2,2,4);
        plot(cumsum(ydata-parameters(5))*dt,cumsum(modelY-parameters(5))*dt,'k')
        
        parameters_x(x,:)=parameters;
    catch me
        parameters_x(x,:)=NaN*ones(1,5);
    end
end

tauRC = parameters_x(:,1)./parameters_x(:,2);
tauRC_median = nanmedian(tauRC);
ZRC = abs(parameters_x(:,1)+1i*parameters_x(:,2));

parameters_median = nanmedian(parameters_x);
parameters_x = [parameters_x;parameters_median]

%% Test: Use parameters from X on segment Y (Pes)
PesToFlow.Parameters = parameters_x;
PesToFlow.Rsquared=[];

x=size(parameters_x,1);
for i=1:size(parameters_x,1)-1
    try
        y=i;
        
        I=find(Time>PesToFlow.ranges(y,1)&Time<PesToFlow.ranges(y,2));
        
        xdata_filtered = Pes(I);
        
        minimum_figs=1;
        [I_edi,~,~,~,~]=Vflowanalysis2(FlowEdi(I),Time(I),dt,minimum_figs);
        baselinePes = median(xdata_filtered(I_edi.starti));
        
        downsamplefactor=1;
        xdata.data = downsample(-xdata_filtered+baselinePes,downsamplefactor); %downsample(Pes(I),downsamplefactor);
        xdata.time = downsample(Time(I),downsamplefactor);
        xdata.dt = dt;
        
        ydata = downsample(Flow(I),downsamplefactor);
        
        modelY = PesToVflowModel(parameters_x(x,:),xdata); %is this use of 'x' right?
        
        SSE = sum(((modelY-mean(modelY))-(ydata-mean(ydata))).^2);
        SStot = sum((modelY-mean(modelY)).^2);
        
        figure()
        ax(1)=subplot(3,2,1)
        plot(xdata.time,ydata-mean(ydata),'k');
        hold('on');
        plot(xdata.time,modelY-mean(modelY),'r');
        hold('off');
        ax(2)=subplot(3,2,3)
        plot(xdata.time,cumsum(ydata-mean(ydata))*dt,'k');
        hold('on')
        plot(xdata.time,cumsum(modelY-mean(modelY))*dt,'r');
        hold('off');
        ax(3)=subplot(3,2,5)
        plot(xdata.time,xdata.data,'k')
        linkaxes(ax,'x')
        
        subplot(2,2,2)
        plot(ydata,modelY,'k');
        
        subplot(2,2,4);
        plot(cumsum(ydata-parameters(5))*dt,cumsum(modelY-parameters(5))*dt,'k')
        
        PesToFlow.Rsquared(i) = 1-SSE/SStot;
    catch me
        PesToFlow.Rsquared(i) = NaN;
    end
    
end

close all

%% Calculate Vdrive (FlowPes) for whole night

PesToFlow.useparameterset=size(parameters_x,1);
parametersOSA=parameters_x(PesToFlow.useparameterset,:)
parametersOSA(5)=0;
Frescale=1;
parametersOSA(1:2)=parametersOSA(1:2)*Frescale;

xdata.data = -Pes1_filtered; %downsample(Pes(I),downsamplefactor);
xdata.time = Time;
xdata.dt = dt;
parametersOSA(3) = -xdata.data(1)/parametersOSA(2);
FlowPes = PesToVflowModel(parametersOSA,xdata);
clear xdata

channels{length(channels)+1}='FlowPes';

if 0
    figure(101); plot(Time,[Flow FlowPes]);
end

%% Save channel FlowPes
if ~isstruct(FlowPes)
    temp = FlowPes;
end
clear FlowPes;
FlowPes=PesToFlow;
FlowPes.values=temp;
FlowPes.interval=dt;
if saveontherun
    save(filenameanddir,'FlowPes','PesToFlow','-append');
end
clear FlowPes
FlowPes=temp;

%% Estimate respiratory impedance (Pmus = Pes + Ecw.Vol)

filter_HFcutoff_butter0 = 8;
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
Pcw = filtfilt(B_butter0,A_butter0,-1/CCWCPAP.Ccw*cumsum(Flow)*dt);

filter_LFcutoff_butterX = 1/30;
filter_orderX = 2;
[B_butterX,A_butterX] = butter(filter_orderX,filter_LFcutoff_butterX/(1/dt/2),'high');
Pcw = filtfilt(B_butterX,A_butterX,Pcw);

parameters_x=[];
for x=1:size(PesToFlow.ranges,1)
    try
        %get data
        I=find(Time>PesToFlow.ranges(x,1)&Time<PesToFlow.ranges(x,2));
        
        filter_HFcutoff_butter0 = 8;
        filter_order0 = 2;
        [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
        
        xdata_filtered = Pes1_filtered(I) + Pcw(I);
        
        minimum_figs=1;
        [I_edi,~,~,~,~]=Vflowanalysis2(FlowEdi(I),Time(I),dt,minimum_figs);
        baselinePes = median(xdata_filtered(I_edi.starti));
        
        downsamplefactor=1;
        xdata.data = downsample(-xdata_filtered+baselinePes,downsamplefactor); %downsample(Pes(I),downsamplefactor);
        xdata.time = downsample(Time(I),downsamplefactor);
        xdata.dt = dt;
        ydata = downsample(Flow(I),downsamplefactor);
        
        set_alpha=0;
        lsqoptions=optimset('display','off');
        parameters=[6 6 0 0 50]; %guess: R E initialVL Vleak baselinePesErr %alpha = Younes exponent for length-tension, 0 is no effect
        
        lower=[0.1 0.1 -0.01 0 0]; %R E initialVL Vleak baselinePesErr
        upper=[60 60 0.01 0.5 100]; %R E initialVL Vleak baselinePesErr
        
        [modelYguess] = PesToVflowModel(parameters,xdata);
        
        for i=1:5
            [parameters,resnormC1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@PesToVflowModel,parameters,xdata,ydata,lower,upper,lsqoptions);
        end
        
        [modelY] = PesToVflowModel(parameters,xdata);
        
        rsquared_x(x) = sum((modelY-ydata).^2)/sum((ydata-mean(ydata)).^2);
        
        figure()
        ax(1)=subplot(3,2,1)
        plot(xdata.time,ydata,'k')
        hold('on')
        plot(xdata.time,modelY,'r')
        hold('off');
        ax(2)=subplot(3,2,3)
        plot(xdata.time,cumsum(ydata-parameters(5))*dt,'k')
        hold('on')
        plot(xdata.time,cumsum(modelY-parameters(5))*dt,'r')
        hold('off');
        ax(3)=subplot(3,2,5)
        plot(xdata.time,xdata.data,'k')
        linkaxes(ax,'x')
        
        subplot(2,2,2)
        plot(ydata,modelY,'k');
        
        subplot(2,2,4);
        plot(cumsum(ydata-parameters(5))*dt,cumsum(modelY-parameters(5))*dt,'k')
        
        parameters_x(x,:)=parameters;
    catch me
        parameters_x(x,:)=NaN*ones(1,5);
    end
end

tauRC = parameters_x(:,1)./parameters_x(:,2);
tauRC_median = nanmedian(tauRC);
ZRC = abs(parameters_x(:,1)+1i*parameters_x(:,2));

parameters_median = nanmedian(parameters_x);
parameters_x = [parameters_x;parameters_median]

PesToFlow.PmusParameters = parameters_x;
%% Test: Use parameters from X on segment Y (Pmus)

PesToFlow.PmusRsquared=[];

x=size(PesToFlow.PmusParameters,1);
for i=1:size(PesToFlow.PmusParameters,1)-1
    try
        y=i;
        
        I=find(Time>PesToFlow.ranges(y,1)&Time<PesToFlow.ranges(y,2));
        
        xdata_filtered = Pes1_filtered(I) + Pcw(I);
        
        minimum_figs=1;
        [I_edi,~,~,~,~]=Vflowanalysis2(FlowEdi(I),Time(I),dt,minimum_figs);
        baselinePes = median(xdata_filtered(I_edi.starti));
        
        downsamplefactor=1;
        xdata.data = downsample(-xdata_filtered+baselinePes,downsamplefactor); %downsample(Pes(I),downsamplefactor);
        xdata.time = downsample(Time(I),downsamplefactor);
        xdata.dt = dt;
        
        ydata = downsample(Flow(I),downsamplefactor);
        
        tempparameters = parameters_x(x,:);
        %tempparameters(2) = tempparameters(2)+1/CCWCPAP.Ccw;
        modelY = PesToVflowModel(parameters_x(x,:),xdata); %is this use of 'x' right?
        
        SSE = sum(((modelY-mean(modelY))-(ydata-mean(ydata))).^2);
        SStot = sum((modelY-mean(modelY)).^2);
        
        
        figure()
        ax(1)=subplot(3,2,1)
        plot(xdata.time,ydata-mean(ydata),'k');
        hold('on');
        plot(xdata.time,modelY-mean(modelY),'r');
        hold('off');
        ax(2)=subplot(3,2,3)
        plot(xdata.time,cumsum(ydata-mean(ydata))*dt,'k');
        hold('on')
        plot(xdata.time,cumsum(modelY-mean(modelY))*dt,'r');
        hold('off');
        ax(3)=subplot(3,2,5)
        plot(xdata.time,xdata.data,'k')
        linkaxes(ax,'x')
        
        subplot(2,2,2)
        plot(ydata,modelY,'k');
        
        subplot(2,2,4);
        plot(cumsum(ydata-parameters(5))*dt,cumsum(modelY-parameters(5))*dt,'k')
        
        PesToFlow.PmusRsquared(i) = 1-SSE/SStot;
    catch me
        PesToFlow.PmusRsquared(i)=NaN;
    end
end

close all



%% Calculate Vdrive (FlowPmus) for whole night

PesToFlow.useparameterset=size(PesToFlow.Parameters,1);
parametersOSA=PesToFlow.PmusParameters(PesToFlow.useparameterset,:)
parametersOSA(5)=0;
%parametersOSA(2)=parametersOSA(2) + 1/CCWCPAP.Ccw;
Frescale=1;
parametersOSA(1:2)=parametersOSA(1:2)*Frescale;

Pmus1_filtered = Pes1_filtered + Pcw; %%%%%%%%%%%%%%%%%% Temp fudge factor

%Pmus1_filtered = filtfilt(B_butter0,A_butter0,Pes - 1/CCWCPAP.Ccw*cumsum(Flow)*dt);

xdata.data = -Pmus1_filtered ; %downsample(Pes(I),downsamplefactor);
xdata.time = Time;
xdata.dt = dt;
%clear Pmus1_filtered
parametersOSA(3) = -xdata.data(1)/parametersOSA(2);
FlowPmus = PesToVflowModel(parametersOSA,xdata);
clear xdata

channels{length(channels)+1}='FlowPmus';

%% Save channel FlowPmus
if ~isstruct(FlowPmus)
    temp = FlowPmus;
end
clear FlowPmus;
FlowPmus=PesToFlow;
FlowPmus.values=temp;
FlowPmus.interval=dt;
if saveontherun
    save(filenameanddir,'FlowPmus','PesToFlow','-append');
end
clear FlowPmus
FlowPmus=temp;

%% END

%% Consider how to employ Pmus best power: challenge is that there is no known baseline:
%might be easier to transform VTs afterwards and correlate the tidal volumes
% PesToFlow.Powers = 0.9:0.1:2;
% [EdiToFlow.meanRsquared] = FindPowerRespSysMech(PesToFlow,PesToFlow.Powers,-Pmus1_filtered,Flow,Time,dt);
% [EdiToFlow.meanRsquaredmax,tempi]=max(EdiToFlow.meanRsquared);
% EdiToFlow.bestPower = EdiToFlow.Powers(tempi);
% figure(345)
% plot(EdiToFlow.Powers,EdiToFlow.meanRsquared,'.'); ylabel('Rsq'); xlabel('Edi exponent');
% hold('on');
%         plot(EdiToFlow.bestPower,EdiToFlow.meanRsquaredmax,'ro');
%         hold('off');
%     if tempi>1&tempi<length(EdiToFlow.Powers)
%         range = (tempi-1:tempi+1);
%         [x_peak,y_peak] = PeakFitQuadratic(EdiToFlow.Powers(range),EdiToFlow.meanRsquared(range));
%         EdiToFlow.bestPower = x_peak;
%         hold('on');
%         plot(x_peak,y_peak,'r.');
%         hold('off');
%     end

%% Find optimum volume coefficient (i.e. Ecw) using OSA segments [slow]

clear normalizedRsquared normalizedRsquaredX SSEXplusY

NfoldCcw=[0 0.5 1 1.5 2:6];
for j=1:length(NfoldCcw)
    VT_edi_=[];
    VT_pmus_=[];
    VT_pes_ = [];
    VT_pmus2_ = [];
    VT_pmusTest_ = [];
    VT_flow_ = [];
    for y=1:size(OSA.ranges,1)
        try
            I=find(Time>OSA.ranges(y,1)&Time<OSA.ranges(y,2));
            FlowEdiTemp = FlowEdi(I);
            PmusTestTemp = Pes1_filtered(I)+NfoldCcw(j)*Pcw(I);
            
            FlowTemp = Flow(I);
            TimeTemp = Time(I);
            %FlowPmusTemp = FlowPmus(I);
            %FlowPmus2Temp = FlowPmus2(I);
            FlowPesTemp = FlowPes(I);
            minimum_figs=1;
            CPAPrange = [-2 2];
            Pmasktemp = Pmask(I);
            
            downsamplefactor=1;
            xdata.data = downsample(-PmusTestTemp,downsamplefactor); %downsample(Pes(I),downsamplefactor);
            xdata.time = downsample(Time(I),downsamplefactor);
            xdata.dt = dt;
            ratio = PesToFlow.Parameters(end,1)/EdiToFlow.Parameters(end,1);
            modelY = PesToVflowModel(ratio*EdiToFlow.Parameters(end,:),xdata); %is this use of 'x' right?
            
            figure(234);
            ax234(1)=subplot(3,1,1); plot(TimeTemp,modelY);
            ax234(2)=subplot(3,1,2); plot(TimeTemp,FlowEdiTemp);
            ax234(3)=subplot(3,1,3); plot(TimeTemp,FlowTemp);
            linkaxes(ax234,'x');
            
            [I_edi,~,VT_edi,~,~]=Vflowanalysis2(FlowEdiTemp,TimeTemp,dt,minimum_figs);
            [~,~,VT_pmusTest]=Vflowanalysis_knownI2(modelY-mean(modelY),TimeTemp,dt,minimum_figs,I_edi);
            
            CPAP_OSA = median(Pmasktemp(I_edi.starti));
            if CPAP_OSA<CPAPrange(1)||CPAP_OSA>CPAPrange(2)
                continue
            end
            %[~,~,VT_pmus]=Vflowanalysis_knownI2(FlowPmusTemp-mean(FlowPmusTemp),TimeTemp,dt,minimum_figs,I_edi);
            %[~,~,VT_pmus2]=Vflowanalysis_knownI2(FlowPmus2Temp-mean(FlowPmus2Temp),TimeTemp,dt,minimum_figs,I_edi);
            %[~,~,VT_pes]=Vflowanalysis_knownI2(FlowPesTemp-mean(FlowPesTemp),TimeTemp,dt,minimum_figs,I_edi);
            %[~,~,VT_flow]=Vflowanalysis_knownI2(FlowTemp-mean(FlowTemp),TimeTemp,dt,minimum_figs,I_edi);
            
            VT_edi_ = [VT_edi_ VT_edi];
            %VT_pmus_ = [VT_pmus_ VT_pmus];
            %VT_pmus2_ = [VT_pmus2_ VT_pmus2];
            %VT_pes_ = [VT_pes_ VT_pes];
            %VT_flow_ = [VT_flow_ VT_flow];
            VT_pmusTest_ = [VT_pmusTest_ VT_pmusTest];
        catch me
        end
        %figure(100); subplot(2,1,1);plot(VT_edi_,VT_pmus_,'.'); hold('on');
        figure(101); subplot(1,length(NfoldCcw),j); plot(VT_edi_,VT_pmusTest_,'.','markersize',1);
    end
    %rangeEdi = [prctile(VT_edi_,1) prctile(VT_edi_,99)];
    %rangePmus = [prctile(VT_pmusTest_,1) prctile(VT_pmusTest_,99)];
    Nciles = 10; clear VTEdiPrctiles VTPmusTestPrctiles
    for i=1:Nciles
        II = (VT_edi_>=prctile(VT_edi_,(i-1)*(100/Nciles))&VT_edi_<=prctile(VT_edi_,i*(100/Nciles)));
        VTEdiPrctiles(i) = median(VT_edi_(II));
        VTPmusTestPrctiles(i) = median(VT_pmusTest_(II));
    end
    
    
    %criteria = VT_edi_<rangeEdi(1)|VT_pmusTest_<rangePmus(1);
    %VT_edi_(criteria)=[];
    %VT_pmusTest_(criteria)=[];
    xdata = VTEdiPrctiles;
    ydata = VTPmusTestPrctiles;
    lsqoptions=optimset('display','off');
    fun = @(x,xdata)x(1)*(xdata);
    [Xmodel(j)]=lsqcurvefit(fun,1,xdata,ydata,0.1,10,lsqoptions);
    ytemp = ydata/Xmodel(j); %%ytemp is the new ydata, xdata = ymodel because the model is now y=x
    normalizedRsquared(j) = 1-sum((ytemp-xdata).^2)/sum((ytemp-mean(ytemp)).^2);
    normalizedRsquaredX(j) = 1-sum((ytemp-xdata).^2)/sum((xdata-mean(xdata)).^2);
    SSEXplusY(j) = sum((ytemp-xdata).^2); %actual y values (ytemp) - model y values = actual x values - model x values
    figure(100); subplot(1,length(NfoldCcw),j); plot(VT_edi_,VT_pmusTest_/Xmodel(j),'.','markersize',1); hold('on');
    plot(xdata,ydata/Xmodel(j),'r.','markersize',6);
    plot([0 2],[0 2],'r:'); hold('off');
    
    if j>4&&SSEXplusY(j)>SSEXplusY(j-1)
        break
    end
end
[~,mini] = min(SSEXplusY);
PesToFlow.NfoldCcw_best = NfoldCcw(mini);

if mini>1&mini<length(SSEXplusY)
    range = (mini-1:mini+1);
    [x_peak,y_peak] = PeakFitQuadratic(NfoldCcw(range),-SSEXplusY(range));
    PesToFlow.NfoldCcw_best = x_peak;
end


%% plotforshow
for y=1:size(OSA.ranges,1)
    
    I=find(Time>OSA.ranges(y,1)&Time<OSA.ranges(y,2));
    FlowEdiTemp = FlowEdi(I);
    PmusTestTemp = Pes1_filtered(I) + PesToFlow.NfoldCcw_best*Pcw(I); %conservative
    
    FlowTemp = Flow(I);
    PesTemp = Pes(I);
    TimeTemp = Time(I);
    FlowPesTemp = FlowPes(I);
    EdiTemp = Edi(I); EdiTemp(EdiTemp<0)=0; EdiTemp=EdiTemp.^EdiToFlow.bestPower;
    CPAPrange = [-2 2];
    Pmasktemp = Pmask(I);
    
    downsamplefactor=1;
    xdata.data = downsample(-PmusTestTemp,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    ratio = PesToFlow.Parameters(end,1)/EdiToFlow.Parameters(end,1)/Xmodel(mini);
    modelY = PesToVflowModel(ratio*EdiToFlow.Parameters(end,:),xdata); %is this use of 'x' right?
    
    figure(234);
    ax234(1)=subplot(3,2,1); plot(TimeTemp,modelY);
    ax234(2)=subplot(3,2,3); plot(TimeTemp,FlowEdiTemp);
    ax234(3)=subplot(3,2,5); plot(TimeTemp,FlowTemp);
    ax235(1)=subplot(3,2,2); plot(TimeTemp,PmusTestTemp);
    ax235(2)=subplot(3,2,4); plot(TimeTemp,-EdiTemp*ratio);
    ax235(3)=subplot(3,2,6); plot(TimeTemp,PmusTestTemp,'r',TimeTemp,PesTemp,'b');
    linkaxes(ax234);
    %linkaxes([ax234 ax235],'x');
    
    %pause
end

%% Estimate respiratory impedance with NfoldPcw (Pmus = Pes + Ecw.Vol)

NfoldPcw = PesToFlow.NfoldCcw_best;

parameters_x=[];
for x=1:size(PesToFlow.ranges,1)
    
    %get data
    I=find(Time>PesToFlow.ranges(x,1)&Time<PesToFlow.ranges(x,2));
    
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    
    xdata_filtered = Pes1_filtered(I) + NfoldPcw*Pcw(I);
    
    minimum_figs=1;
    [I_edi,~,~,~,~]=Vflowanalysis2(FlowEdi(I),Time(I),dt,minimum_figs);
    baselinePes = median(xdata_filtered(I_edi.starti));
    
    downsamplefactor=1;
    xdata.data = downsample(-xdata_filtered+baselinePes,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    ydata = downsample(Flow(I),downsamplefactor);
    
    set_alpha=0;
    lsqoptions=optimset('display','off');
    parameters=[6 20 0 0 0]; %guess: R E initialVL Vleak baselinePesErr %alpha = Younes exponent for length-tension, 0 is no effect
    
    lower=[0.1 0.1 -0.01 0 -0.0001]; %R E initialVL Vleak baselinePesErr
    upper=[60 300 0.01 0.5 0]; %R E initialVL Vleak baselinePesErr
    
    [modelYguess] = PesToVflowModel(parameters,xdata);
    
    for i=1:5
        [parameters,resnormC1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@PesToVflowModel,parameters,xdata,ydata,lower,upper,lsqoptions);
    end
    
    [modelY] = PesToVflowModel(parameters,xdata);
    
    rsquared_x(x) = sum((modelY-ydata).^2)/sum((ydata-mean(ydata)).^2);
    
    figure()
    ax(1)=subplot(3,2,1)
    plot(xdata.time,ydata,'k')
    hold('on')
    plot(xdata.time,modelY,'r')
    hold('off');
    ax(2)=subplot(3,2,3)
    plot(xdata.time,cumsum(ydata-parameters(5))*dt,'k')
    hold('on')
    plot(xdata.time,cumsum(modelY-parameters(5))*dt,'r')
    hold('off');
    ax(3)=subplot(3,2,5)
    plot(xdata.time,xdata.data,'k')
    linkaxes(ax,'x')
    
    subplot(2,2,2)
    plot(ydata,modelY,'k');
    
    subplot(2,2,4);
    plot(cumsum(ydata-parameters(5))*dt,cumsum(modelY-parameters(5))*dt,'k')
    
    parameters_x(x,:)=parameters;
    
end

tauRC = parameters_x(:,1)./parameters_x(:,2);
tauRC_median = median(tauRC);
ZRC = abs(parameters_x(:,1)+1i*parameters_x(:,2));

parameters_median = median(parameters_x);
parameters_x = [parameters_x;parameters_median]

PesToFlow.PmusParametersOpt = parameters_x;
%% Test: Use parameters from X on segment Y (PmusOpt)

PesToFlow.PmusRsquaredOpt=[];

for i=1:size(PesToFlow.PmusParametersOpt,1)-1
    y=i;
    
    I=find(Time>PesToFlow.ranges(y,1)&Time<PesToFlow.ranges(y,2));
    
    xdata_filtered = Pes1_filtered(I) + NfoldPcw*Pcw(I);
    
    minimum_figs=1;
    [I_edi,~,~,~,~]=Vflowanalysis2(FlowEdi(I),Time(I),dt,minimum_figs);
    baselinePes = median(xdata_filtered(I_edi.starti));
    
    downsamplefactor=1;
    xdata.data = downsample(-xdata_filtered+baselinePes,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    
    ydata = downsample(Flow(I),downsamplefactor);
    
    %tempparameters = PesToFlow.PmusParameters2(x,:);
    %tempparameters(2) = tempparameters(2)+1/CCWCPAP.Ccw;
    modelY = PesToVflowModel(PesToFlow.PmusParametersOpt(end,:),xdata); %is this use of 'x' right?
    
    SSE = sum(((modelY-mean(modelY))-(ydata-mean(ydata))).^2);
    SStot = sum((modelY-mean(modelY)).^2);
    
    
    figure()
    ax(1)=subplot(3,2,1)
    plot(xdata.time,ydata-mean(ydata),'k');
    hold('on');
    plot(xdata.time,modelY-mean(modelY),'r');
    hold('off');
    ax(2)=subplot(3,2,3)
    plot(xdata.time,cumsum(ydata-mean(ydata))*dt,'k');
    hold('on')
    plot(xdata.time,cumsum(modelY-mean(modelY))*dt,'r');
    hold('off');
    ax(3)=subplot(3,2,5)
    plot(xdata.time,xdata.data,'k')
    linkaxes(ax,'x')
    
    subplot(2,2,2)
    plot(ydata,modelY,'k');
    
    subplot(2,2,4);
    plot(cumsum(ydata-parameters(4))*dt,cumsum(modelY-parameters(4))*dt,'k')
    
    PesToFlow.PmusRsquaredOpt(i) = 1-SSE/SStot;
end

close all



%% Calculate Vdrive (FlowPmusOpt) for whole night

PesToFlow.useparametersetOpt=size(PesToFlow.PmusParametersOpt,1);
parametersOSA=PesToFlow.PmusParametersOpt(PesToFlow.useparametersetOpt,:)
parametersOSA(5)=0;
%parametersOSA(2)=parametersOSA(2) + 1/CCWCPAP.Ccw;
Frescale=1;
parametersOSA(1:2)=parametersOSA(1:2)*Frescale;

Pmus1_filtered2 = Pes1_filtered + NfoldPcw*Pcw; %%%%%%%%%%%%%%%%%% Temp fudge factor

%Pmus1_filtered = filtfilt(B_butter0,A_butter0,Pes - 1/CCWCPAP.Ccw*cumsum(Flow)*dt);

xdata.data = -Pmus1_filtered2 ; %downsample(Pes(I),downsamplefactor);
xdata.time = Time;
xdata.dt = dt;
%clear Pmus1_filtered
parametersOSA(3) = -xdata.data(1)/parametersOSA(2);
FlowPmusOpt = PesToVflowModel(parametersOSA,xdata);
clear xdata



%%
save(filenameanddir,'PesToFlow','-append');
%% Plot all
if 1
    figure(101);
    if 1
        ax101(1)=subplot(3,1,1); plot(Time,[Flow FlowPes FlowPmus FlowPmusOpt FlowEdi]);
        magicratio = PesToFlow.PmusParameters(end,1)/EdiToFlow.ParametersSqrt(end,1);
        ax101(2)=subplot(3,1,2); plot(Time,[0*Flow Pes Pmus1_filtered -magicratio*(Edi.^0.5)]); %EdiToFlow.ParametersSqrt(end,1)
    elseif 0
        ax101(1)=subplot(3,1,1); plot(Time,[Flow FlowPes FlowPmus2 FlowEdiPower75]);
        magicratio = PesToFlow.PmusParameters2(end,1)/EdiToFlow.ParametersPower75(end,1);
        ax101(2)=subplot(3,1,2); plot(Time,[0*Flow Pes Pmus1_filtered -magicratio*(Edi.^0.75)]); %EdiToFlow.ParametersSqrt(end,1)
    elseif 0
        ax101(1)=subplot(3,1,1); plot(Time,[Flow FlowPes FlowPmus FlowEdiPower75]);
        magicratio = PesToFlow.PmusParameters2(end,1)/EdiToFlow.ParametersPower75(end,1);
        ax101(2)=subplot(3,1,2); plot(Time,[0*Flow Pes Pmus1_filtered -magicratio*(Edi.^0.75)]); %EdiToFlow.ParametersSqrt(end,1)
    else
        ax101(1)=subplot(3,1,1); plot(Time,[Flow FlowPes FlowPmus2 FlowEdiSqrt]);
        magicratio = PesToFlow.PmusParameters2(end,1)/EdiToFlow.ParametersSqrt(end,1);
        ax101(2)=subplot(3,1,2); plot(Time,[0*Flow Pes Pmus1_filtered -magicratio*(Edi.^0.5)]); %EdiToFlow.ParametersSqrt(end,1)
    end
    ax101(3)=subplot(3,1,3); plot(Time,Pmask);
    linkaxes(ax101,'x');
end

%%
clear Pmus1_filtered Pcw Pmus2_filtered
%% Compare and Plot Pmus vs Edi

figure();
x=size(parameters_x,1);
ydata = [];
xdata = [];
zdata = [];

for i=1:size(parameters_x,1)-1
    y=i;
    
    I=find(Time>PesToFlow.ranges(y,1)&Time<PesToFlow.ranges(y,2));
    [I_edi,~,~,~,~]=Vflowanalysis2(FlowEdi(I),Time(I),dt,minimum_figs);
    tempPes = Pes(I);
    baselinePes = median(tempPes(I_edi.starti-1));
    
    tempPmus = Pes(I) - 1/CCWCPAP.Ccw*cumsum(Flow(I)-mean(Flow(I)))*dt;
    baselinePmus = median(tempPmus(I_edi.starti-1));
    % subplot(1,2,1); plot(Edi(I),Pes(I)-baselinePes,'.','markersize',1); hold('on');
    % subplot(1,2,2); plot(Edi(I),tempPmus-baselinePmus,'.','markersize',1); hold('on');
    % pause
    xdata = [xdata;Pes(I)-baselinePes];
    ydata = [ydata;tempPmus-baselinePmus];
    zdata = [zdata;Edi(I)];
end
subplot(2,2,1); plot(zdata,xdata,'.','markersize',1);
subplot(2,2,2); plot(zdata,ydata,'.','markersize',1);
subplot(2,2,3); plot(zdata.^0.5,xdata,'.','markersize',1);
subplot(2,2,4); plot(zdata.^0.5,ydata,'.','markersize',1);

powers = [0.5:0.01:1];

for i=1:length(powers)
    Rpes(i) = corr(zdata.^powers(i),xdata);
    Rpmus(i) = corr(zdata.^powers(i),ydata);
end

[~,temp] = max(abs(Rpes)); power_maxcorrPes = powers(temp);
[~,temp] = max(abs(Rpmus)); power_maxcorrPmus = powers(temp);

EdiToFlow.power_maxcorrPes = power_maxcorrPes;
EdiToFlow.power_maxcorrPmus = power_maxcorrPmus;

save(filenameanddir,'EdiToFlow','-append');
save(filenameanddir,'PesToFlow','-append');

%% Find Spontaneous OSA with Edi

% Move left and right, select range with left click inside window
figure(2);
Xminute=7;
Yminute=0.5;
dsf=5;
X=4;
ax2(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(EpochsXHz,dsf),'b',downsample(Time,dsf),downsample(EventsArXHz,dsf),'r');
ax2(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Pmask,dsf));
ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Edi,dsf));
ax2(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(Flow,dsf)); ylim([-1.1 1.1]);
linkaxes(ax2,'x');
upperX=Xminute*60;
lowerX=0;
xlim([lowerX upperX]);
rangedata=[];
while 1
    clear x1 x2 y1 y2 button
    [x1,y1,button]=ginput(1);
    if button==3
        break
    end
    y1
    if y1>0
        deltaX=Xminute*60;
    else
        deltaX=Yminute*60;
    end
    if x1>upperX
        upperX=upperX+deltaX;
        lowerX=lowerX+deltaX;
    elseif x1<lowerX
        upperX=upperX-deltaX;
        lowerX=lowerX-deltaX;
    else
        rangedata=[rangedata;lowerX,upperX];
    end
    xlim([lowerX upperX]);
end
%%
OSA.ranges=rangedata;
%OSA.ranges=[rangedata;OSA.ranges];

%% Get CPAP and sleep proportions for each OSA window

OSA_rangesi = round(((OSA.ranges-Time(1))/dt +1));
OSA.CPAP=[];
for i=1:size(OSA_rangesi,1)
    I = OSA_rangesi(i,1):OSA_rangesi(i,2);
    OSA.CPAP(i,:) = prctile(Pmask(I),50);
    OSA.FREM(i,:) = sum(EpochsXHz(I)==5)/length(I);
    OSA.FWake(i,:) = sum(EpochsXHz(I)==0)/length(I);
    OSA.FAr(i,:) = sum(EventsArXHz(I)==0)/length(I);
    OSA.FNREM(i,:) = sum(EpochsXHz(I)>0&EpochsXHz(I)<5)/length(I);
    OSA.FN1(i,:) = sum(EpochsXHz(I)>0&EpochsXHz(I)<5)/length(I);
    OSA.FN2(i,:) = sum(EpochsXHz(I)>0&EpochsXHz(I)<5)/length(I);
    OSA.FN3(i,:) = sum(EpochsXHz(I)>0&EpochsXHz(I)<5)/length(I);
end

%% Remove duplicate channels
for i=length(channels):-1:1
    for j=1:i-1
        if strcmp(channels{i},channels{j})|isempty(channels{i})
            channels(i)=[];
            break
        end
    end
end

%% ---------------- Export segments for further analysis

directory2 = [directory 'OSAsegments'];
currentdir = cd;
directory2(1)=currentdir(1);
analysis_list = {'CCWCPAP','SpO2wake','Pcrit','ArThresPes','ArThresEpi','ArthresEdi','Vactive','Vpassive','Varousal','LGediVarousal'} %Add more to this list
savestring2 = [];
for i=1:length(analysis_list)
    if exist(analysis_list{i})
        savestring2 = [savestring2 analysis_list{i} ' '];
    end
end

for i=1:size(OSA.ranges,1);
    segmentfilename = [filename '_' num2str(OSA.ranges(i,1))];
    left = round((OSA.ranges(i,1)-Time(1))/dt)+1;
    right = round((OSA.ranges(i,2)-Time(1))/dt)+1;
    I=left:right;
    savestring=['channels Time1 Veupnea EdiToFlow OSA ' savestring2];
    for m=1:length(channels)
        if ~isempty(channels{m})
            eval([channels{m} '1=' channels{m} '(I);'])
            savestring = [savestring ' ' channels{m} '1'];
        end
    end
    Time1=Time(I);
    cd(directory2)
    eval(['save ' segmentfilename ' ' savestring]);
end
cd(currentdir);

%% Save
if saveontherun
    %     save(filenameanddir,'OSA','-append');
    save(filenameanddir,'OSA','EdiToFlow','-append');
end

%% obstructedXHz and centralapneasXHz for single window
% obstructedXHz = 0*xdata.time;
% centralapneasXHz = 0*xdata.time;
% for i=1:length(Evts.starttimes)
%     if Evts.starttimes(i)>xdata.time(end)||Evts.endtimes(i)<xdata.time(1)
%         continue
%     end
%     if sum(Evts.codes(i)==[2 4 5 7 9])==0
%         continue
%     end
%     lefti=round((Evts.starttimes(i)-xdata.time(1))/dt+1);
%     righti=round((Evts.endtimes(i)-xdata.time(1))/dt+1);
%     if lefti<i&&righti>=1
%         lefti=1;
%     end
%     if lefti<=length(xdata.time)&&righti>length(xdata.time)
%         righti=length(xdata.time);
%     end
%     obstructedXHz(lefti:righti)=1;
% end
% for i=1:length(Evts.starttimes)
%     if Evts.starttimes(i)>xdata.time(end)||Evts.endtimes(i)<xdata.time(1)
%         continue
%     end
%     if Evts.codes(i)~=3
%         continue
%     end
%     lefti=round((Evts.starttimes(i)-xdata.time(1))/dt+1);
%     righti=round((Evts.endtimes(i)-xdata.time(1))/dt+1);
%     if lefti<i&&righti>=1
%         lefti=1;
%     end
%     if lefti<=length(xdata.time)&&righti>length(xdata.time)
%         righti=length(xdata.time);
%     end
%     centralapneasXHz(lefti:righti)=1;
% end
%
% %    plot(xdata.time,1-obstructedXHz,'r');
% %    linkaxes(ax6,'x')

%% Get FlowEdi (Vdrive) XHz within OSA ranges

x=2
%1710 #2 is great!
%get data
I=find(Time>OSA.ranges(x,1)&Time<OSA.ranges(x,2));

figure(5)
X=5;
ax5(1)=subplot(X,1,1); plot(Time(I),EpochsXHz(I),'b',Time(I),EventsArXHz(I),'r');
ax5(2)=subplot(X,1,2); plot(Time(I),Pmask(I));
ax5(3)=subplot(X,1,3); plot(Time(I),Edi(I));

modelY = FlowEdi(I); %%%%%%%%%%%%%
ax5(4)=subplot(X,1,4); plot(Time(I),modelY,'r'); hold('on');
ax5(4)=subplot(X,1,4); plot(Time(I),Flow(I),'b');

linkaxes(ax5,'x')

%% Estimate Vdrive using modelling

xdata.time = Time(I);

minimum_figs=0
[I_,VI,VTi,Ti,Te]=Vflowanalysis1(modelY,Time(I),dt,minimum_figs)

CPAP_OSA = median(Pmask(I(1)+I_.starti));
%figure(6);
%ax6(2)=subplot(4,1,2); plot(xdata.time(I_.starti),modelY(I_.starti),'r.'); hold('on');

[I2,VI2,VTi2,Ti2,Te2,leak]=Vflowanalysis_knownI(Flow(I)-mean(Flow(I)),Time(I),dt,minimum_figs,I_);

%figure(6);
%ax6(2)=subplot(4,1,2); plot(xdata.time(I2.starti),ydata(I2.starti),'ro');

%ax6(3)=subplot(4,1,3);

figure(601)
stairs(xdata.time(I_.starti),VI2,'k');  hold('on');
stairs(xdata.time(I_.starti),VI,'r');
plot(xdata.time(I_.starti),VI*0+mean(VI2),'k:');
plot(xdata.time(I_.starti),VI*0+Veupnea.mean/60,'g:');
%ax6(4)=subplot(4,1,4);



%% Events in breath domain
obstructedXHz = ObstructedXHz(I);
centralapneasXHz = CentralapneasXHz(I);

E1=0*xdata.time(I_.starti);
Ar=0*E1;
CentralApneas=0*E1
for i=1:length(VI)
    E1(i)=1-min(obstructedXHz(I_.starti(i):I_.endi(i)));
    Ar(i)=max(EventsArXHz(I(1)-1+(I_.starti(i):I_.endi(i))));
    CentralApneas(i)=min(centralapneasXHz(I_.starti(i):I_.endi(i)));
end
Ar(E1==0|CentralApneas==1)=0;

Ttot_n = Ti+Te;
Ttot_previous = [median(Ttot_n) Ttot_n(1:end-1)];

%% Flowlimitationindex: Teschler
ydata = Flow(I);
figure(99);
subplot(1,1,1);
hold('on');
plot(xdata.time,ydata-leak,'color',[0.7 0.7 0.7]);
clear fli
for i=1:length(I2.starti)
    Nrows = ceil((length(I2.starti)*0.5)^0.5);
    Ncols = 2*Nrows;
    rangetemp = I2.starti(i):I2.midi(i);
    %meanflowtemp = mean(ydata((I2.starti(i):I2.endi(i))));
    FlowData = ydata(rangetemp)-leak;
    fli(i) = FlatteningIndexTeschler96(FlowData);
    if fli(i)<0.15
        colorplot = [1 0 0];
    else
        colorplot = [0 1 0];
    end
    plot(xdata.time(rangetemp),ydata(rangetemp)-leak,'color',colorplot)
    plot(xdata.time(rangetemp),mean(ydata(rangetemp))+0*ydata(rangetemp),'k:');
    text(xdata.time(rangetemp(1)),-0.05,num2str(fli(i),2),'fontsize',7)
end
hold('off');

%% Flowlimitationindex: inverted parabola

figure(95);
subplot(1,1,1);
hold('on');
plot(xdata.time,ydata-leak,'color',[0.7 0.7 0.7]);
clear fli_ip
for i=1:length(I2.starti)
    rangetemp = (I2.starti(i)):(I2.midi(i));
    %meanflowtemp = mean(ydata((I2.starti(i):I2.endi(i))));
    FlowData = ydata(rangetemp)-leak;
    Thres = 95;
    [c,indices,y] = FindBestInvertedParabola(FlowData,Thres);
    fli_ip(i,:) = indices;
    %fli_ip(i,1)>2; fli_ip(i,1)<0.30
    if fli_ip(i,1)<0.9
        colorplot = [1 0 0]; %red for FL
    else
        colorplot = [0 1 0]; %green for nFL
    end
    plot(xdata.time(rangetemp),ydata(rangetemp)-leak,'color',colorplot);
    plot(xdata.time(rangetemp),y,':','color',colorplot);
    text(xdata.time(rangetemp(1)),-0.05,num2str(fli_ip(i,1),2),'fontsize',7);
    %text(xdata.time(rangetemp(1)),-0.10,num2str(100*(VI(i)-VI2(i))/tempeupnea,2),'fontsize',7);
    text(xdata.time(rangetemp(1)),-0.10,num2str(100*(VI2(i)/VI(i)),2),'fontsize',7);
end
hold('off');

% Evaluate FL indices
%ydataplot = 100*(VI-VI2)/tempeupnea;
ydataplot = 100*(VI2./VI);

figure(100);
tempeupnea = Veupnea.mean/60;
subplot(1,5,1); plot(fli,ydataplot,'.'); %intended - actual flow; higher = greater discrepancy
%>30% should be an event.
ylabel('VE/Vdrive(%)'); xlabel('Teschler FLI');

subplot(1,5,2); plot(Ti2./(Ti2+Te2),ydataplot,'.');
ylabel('VE/Vdrive(%)'); xlabel('Ti/Ttot');

subplot(1,5,3); plot(fli_ip(:,1),ydataplot,'.');
ylabel('VE/Vdrive(%)'); xlabel('Farea');

subplot(1,5,4); plot(fli_ip(:,2),ydataplot,'.');
ylabel('VE/Vdrive(%)'); xlabel('Fmidinsp');

subplot(1,5,5); plot(fli_ip(:,3),ydataplot,'.');
ylabel('VE/Vdrive(%)'); xlabel('Fmindelta');

% figure(105);
% plot(fli_ip(:,4),ydataplot,'.');
% ylabel('VE/Vdrive(%)'); xlabel('Rsquared');

%% LG from measured Vdrive (Edi)
warning('off');
polyfitorder=3;
if 0
    Veupnea2=mean(VI2);%Veupnea.mean/60;
else
    Veupnea2=Veupnea.mean/60;%;
end

VraOn=1;
VEinput=(VI2'-Veupnea2)/Veupnea2;
Veupnea_input=1;
Vdriveknowninput=(VI'-Veupnea2)/Veupnea2;

Data2=[VEinput E1 Ar xdata.time(I_.starti) CentralApneas Ttot_n' Ttot_previous' Vdriveknowninput];
[Error2,Error2prepoly,Vchem2,Varousal2,LoopGainParameters2,BestParameters2,BestSSres2,FitQuality2,i_12,i_end2] = FindBestModelParametersEdi(Data2,VraOn,Veupnea_input,polyfitorder)


LGss_from_means = -mean(Vdriveknowninput)/mean(VEinput);

f=0;
%X=(LoopGainParameters2(7)/LoopGainParameters2(6))^2;
%tau_min = ((X-1)/(4*pi^2*(1-4*X)))^0.5;
LGEdif = LoopGainParameters2(1)/(1+(2*pi*LoopGainParameters2(2)*f/60)^2)^0.5


%% Real Vdrive figure

figure(55);
set(gcf,'color',[1 1 1])
ax55(1)=subplot(2,1,1); hold('on');
set(gca,'XColor',[1 1 1]);
set(gca,'tickdir','out');

%make traces for shaded regions

dt_new=0.25;
startT=ceil(xdata.time(I_.starti(1)));
endT=ceil(xdata.time(I_.starti(end)));
time_dt=startT:dt_new:endT;
E1_rs=0*time_dt;
for i=1:length(time_dt);
    E1_rs(i) = 1-E1(find(xdata.time(I_.starti)<=time_dt(i),1,'last'));
end
A1_rs=0*time_dt;
for i=1:length(time_dt);
    A1_rs(i) = Ar(find(xdata.time(I_.starti)<=time_dt(i),1,'last'));
end
C1_rs=0*time_dt;
for i=1:length(time_dt);
    C1_rs(i) = CentralApneas(find(xdata.time(I_.starti)<=time_dt(i),1,'last'));
end
E1_rs=E1_rs&(~C1_rs);

tempsize=1; tempoffset=-tempsize/2;
tempy=tempsize*(E1_rs)+tempoffset;
tempxx=[time_dt fliplr(time_dt)];
tempyy=[tempy tempoffset+zeros(1,length(tempy))];
fill(tempxx,tempyy,[210 210 255]/255,'EdgeColor','none');
hold('on')

tempsize=1; tempoffset=-tempsize/2;
tempy=tempsize*(C1_rs)+tempoffset;
tempxx=[time_dt fliplr(time_dt)];
tempyy=[tempy tempoffset+zeros(1,length(tempy))];
fill(tempxx,tempyy,[255 198 198]/255,'EdgeColor','none');

tempsize=0.5; tempoffset=2;
tempy=tempsize*(A1_rs)+tempoffset;
tempxx=[time_dt fliplr(time_dt)];
tempyy=[tempy tempoffset+zeros(1,length(tempy))];
fill(tempxx,tempyy,[140 255 140]/255,'EdgeColor','none');

temp_leak=[];
for i=1:length(I2.starti)
    temp_leak(i)=mean(Flow(I(1)+(I2.starti:I2.endi)));
end
leak=median(temp_leak);


plot(Time(I),Flow(I)-leak);
distance = -2;
plot(Time(I),modelY+distance,'r');
plot(Time(I),Edi(I)/15+distance*2.5,'r');

ax55(2)=subplot(2,1,2); hold('on');


tempsize=0.5; tempoffset=0.375;
tempy=tempsize*(E1_rs)+tempoffset;
tempxx=[time_dt fliplr(time_dt)];
tempyy=[tempy tempoffset+zeros(1,length(tempy))];
fill(tempxx,tempyy,[210 210 255]/255,'EdgeColor','none');
hold('on')
tempsize=0.5; tempoffset=0.375;
tempy=tempsize*(C1_rs)+tempoffset;
tempxx=[time_dt fliplr(time_dt)];
tempyy=[tempy tempoffset+zeros(1,length(tempy))];
fill(tempxx,tempyy,[255 198 198]/255,'EdgeColor','none');

tempsize=0.25; tempoffset=2;
tempy=tempsize*(A1_rs)+tempoffset;
tempxx=[time_dt fliplr(time_dt)];
tempyy=[tempy tempoffset+zeros(1,length(tempy))];
fill(tempxx,tempyy,[140 255 140]/255,'EdgeColor','none');

stairs(xdata.time(I_.starti),VI2/Veupnea2,'color',[0.4 0.4 0.7]);  hold('on');
stairs(xdata.time(I_.starti),VI/Veupnea2,'r');
plot(xdata.time(I_.starti),VI*0+mean(VI2)/Veupnea2,'k:');
plot(xdata.time(I_.starti),VI*0+Veupnea.mean/60/Veupnea2,'g:');


plot(xdata.time(I_.starti),Vchem2+1,'k--');
plot(xdata.time(I_.starti),Vchem2+1+Varousal2,'k');

%plot(xdata.time(I_.starti),Varousal2,'r');
%plot(xdata.time(I_.starti),Error2,'b:');
set(gca,'tickdir','out');
linkaxes(ax55,'x')


%% Use VE to get Vdrive from OSA only
manualscoringtouchups=0;
Veupnea_=mean(VI2);
%Veupnea_=Veupnea.mean/60;
h34=figure(34);
while 1
    
    % LG from OSA
    warning('off');
    polyfitorder=3;
    VraOn=1;
    Data=[(VI2'-Veupnea_)/Veupnea_ E1   Ar  xdata.time(I_.starti)  CentralApneas  Ttot_n'  Ttot_previous' ];
    [Error,Vchem,Var,LoopGainParameters,BestParameters,BestSSres,FitQuality,i_1,i_end] = FindBestModelParameters(Data,VraOn,1,polyfitorder)
    
    %     figure(601)
    %     hold('on');
    %     plot(xdata.time(I_.starti),(Vchem+1)*Veupnea_,'b');  hold('on');
    %     plot(xdata.time(I_.starti),(Vchem+Varousal+1)*Veupnea_,'g');
    %
    f=0;
    %X=(LoopGainParameters2(7)/LoopGainParameters2(6))^2;
    %tau_min = ((X-1)/(4*pi^2*(1-4*X)))^0.5;
    LGosaf = LoopGainParameters(1)/(1+(2*pi*LoopGainParameters(2)*f/60)^2)^0.5
    
    delete h34
    h34=figure(34);
    set(gcf,'color',[1 1 1])
    ax55(1)=subplot(2,1,1); hold('on');
    set(gca,'XColor',[1 1 1]);
    set(gca,'tickdir','out');
    
    %make traces for shaded regions
    
    dt_new=0.25;
    startT=ceil(xdata.time(I_.starti(1)));
    endT=ceil(xdata.time(I_.starti(end)));
    time_dt=startT:dt_new:endT;
    E1_rs=0*time_dt;
    for i=1:length(time_dt);
        E1_rs(i) = 1-E1(find(xdata.time(I_.starti)<=time_dt(i),1,'last'));
    end
    A1_rs=0*time_dt;
    for i=1:length(time_dt);
        A1_rs(i) = Ar(find(xdata.time(I_.starti)<=time_dt(i),1,'last'));
    end
    C1_rs=0*time_dt;
    for i=1:length(time_dt);
        C1_rs(i) = CentralApneas(find(xdata.time(I_.starti)<=time_dt(i),1,'last'));
    end
    E1_rs=E1_rs&(~C1_rs);
    
    tempsize=1; tempoffset=-tempsize/2;
    tempy=tempsize*(E1_rs)+tempoffset;
    tempxx=[time_dt fliplr(time_dt)];
    tempyy=[tempy tempoffset+zeros(1,length(tempy))];
    fill(tempxx,tempyy,[210 210 255]/255,'EdgeColor','none');
    hold('on')
    
    tempsize=1; tempoffset=-tempsize/2;
    tempy=tempsize*(C1_rs)+tempoffset;
    tempxx=[time_dt fliplr(time_dt)];
    tempyy=[tempy tempoffset+zeros(1,length(tempy))];
    fill(tempxx,tempyy,[255 198 198]/255,'EdgeColor','none');
    
    tempsize=0.5; tempoffset=2;
    tempy=tempsize*(A1_rs)+tempoffset;
    tempxx=[time_dt fliplr(time_dt)];
    tempyy=[tempy tempoffset+zeros(1,length(tempy))];
    fill(tempxx,tempyy,[140 255 140]/255,'EdgeColor','none');
    
    temp_leak=[];
    for i=1:length(I2.starti)
        temp_leak(i)=mean(Flow(I(1)+(I2.starti:I2.endi)));
    end
    leak=median(temp_leak);
    
    
    plot(Time(I),Flow(I)-leak);
    distance = -2;
    plot(Time(I),modelY+distance,'r');
    plot(Time(I),Edi(I)/15+distance*2.5,'r'); hold('off');
    
    ax55(2)=subplot(2,1,2); hold('on');
    
    
    tempsize=0.5; tempoffset=0.375;
    tempy=tempsize*(E1_rs)+tempoffset;
    tempxx=[time_dt fliplr(time_dt)];
    tempyy=[tempy tempoffset+zeros(1,length(tempy))];
    fill(tempxx,tempyy,[210 210 255]/255,'EdgeColor','none');
    hold('on')
    tempsize=0.5; tempoffset=0.375;
    tempy=tempsize*(C1_rs)+tempoffset;
    tempxx=[time_dt fliplr(time_dt)];
    tempyy=[tempy tempoffset+zeros(1,length(tempy))];
    fill(tempxx,tempyy,[255 198 198]/255,'EdgeColor','none');
    
    tempsize=0.25; tempoffset=2;
    tempy=tempsize*(A1_rs)+tempoffset;
    tempxx=[time_dt fliplr(time_dt)];
    tempyy=[tempy tempoffset+zeros(1,length(tempy))];
    fill(tempxx,tempyy,[140 255 140]/255,'EdgeColor','none');
    
    stairs(xdata.time(I_.starti),VI2/Veupnea2,'color',[0.4 0.4 0.7]);  hold('on');
    stairs(xdata.time(I_.starti),VI/Veupnea2,'r');
    plot(xdata.time(I_.starti),VI*0+mean(VI2)/Veupnea2,'k:');
    plot(xdata.time(I_.starti),VI*0+Veupnea.mean/60/Veupnea2,'g:');
    
    
    plot(xdata.time(I_.starti),Vchem+1,'k--');
    plot(xdata.time(I_.starti),Vchem+1+Var,'k');  hold('off');
    
    set(gca,'tickdir','out');
    linkaxes(ax55,'x')
    
    if manualscoringtouchups
        disp('click left of data if no touch ups needed, or left then right to add events in range, or right then left to remove events in range');
        [t(1),y(1),~]=ginput(1);
        if t(1)<min(xdata.time(I_.starti(1)))
            manualscoringtouchupsincomplete=0;
        else
            [t(2),y(1),button]=ginput(1);
            manualscoringtouchupsincomplete=1;
            if y(1)<2
                if button==1
                    if t(1)<t(2)
                        E1(xdata.time(I_.starti)>t(1)&xdata.time(I_.starti)<t(2))=0; %add events...
                    elseif t(1)>t(2)
                        E1(xdata.time(I_.starti)>t(2)&xdata.time(I_.starti)<t(1))=1; %remove events...
                    end
                elseif button>1 %right click
                    if t(1)<t(2)
                        CentralApneas(xdata.time(I_.starti)>t(1)&xdata.time(I_.starti)<t(2))=1; %add events...
                    elseif t(1)>t(2)
                        CentralApneas(xdata.time(I_.starti)>t(2)&xdata.time(I_.starti)<t(1))=0; %remove events...
                    end
                end
            else
                if t(1)<t(2)
                    Ar(xdata.time(I_.starti)>t(1)&xdata.time(I_.starti)<t(2))=1; %add events...
                elseif t(1)>t(2)
                    Ar(xdata.time(I_.starti)>t(2)&xdata.time(I_.starti)<t(1))=0; %remove events...
                end
                
                
            end
        end
    else
        manualscoringtouchupsincomplete=0;
    end
    
    if manualscoringtouchupsincomplete==0
        break
    end
    
end

LGsummary = [LoopGainParameters LGosaf FitQuality(2) LoopGainParameters2 LGEdif FitQuality2(2)]';


%% Real Vdrive vs VE

%Veupnea2=Veupnea.mean/60;

figure(110)
ArNext1=[Ar(2:end);NaN];
ArNext2=[ArNext1(2:end);NaN];
ArPrev1=[NaN;Ar(1:end-1)];
ArPrev2=[NaN;ArPrev1(1:end-1)];
postAr1=(Ar==0&ArPrev1==1);
postAr2=(ArPrev1==0&ArPrev2==1);
preAr1=(Ar==0&ArNext1==1);
preAr2=(ArNext1==0&ArNext2==1);
startAr=(Ar==1&ArPrev1==0);
plot(VI/Veupnea2,VI2/Veupnea2,'.','markersize',16); hold('on');
plot(VI(Ar==1)/Veupnea2,VI2(Ar==1)/Veupnea2,'r.','markersize',16);
plot(VI(startAr==1)/Veupnea2,VI2(startAr==1)/Veupnea2,'.','markersize',16,'color',[0.8 0 0]);
plot(VI(postAr1==1)/Veupnea2,VI2(postAr1==1)/Veupnea2,'.','markersize',16,'color',[1 0.4 0.4]);
plot(VI(postAr2==1)/Veupnea2,VI2(postAr2==1)/Veupnea2,'.','markersize',16,'color',[1 0.6 0.6]);
plot(VI(preAr1==1)/Veupnea2,VI2(preAr1==1)/Veupnea2,'.','markersize',16,'color',[0.6 1 0.6]);
plot(VI(preAr2==1)/Veupnea2,VI2(preAr2==1)/Veupnea2,'.','markersize',16,'color',[0.4 1 0.4]);
plot([0 max([VI VI2])]/Veupnea2,[0 max([VI VI2])]/Veupnea2,'k:');
xlim([0 max([VI VI2])]/Veupnea2);
ylim([0 max([VI VI2])]/Veupnea2);

%% VE, Vdrive using Edi
Feupneaerror = (Veupnea2*60)/Veupnea.mean;

Xdata2=1+Vchem2';
XdataGrad=gradient(Xdata2);

Xdata=VI/Veupnea2*Feupneaerror;
%Xdata=1+Vchem2';
%Xdata=1+Vchem2'+Varousal2';
Ydata=VI2/Veupnea2*Feupneaerror;
%criteria=1-((Ar==1)|(VI2'/Veupnea2>1));
%criteria=1-((Ar==1)|(VI2'/Veupnea2>1)|(E1==1));

%criteria=(1-((E1)|(VI2'/Veupnea2>1)))';

%criteria=(1-((Ar'==1)|(Ydata>1.2))|isnan(Xdata));
%criteria=(1-((Ar'==1)|(Ydata>1.2)|isnan(Xdata

%criteria=(1-((Ar'==1)|(Ydata>1.2)|isnan(Xdata)));

%criteria=(1-((Ar'==1)|isnan(Xdata)));
%criteria=(1-(isnan(Xdata)));
criteria=(1-((Ar'==1)|isnan(Xdata)|(Ydata>1)|XdataGrad<0)); %Excluding supra-eupneic VE and falling Vdrive

%criteria=(1-((Ar'==1)|isnan(Xdata)|(Ydata>1)));

%criteria=1-(isnan(Xdata));

%criteria=1-(isnan(Xdata)|(Ar'==1)|XdataGrad>0);

%criteria=1-((Ar'==1)|isnan(Xdata)|(Ydata>2)|(Xdata>5));

bounds = 0.1;
Vpassivemedianneareupnea=median(Ydata(criteria==1&Xdata<1+bounds&Xdata>1-bounds));
VpassivemedianneareupneaN=length(Ydata(criteria==1&Xdata<1+bounds&Xdata>1-bounds));

plotpoints=1;

ArthresEdiA=median(VI(preAr1==1)/Veupnea2);
ArthresEdiB=median(VI(startAr==1)/Veupnea2);

%[ArThresEdiOSA,i] = max([ArthresEdiA ArthresEdiB]);
%ArThresEdiOSA=(ArthresEdiA+ArthresEdiB)/2;
ArThresEdiOSA=ArthresEdiA;

global fixedslope fixedArthres
fixedslope=[];
fixedArthres=[];
[UAGSEMEdiOSA,UAGvalueEdiOSA,VcritSEMEdiOSA,VcritvalueEdiOSA]=VdriveVEmodel_run_equi(Xdata(criteria==1),Ydata(criteria==1),plotpoints);

hold('on');
criteria2 = 1-(E1|(Ar==1));
plot(Xdata(criteria==1&criteria2'==1),Ydata(criteria==1&criteria2'==1),'.','markersize',18,'color',[0.4 0.4 1]);
plot(Xdata(criteria==1&Ar'==1),Ydata(criteria==1&Ar'==1),'.','markersize',18,'color',[1 0 0]);
plot([ArThresEdiOSA ArThresEdiOSA],[0 1],'g','color',[0.4 1 0.4],'linewidth',2);
plot([1 ArThresEdiOSA],[1 1-(ArThresEdiOSA-1)/LGEdif],'g.-','markersize',18);

if 0
    criteria3=criteria;
    while 1
        plot(Xdata(criteria==1&criteria3==0),Ydata(criteria==1&criteria3==0),'ro','markersize',12);
        
        [x,y,button]=ginput(1)
        if button==3
            break
        end
        distance = (((Xdata)-x).^2 + ((Ydata)-y).^2).^0.5;
        [temp,i]=min(distance)
        criteria3(i)=0;
    end
    
    criteria=criteria3
    [UAGSEMEdiOSA,UAGvalueEdiOSA,VcritSEMEdiOSA,VcritvalueEdiOSA]=VdriveVEmodel_run(Xdata(criteria==1),Ydata(criteria==1),plotpoints);
    hold('on');
    criteria2 = 1-(E1|(Ar==1));
    plot(Xdata(criteria2==1),Ydata(criteria2==1),'.','markersize',18,'color',[0.4 0.4 1]);
    plot(Xdata(criteria==1&Ar==1),Ydata(criteria==1&Ar==1),'.','markersize',18,'color',[1 0 0]);
    plot([ArThresEdiOSA ArThresEdiOSA],[0 1],'g','color',[0.4 1 0.4],'linewidth',2);
    plot([1 ArThresEdiOSA],[1 1-(ArThresEdiOSA-1)/LGEdif],'g.-','markersize',18);
end

VactiveEdiOSA = VcritvalueEdiOSA + UAGvalueEdiOSA*(ArThresEdiOSA-1);
VarousalEdiOSA = 1-(ArThresEdiOSA-1)/LGEdif;
gapEdiOSA = VarousalEdiOSA-VactiveEdiOSA;
TraitsSummaryEdiOSA = [VcritvalueEdiOSA UAGvalueEdiOSA ArThresEdiOSA VactiveEdiOSA VarousalEdiOSA gapEdiOSA Vpassivemedianneareupnea]';



%% VE, Vdrive using modelling
Feupneaerror = (Veupnea2*60)/Veupnea.mean;

Ydata=VI2/Veupnea2;

Xdata=1+Vchem';
%Xdata=(1+Vchem')*(1+LoopGainParameters(8)); %assuming VRA is a fractional error, not physiological.
Xdata=1+Vchem'+LoopGainParameters(8); %assuming VRA is a fractional error, not physiological.


XdataGrad=gradient(Xdata);

%Xdata=1+Vchem';
%Xdata=1+Vchem'+Var';
%Ydata=VI2/Veupnea2;
%criteria=1-((Ar==1)|(VI2'/Veupnea2>1));
%criteria=1-((Ar==1)|(VI2'/Veupnea2>1)|(E1==1));

%criteria=(1-((E1)|(VI2'/Veupnea2>1)))';

%criteria=(1-((Ar'==1)|(Ydata>1.2))|isnan(Xdata));
%criteria=(1-((Ar'==1)|(Ydata>1.2)|isnan(Xdata

%criteria=(1-((Ar'==1)|(Ydata>1.2)|isnan(Xdata)));
criteria=(1-((Ar'==1)|isnan(Xdata)));
%criteria=(1-((Ar'==1)|isnan(Xdata)|(E1'==1)));
%criteria=(1-((Ar'==1)|isnan(Xdata)|(Ydata>1)|XdataGrad<0));
%criteria=1-(isnan(Xdata));
%criteria=1-((Ar'==1)|isnan(Xdata)|(Ydata>1)|(Xdata>2.2));

bounds = 0.1;
VpassivemedianneareupneaOSA=median(Ydata(criteria==1&Xdata<1+bounds&Xdata>1-bounds));
VpassivemedianneareupneaOSAN=length(Ydata(criteria==1&Xdata<1+bounds&Xdata>1-bounds));


plotpoints=1;

ArthresOSAA=median(Xdata(preAr1'==1&~isnan(Xdata)));
ArthresOSAB=median(Xdata(startAr'==1&~isnan(Xdata)));

%[ArThresOSA,i] = max([ArthresOSAA ArthresOSAB]);
ArThresOSA=ArthresOSAA;

global fixedslope fixedArthres
fixedslope=[];
fixedArthres=[];
[UAGSEMOSA,UAGvalueOSA,VcritSEMOSA,VcritvalueOSA]=VdriveVEmodel_run_equi(Xdata(criteria==1),Ydata(criteria==1),plotpoints);

hold('on');
criteria2 = 1-(E1|(Ar==1));
plot(Xdata(criteria2==1),Ydata(criteria2==1),'.','markersize',18,'color',[0.4 0.4 1]);
plot(Xdata(criteria==1&Ar'==1),Ydata(criteria==1&Ar'==1),'.','markersize',18,'color',[1 0 0]);
plot([ArThresOSA ArThresOSA],[0 1],'g','color',[0.4 1 0.4],'linewidth',2);
plot([1 ArThresOSA],[1 1-(ArThresOSA-1)/LGosaf],'g.-','markersize',18)

VactiveOSA = VcritvalueOSA + UAGvalueOSA*(ArThresOSA-1);
VarousalOSA = 1-(ArThresOSA-1)/LGosaf;
gapOSA = VarousalOSA-VactiveOSA;
TraitsSummaryOSA = [VcritvalueOSA UAGvalueOSA ArThresOSA VactiveOSA VarousalOSA gapOSA VpassivemedianneareupneaOSA]';

TraitsSummaryCombined = [LGsummary;TraitsSummaryOSA;TraitsSummaryEdiOSA;CPAP_OSA];

if ~exist('Varousal')
    Varousal.Feupnea=NaN;
end

if ~exist('ArThres')
    ArThres.Direct_Feupnea=NaN;
    ArThres.Indirect_Feupnea=NaN;
end
TraitsCPAPdrops = [Vpassive.Feupnea Vactive.Vcrit_Feupnea Varousal.Feupnea Veupnea.mean  ArThres.Direct_Feupnea ArThres.Indirect_Feupnea]';

%% Find ArthresEdi


%Right-click to exclude and move on
%Click left and right WITHIN the 5 breaths prior to arousal to include these breaths

clear VARx1 VARy1
figure(2)
Xminute=3;
Yminute=1;
dsf=5;
X=4;


ax2(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(EpochsXHz,dsf),'b',downsample(Time,dsf),downsample(EventsArXHz,dsf),'r');
if exist('WakeSleep'), hold('on'); plot(downsample(Time,dsf),downsample(WakeSleep,dsf),'color',[1 0.3 0.3]); hold('off'); end
ax2(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Pmask,dsf));
ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Flow,dsf)); ylim([-1.2 1.2]);
%ax2(5)=subplot(X,1,5); plot(downsample(Time,dsf),downsample(Edi,dsf));
ax2(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(FlowEdi,dsf)); ylim([-1.2 1.2]);
linkaxes(ax2,'x')

%arousals are at temp==1:
temp=diff(EventsArXHz);
I=find(temp==1);
i=1;
while i<=length(I)
    lowerX=Time(I(i))-Yminute*60;
    upperX=Time(I(i))+Yminute*60;
    xlim([lowerX upperX])
    
    lefti=(lowerX-Time(1))/dt+1;
    righti=(upperX-Time(1))/dt+1;
    
    xdata.data = Edi(lefti:righti); %downsample(Pes(I),downsamplefactor);
    xdata.time = Time(lefti:righti);
    xdata.dt = dt;
    
    linkaxes(ax2,'x')
    %[ymin,temptemp]=get(ax2(4),'ylim')
    %set(ax2(4),'ylim',[min([Pepi1trend(Time>lowerX&Time<upperX)+ArThresEpi.mean]-5) prctile(Pepi(Time>lowerX&Time<upperX),99.5)]);
    [VARx1(i,1),~,button(i)]=ginput(1);
    %     if button(i)~=3
    %         [VARx1(i,2),VARy1(i,2),button(i)]=ginput(1);
    %     end
    i=i+1;
end

VARx1(button==3)=[]; %keep 'approved' x values for left and right of flow trace.
VARx1 = [VARx1 VARx1+0.1];
ArthresEdi.ranges=VARx1;
ArthresEdi.button = button;
%% --save
if saveontherun
    save(filenameanddir,'ArthresEdi','-append');
end

%% Analyze ArThres using Edi and model
%load(filenameanddir,'ArthresEdi');
button=[];
chooseArs=0;
%ArthresEdi.comments = 'variable calibration';
ArthresEdi.Frescale = 1;
for i=1:size(ArthresEdi.ranges,1)
    
    deltaX=8;
    lowerX=ArthresEdi.ranges(i,1)-deltaX;
    upperX=ArthresEdi.ranges(i,2)+deltaX;
    lefti=(lowerX-Time(1))/dt+1;
    righti=(upperX-Time(1))/dt+1;
    
    modelY = ArthresEdi.Frescale*FlowEdi(lefti:righti); %downsample(Pes(I),downsamplefactor);
    
    xdata.time = Time(lefti:righti);
    xdata.dt = dt;
    %     if 0
    %         modelY = EdiToVflowModel(EdiToFlow.Parameters(end,:),xdata); % parameters_x(end,:)
    %     else
    %         modelY = FlowEdi(lefti:righti);
    %     end
    
    figure(3)
    %I=find(Time>lowerX&Time<upperX);
    subplot(2,1,1); plot(xdata.time,modelY,'k',[ArthresEdi.ranges(i,1) ArthresEdi.ranges(i,1)],[-1 1],'r',[ArthresEdi.ranges(i,2) ArthresEdi.ranges(i,2)],[-1 1],'r');
    vol1=cumsum(modelY-mean(modelY))*dt;
    thresh=std(vol1)/10;
    [min_list, max_list] = peakdet(-vol1,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
    %[indVI, indVE] = CheckVind(max_list, min_list);
    %TidalVol = V(indVI)-V(indVE);
    subplot(2,1,2); plot(xdata.time,vol1,'k',[ArthresEdi.ranges(i,1) ArthresEdi.ranges(i,1)],[-1 1],'r',[ArthresEdi.ranges(i,2) ArthresEdi.ranges(i,2)],[-1 1],'r');
    hold('on')
    subplot(2,1,2); plot(xdata.time(max_list(:,1)),vol1(max_list(:,1)),'k.');
    subplot(2,1,2); plot(xdata.time(min_list(:,1)),vol1(min_list(:,1)),'r.');
    %[ymin,temptemp]=get(ax2(4),'ylim')
    
    lefti=find(xdata.time(min_list(:,1))<ArthresEdi.ranges(i,1),1,'last');
    righti=find(xdata.time(min_list(:,1))>ArthresEdi.ranges(i,2),1,'first')-1;
    if isempty(lefti)
        lefti=1;
    end
    if isempty(righti)
        righti=lefti;
    end
    if righti<lefti
        righti=lefti;
    end
    subplot(2,1,2); plot([xdata.time(min_list(lefti,1)) xdata.time(min_list(lefti,1))],[-1 1],'r:',[xdata.time(min_list(righti+1,1)) xdata.time(min_list(righti+1,1))],[-1 1],'r:');
    
    xpoly = xdata.time(min_list(lefti:(righti+1),1));
    ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
    ypoly = vol1(min_list(lefti:(righti+1),1));
    ppoly = polyfit(xpoly,ypoly,1);
    Yvoldrift = polyval(ppoly,xdata.time);
    
    %repeat find end-exp and end-insp now that we have done some leak correction:
    vol3=cumsum(modelY-mean(modelY)-ppoly(1))*dt;
    thresh=std(vol3)/10;
    [min_list, max_list] = peakdet(-vol3,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
    
    subplot(2,1,2); plot(xdata.time(max_list(:,1)),vol1(max_list(:,1)),'ko');
    subplot(2,1,2); plot(xdata.time(min_list(:,1)),vol1(min_list(:,1)),'ro');
    
    lefti=find(xdata.time(min_list(:,1))<ArthresEdi.ranges(i,1),1,'last');
    righti=find(xdata.time(min_list(:,1))>ArthresEdi.ranges(i,2),1,'first')-1;
    if isempty(lefti)
        lefti=1;
    end
    if isempty(righti)||righti==0
        righti=size(min_list,1)-1;
    end
    subplot(2,1,2); plot([xdata.time(min_list(lefti,1)) xdata.time(min_list(lefti,1))],[-1 1],'r:',[xdata.time(min_list(righti+1,1)) xdata.time(min_list(righti+1,1))],[-1 1],'r:');
    
    xpoly = xdata.time(min_list(lefti:(righti+1),1));
    ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
    ypoly = vol1(min_list(lefti:(righti+1),1));
    ppoly = polyfit(xpoly,ypoly,1);
    Yvoldrift = polyval(ppoly,xdata.time);
    
    Yvol1_exp = ypoly(1:end-1); %end exp
    Yvol1_insp = vol1(max_list(lefti:(righti),1)); %end insp
    Yvoldrift_exp = polyval(ppoly,xpoly); Yvoldrift_exp(end)=[];
    Yvoldrift_insp = polyval(ppoly,xdata.time(max_list(lefti:(righti),1)));
    
    subplot(2,1,2); plot(xdata.time,Yvoldrift);
    hold('off');
    
    endexpvol=Yvol1_exp-Yvoldrift_exp;
    endinspvol=Yvol1_insp-Yvoldrift_insp;
    volbreath=endinspvol-endexpvol;
    ArthresEdi.data(i)=mean(volbreath)/ttot*60;
    ArthresEdi.CPAP(i)=mean(Pmask(I(1)-1+min_list(lefti:(righti+1),1)));
    if chooseArs
        [~,~,button(i)]=ginput(1);
    else
        button(i)=1;
    end
end

II = ArthresEdi.CPAP>-1&ArthresEdi.CPAP<1;

ArthresEdi.mean = mean(ArthresEdi.data(button==1));
ArthresEdi.median = median(ArthresEdi.data(button==1));
ArthresEdi.Feupnea = ArthresEdi.median/Veupnea.mean;
ArthresEdi.F_sem = std(ArthresEdi.data(button==1))/length(button==3)^0.5/Veupnea.mean;

%% --save
if saveontherun
    save(filenameanddir,'ArthresEdi','-append');
end

%% Analyze LG from Edi during Varousal
%load(filenameanddir,'VARx1');

%right click to exclude
for i=1:size(Varousal.ranges,1)
    try
        deltaX=4;
        lowerX=Varousal.ranges(i,1)-deltaX;
        upperX=Varousal.ranges(i,2)+deltaX;
        figure(3)
        I=find(Time>lowerX&Time<upperX);
        %     subplot(2,1,1); plot(Time(I),Flow(I),'k',[Varousal.ranges(i,1) Varousal.ranges(i,1)],[-1 1],'r',[Varousal.ranges(i,2) Varousal.ranges(i,2)],[-1 1],'r');
        
        %     lefti=(lowerX-Time(1))/dt+1;
        %     righti=(upperX-Time(1))/dt+1;
        
        xdata.data = Edi(I); %downsample(Pes(I),downsamplefactor);
        xdata.time = Time(I);
        xdata.dt = dt;
        
        modelY = EdiToVflowModel(parameters_x(1,:),xdata);
        
        %     figure(3)
        %I=find(Time>lowerX&Time<upperX);
        ax3(1)=subplot(2,1,1); plot(xdata.time,modelY,'k',[Varousal.ranges(i,1) Varousal.ranges(i,1)],[-1 1],'r',[Varousal.ranges(i,2) Varousal.ranges(i,2)],[-1 1],'r');
        %vol1=cumsum(modelY-mean(modelY))*dt;
        vol1=cumsum(modelY)*dt;
        
        %     vol1=cumsum(Flow(I)-mean(Flow(I)))*dt;
        
        thresh=std(vol1)/10;
        [min_list, max_list] = peakdet(-vol1,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
        %[indVI, indVE] = CheckVind(max_list, min_list);
        %TidalVol = V(indVI)-V(indVE);
        ax3(2)=subplot(2,1,2); plot(Time(I),vol1,'k',[Varousal.ranges(i,1) Varousal.ranges(i,1)],[-1 1],'r',[Varousal.ranges(i,2) Varousal.ranges(i,2)],[-1 1],'r');
        hold('on')
        ax3(2)=subplot(2,1,2); plot(Time(I(1)-1+max_list(:,1)),vol1(max_list(:,1)),'k.');
        ax3(2)=subplot(2,1,2); plot(Time(I(1)-1+min_list(:,1)),vol1(min_list(:,1)),'r.');
        %[ymin,temptemp]=get(ax2(4),'ylim')
        
        lefti=find(Time(I(1)-1+min_list(:,1))<Varousal.ranges(i,1),1,'last');
        righti=find(Time(I(1)-1+min_list(:,1))>Varousal.ranges(i,2),1,'first')-1;
        if isempty(lefti)
            lefti=1;
        end
        if isempty(righti)
            righti=lefti;
        end
        if righti<lefti
            righti=lefti;
        end
        ax3(2)=subplot(2,1,2); plot([Time(I(1)-1+min_list(lefti,1)) Time(I(1)-1+min_list(lefti,1))],[-1 1],'k:',[Time(I(1)-1+min_list(righti+1,1)) Time(I(1)-1+min_list(righti+1,1))],[-1 1],'k:');
        
        xpoly = Time(I(1)-1+min_list(lefti:(righti+1),1));
        ttot = (xpoly(end)-xpoly(1))/(size(xpoly,1)-1);
        ypoly = vol1(min_list(lefti:(righti+1),1));
        ppoly = polyfit(xpoly,ypoly,1);
        Yvoldrift = polyval(ppoly,Time(I));
        
        %repeat find end-exp and end-insp now that we have done some leak correction:
        %     vol3=cumsum(Flow(I)-mean(Flow(I))-ppoly(1))*dt;
        vol3=cumsum(modelY-ppoly(1))*dt;
        
        thresh=std(vol3)/10;
        [min_list, max_list] = peakdet(-vol3,thresh); %in this configuration, max_list(i,1) > minlist(i,1).
        
        ax3(2)=subplot(2,1,2); plot(Time(I(1)-1+max_list(:,1)),vol1(max_list(:,1)),'ko');
        ax3(2)=subplot(2,1,2); plot(Time(I(1)-1+min_list(:,1)),vol1(min_list(:,1)),'ro');
        
        lefti=find(Time(I(1)-1+min_list(:,1))<Varousal.ranges(i,1),1,'last');
        righti=find(Time(I(1)-1+min_list(:,1))>Varousal.ranges(i,2),1,'first')-1;
        if isempty(lefti)
            lefti=1;
        end
        if isempty(righti)||righti==0
            righti=size(min_list,1)-1;
        end
        ax3(2)=subplot(2,1,2); plot([Time(I(1)-1+min_list(lefti,1)) Time(I(1)-1+min_list(lefti,1))],[-1 1],'r:',[Time(I(1)-1+min_list(righti+1,1)) Time(I(1)-1+min_list(righti+1,1))],[-1 1],'r:');
        
        xpoly = Time(I(1)-1+min_list(lefti:(righti+1),1));
        ttot = (xpoly(end)-xpoly(end-1)); %last breath only
        ypoly = vol1(min_list(lefti:(righti+1),1));
        ppoly = polyfit(xpoly,ypoly,1);
        Yvoldrift = polyval(ppoly,Time(I));
        
        Yvol1_exp = ypoly(1:end-1); %end exp
        Yvol1_insp = vol1(max_list(lefti:(righti),1)); %end insp
        Yvoldrift_exp = polyval(ppoly,xpoly); Yvoldrift_exp(end)=[];
        Yvoldrift_insp = polyval(ppoly,Time(I(1)-1+max_list(lefti:(righti),1)));
        
        subplot(2,1,2); plot(Time(I),Yvoldrift);
        hold('off');
        
        endexpvol=Yvol1_exp-Yvoldrift_exp;
        endinspvol=Yvol1_insp-Yvoldrift_insp;
        volbreath=endinspvol-endexpvol;
        LGediVarousal.responsedata(i)=volbreath(end)/ttot*60;%mean(volbreath)/ttot*60; %Edi of last breath only
        LGediVarousal.CPAP(i)=mean(Pmask(I(1)-1+min_list(lefti:(righti+1),1)));
        
        [~,~,button(i)]=ginput(1);
        LGediVarousal.include(i)=button(i)==1;
        linkaxes(ax3,'x');
    catch me
        LGediVarousal.include(i)=0;
    end
end


LGediVarousal.disturbancedata = Varousal.data;

mindist=Veupnea.mean/20;
LGediVarousal.include(LGediVarousal.disturbancedata>(Veupnea.mean-mindist)|LGediVarousal.responsedata<Veupnea.mean)=0;

LGediVarousal.x = Veupnea.mean-LGediVarousal.disturbancedata;
LGediVarousal.y = LGediVarousal.responsedata-Veupnea.mean;
LGediVarousal.value = LGediVarousal.y(LGediVarousal.include)./LGediVarousal.x(LGediVarousal.include);
LGediVarousal.median = median(LGediVarousal.value);
LGediVarousal.arthresmedianF = median(LGediVarousal.responsedata(LGediVarousal.include))/Veupnea.mean;

% LGediVarousal.mean = mean(LGediVarousal.data);
% LGediVarousal.median = median(LGediVarousal.data);
% LGediVarousal.Feupnea = LGediVarousal.median/Veupnea.mean;

%% --save
if saveontherun
    save(filenameanddir,'LGediVarousal','-append');
end

%% --summarydata
if ~exist('LG')
    LG.median=NaN;
end
TraitsSummaryCombined2 = [TraitsSummaryCombined;TraitsCPAPdrops;ArthresEdi.Feupnea;LG.median;Pcrit.PVSlope_VE;LGediVarousal.median;LGediVarousal.arthresmedianF];

%% Tools if needed
if 0
    % Flow analysis
    trange=[Evts.starttimes(I2(m))-60-extra Evts.starttimes(I2(m))+60+extra];
    Flow=Flow.values(Time>trange(1)&Time<=trange(2));
    Time=Time(Time>trange(1)&Time<=trange(2));
    Pmask=Pmask.values(Time>trange(1)&Time<=trange(2));
    epochs1=EpochsXHz(Time>trange(1)&Time<=trange(2));
    EventsAr1=EventsArXHz(Time>trange(1)&Time<=trange(2));
    PCO2=PETCO2(Time>trange(1)&Time<=trange(2));
    minimum_figs_flow=0;
    detrend_flow=0;
    [Bapnea,Te,Ti,BB,BB_t,BB_i_start,BB_i_end,BB_i_mid,VTi,VTe,VE,VI]=Vflowanalysis1(Flow,Time,dt,minimum_figs_flow,detrend_flow);
    
    %PETCO2 analysis
    PETCO2lag_estimated=1;
    minPETCO2fraction=0.4;
    minVTefraction=0.2;
    delay_variable=0;
    Plowerthres=20;
    RemoveErraticHighPICO2=0;
    minimum_figs=0;
    [PETCO2x,PETO2,BB_t_PETCO2fix,PETCO2_Bfix,PICO2_Bfix,PETO2_Bfix,Vflow_to_PETCO2_delay_measured,delaycorr]=PETCO2analysis1(PCO2,PETO2,PETCO2lag_estimated,dt,BB_i_start,Time,Flow,BB_i_mid,analyseO2,VI,VTe,minimum_figs,Plowerthres,RemoveErraticHighPICO2,delay_variable,minPETCO2fraction,minVTefraction,Bapnea,BB_i_end);
end
