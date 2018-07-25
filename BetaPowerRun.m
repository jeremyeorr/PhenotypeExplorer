function [WakeSleep,WakeSleepInfo,EEGsignals,EEGsignalsOut,Powersignals,PowersignalsOut] = BetaPowerRun(EEGsignals,dt,BetaPowerSettings)

Time = evalin('caller',BetaPowerSettings.timestr);
%EEGsignals = {'EEG1','EEG2','EEG3','EEG4'};
Nsignals=length(EEGsignals);
n=1;

for i=1:Nsignals
    try
        signal{n}=evalin('caller',[EEGsignals{i} BetaPowerSettings.suf1]); %isempty(eval(signallist{i}))
        eval([EEGsignals{i} '=signal{n};']);
        n=n+1;
    catch me
        disp(me.message);
        EEGsignals{i}=[];
    end
end

processEEG=1;
if processEEG
    try
        EKG=evalin('caller',['EKG' BetaPowerSettings.suf1]);
    catch me
    end
end


if exist('EKG','var')&&processEEG
    % Create new "clean" EEG signals
    lefttemplatetime = 0.05; righttemplatetime = 0.05;
    contaminationmagnitudethreshold = 4;
    leftdi=round(1/dt*lefttemplatetime); rightdi=round(1/dt*righttemplatetime);
    clear temp contaminationmagnitude contaminationmagnitudepost
   
    
    if processEEG
        for i=1:length(EEGsignals)
            eval([EEGsignals{i} '_clean=' EEGsignals{i} ';'])
        end
        
        % EKG artifact removal
        ECG_peak_i = EKGpeakdetection(EKG,Time,dt);
        
        for i=1:length(EEGsignals)
            [~,template] = crosscontaminated(eval([EEGsignals{i}]),ECG_peak_i,leftdi,rightdi,1,BetaPowerSettings.polyorder,1);
            contaminationmagnitude(i)=max(abs(template));
        end
        disp(['Decontaminating ' num2str(sum(contaminationmagnitude>contaminationmagnitudethreshold)) ' EEG signals']);
        for i=1:length(EEGsignals)
            if contaminationmagnitude(i)>contaminationmagnitudethreshold
                eval([EEGsignals{i} '_clean=crosscontaminated(' EEGsignals{i} ',ECG_peak_i,leftdi,rightdi,0,BetaPowerSettings.polyorder,1);']);
            end
        end
    end
        % after
        for i=1:length(EEGsignals)
            [~,template] = crosscontaminated(eval([EEGsignals{i} '_clean']),ECG_peak_i,leftdi,rightdi,1,BetaPowerSettings.polyorder,1);
            contaminationmagnitudepost(i)=max(abs(template));
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

Flow = evalin('caller',[BetaPowerSettings.flowstr BetaPowerSettings.suf1]);
SaO2 = evalin('caller',['SaO2' BetaPowerSettings.suf1]);
EventsArXHz = evalin('caller',[BetaPowerSettings.arstr BetaPowerSettings.suf1]);
EpochsXHz = evalin('caller',[BetaPowerSettings.epochsstr BetaPowerSettings.suf1]);

[WakeSleep,WakeSleepInfo]=BetaPower(WakeSleepInfo,EpochsXHz,EventsArXHz,Time,SaO2,Flow,plotfigs);
[~,bestEEG]=max(WakeSleepInfo.AUC_M);
bestEEGstr = [WakeSleepInfo.EEGoptions{bestEEG} '_clean'];

for i=1:length(EEGsignals)
    EEGsignalsOut{i}=eval([EEGsignals{i} '_clean;']);
end

%% Power in other bands
temp = WakeSleepInfo;
temp.beta=[1 4];
[~,Pdeltalogfilt,Pbetat] = PowerEEG(temp,Time,eval(bestEEGstr));
Pdeltalogfilt = interp1(Pbetat,Pdeltalogfilt,Time,'nearest'); %interp1 doesn't handle NaN

temp.beta=[4 7];
[~,Pthetalogfilt,~] = PowerEEG(temp,Time,eval(bestEEGstr));
Pthetalogfilt = interp1(Pbetat,Pthetalogfilt,Time,'nearest'); %interp1 doesn't handle NaN

temp.beta=[8 12];
[~,Palphalogfilt,~] = PowerEEG(temp,Time,eval(bestEEGstr));
Palphalogfilt = interp1(Pbetat,Palphalogfilt,Time,'nearest'); %interp1 doesn't handle NaN

temp.beta=[12 16];
[~,Psigmalogfilt,~] = PowerEEG(temp,Time,eval(bestEEGstr));
Psigmalogfilt = interp1(Pbetat,Psigmalogfilt,Time,'nearest'); %interp1 doesn't handle NaN

temp.beta=[16 32];
[~,Pbetalogfilt,~] = PowerEEG(temp,Time,eval(bestEEGstr));
Pbetalogfilt = interp1(Pbetat,Pbetalogfilt,Time,'nearest'); %interp1 doesn't handle NaN

Powersignals = {'Pdeltalogfilt','Pthetalogfilt','Palphalogfilt','Psigmalogfilt','Pbetalogfilt'};

for i=1:length(Powersignals)
    PowersignalsOut{i}=eval(Powersignals{i});
end

if 1
figure()
power=[Pdeltalogfilt Pthetalogfilt Palphalogfilt Psigmalogfilt Pbetalogfilt EventsArXHz WakeSleep];
plot(Time,power)
legend('Pdeltalogfilt', 'Pthetalogfilt', 'Palphalogfilt', 'Psigmalogfilt', 'Pbetalogfilt', 'EventsArXHz', 'WakeSleep')
end

