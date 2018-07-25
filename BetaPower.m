function [WakeSleep,WakeSleepInfo]=BetaPower(WakeSleepInfo,EpochsXHz,EventsArXHz,Time,SaO2,Flow,plotfigs)
%% Beta Power wake versus sleep detection
%Options
filterEEGsignal=0;
%EEG_O2_A1,EEG_F3_A2,EEG_C4_A1,EEG_C3_A2
signallist=WakeSleepInfo.EEGoptions;
suf = '_clean';
signalnothere=[];
Nsignals=length(signallist);
n=1;
for i=1:Nsignals
    try
        signal{n}=evalin('caller',[signallist{i} suf]); %isempty(eval(signallist{i}))
        n=n+1;
    catch me
        disp(me.message);
        signallist{i}=[];
    end
end

count=0;
for m=1:length(WakeSleepInfo.EEGoptions)
    try
        EEGselect1=signal{m}; %EEG_C3_A2; EEG_C4_A1; EEG_O2_A1; EEG_F3_A2;
        if filterEEGsignal
            EEGselect1filt = filterEEG(EEGselect1);
            [PbetaLog,PbetaLogFilt,Pbetat] = PowerEEG(WakeSleepInfo,Time,EEGselect1filt);
        else
            [PbetaLog,PbetaLogFilt,Pbetat] = PowerEEG(WakeSleepInfo,Time,EEGselect1);    
        end      
        if count==0
            [WakeSleepScored] = WakeOrSleep(EpochsXHz,EventsArXHz,Time,WakeSleepInfo,Pbetat,SaO2); 
            WakeSleepInfo.Fwake = nanmean(WakeSleepScored);
            count=1;
        end
        
        [PbetaPrWake,PbetaPrWakeFilt,TempX,PrWake,PrWakeSEM,ydatamodel,Xmodel,WakeSleepInfo.AUC,WakeSleepInfo.OptimalThresholdYoudenJ,WakeSleepInfo.PrOptimalThresholdYoudenJ] = ProbabilityOfWake(PbetaLog,PbetaLogFilt,WakeSleepScored,WakeSleepInfo);
            WakeSleepInfo.logpower50=Xmodel(1); 
            WakeSleepInfo.sigmoidcoefficient=Xmodel(2);
        figure(31);
        subplot(1,1,1); plot(TempX,PrWake,'k',TempX,PrWake+1.96*PrWakeSEM,'k:',TempX,PrWake-1.96*PrWakeSEM,'k:',WakeSleepInfo.logpower50,0.5,'.',WakeSleepInfo.OptimalThresholdYoudenJ,WakeSleepInfo.PrOptimalThresholdYoudenJ,'r.');
        %subplot(1,1,1); plot(TempX,PrWake,'k',TempXp25,PrWake,'k:',TempXp75,PrWake,'k:',WakeSleepInfo.logpower50,0.5,'.');
        box('off');
        xlabel('Power');
        ylabel('Probability of Wake');
        WakeSleepInfo.AUC_M(m)=WakeSleepInfo.AUC;
    catch me
        disp(me);
        WakeSleepInfo.AUC_M(m)=NaN;
    end
end

        
[~,bestm] = max(WakeSleepInfo.AUC_M);

for m=bestm
    try
        EEGselect1=signal{m}; %EEG_C3_A2; EEG_C4_A1; EEG_O2_A1; EEG_F3_A2;
        if filterEEGsignal
            EEGselect1filt = filterEEG(EEGselect1);
            [PbetaLog,PbetaLogFilt,Pbetat] = PowerEEG(WakeSleepInfo,Time,EEGselect1filt);
        else
            [PbetaLog,PbetaLogFilt,Pbetat] = PowerEEG(WakeSleepInfo,Time,EEGselect1);    
        end
        %[WakeSleepScored] = WakeOrSleep(EpochsXHz,EventsArXHz,Time,WakeSleepInfo,Pbetat,SaO2);
        [PbetaPrWake,PbetaPrWakeFilt,TempX,PrWake,PrWakeSEM,ydatamodel,Xmodel,WakeSleepInfo.AUC,WakeSleepInfo.OptimalThresholdYoudenJ,WakeSleepInfo.PrOptimalThresholdYoudenJ] = ProbabilityOfWake(PbetaLog,PbetaLogFilt,WakeSleepScored,WakeSleepInfo);
            WakeSleepInfo.logpower50=Xmodel(1); 
            WakeSleepInfo.sigmoidcoefficient=Xmodel(2);
        figure(31);
        subplot(1,1,1); plot(TempX,PrWake,'k',TempX,PrWake+1.96*PrWakeSEM,'k:',TempX,PrWake-1.96*PrWakeSEM,'k:',WakeSleepInfo.logpower50,0.5,'.',WakeSleepInfo.OptimalThresholdYoudenJ,WakeSleepInfo.PrOptimalThresholdYoudenJ,'r.');
        %subplot(1,1,1); plot(TempX,PrWake,'k',TempXp25,PrWake,'k:',TempXp75,PrWake,'k:',WakeSleepInfo.logpower50,0.5,'.');
        WakeSleepInfo.AUC_M(m)=WakeSleepInfo.AUC;
    catch me
        WakeSleepInfo.AUC_M(m)=NaN;
    end
end

figure(33);
TempX_transformed = 1./(1+exp(-WakeSleepInfo.sigmoidcoefficient*(TempX-WakeSleepInfo.logpower50)));
subplot(1,1,1); plot(TempX_transformed,PrWake);
xlabel('TempX_transformed');
ylabel('Probability of Wake');
WakeSleep = interp1(Pbetat,PbetaPrWakeFilt,Time,'linear');


if plotfigs
figure(30); set(gcf,'color',[1 1 1]);
ax30(1)=subplot(4,1,1); plot(Time,EventsArXHz,Time,EpochsXHz,Pbetat,isnan(WakeSleepScored)); ylim([-0.1 5.1]);
box('off');
ax30(2)=subplot(4,1,2);
plot(Pbetat,PbetaPrWake,'color',[0.8 0.8 0.8]); hold('on');
plot(Pbetat,PbetaPrWakeFilt,'k'); hold('off'); box('off'); ylim([-0.1 1.1]); ylabel('Prwake');


ax30(3)=subplot(4,1,3); plot(Time,EEGselect1,'k');  ylim([-50 50]); box('off');

ax30(4)=subplot(4,1,4);
plot(Time,Flow,'k'); ylim([-1 1]); box('off'); hold('on');

tempT=Time; tempT(WakeSleep<-0.5|WakeSleep>=0)=NaN;
tempV=Flow; tempV(WakeSleep<-0.5|WakeSleep>=0)=NaN;
plot(tempT,tempV,'color',[0.33 0 0]); ylim([-1 1]); hold('on');
clear tempT tempV

tempT=Time; tempT(WakeSleep<0|WakeSleep>=0.5)=NaN;
tempV=Flow; tempV(WakeSleep<0|WakeSleep>=0.5)=NaN;
plot(tempT,tempV,'color',[0.67 0 0]); ylim([-1 1]); hold('on');
clear tempT tempV

tempT=Time; tempT(WakeSleep<0.5)=NaN;
tempV=Flow; tempV(WakeSleep<0.5)=NaN;
plot(tempT,tempV,'color',[1 0 0]); ylim([-1 1]); hold('off');
clear tempT tempV

linkaxes(ax30,'x');
end
%channels{length(channels)+1}='WakeSleep';



%% WakeOrSleep
function [WakeSleep] = WakeOrSleep(EpochsXHz,EventsArXHz,Time,WakeSleepInfo,Pbetat,SaO2)
dt = Time(2)-Time(1);

%Use wake only within 1 min of sleep, also exclude times when SpO2 is off
downsamplefactor = round(1/dt);
Temp2 = downsample(EpochsXHz,downsamplefactor);
Temp2_time = downsample(Time,downsamplefactor);
%figure(32); stairs(Temp2_time,Temp2);
nearestneighbor = round(WakeSleepInfo.nearbyduration*downsamplefactor*dt); %1 min
Temp3=0*Temp2; %Are we within "nearestneigbour" epochs of sleep?
for i=(1+nearestneighbor):(length(Temp2_time)-nearestneighbor);
    I = (i-nearestneighbor):(i+nearestneighbor);
    Temp3(i) = max(Temp2(I)>0&Temp2(I)<=5); %test if scored sleep is in range +/- neighrestneighbor
end
% figure(32); stairs(Temp2_time,[Temp2,Temp3]);
NearSleep = interp1(Temp2_time,Temp3,Pbetat,'nearest');
% figure(32); plot(Time,EpochsXHz,Pbetat,NearSleep);
TempSpO2 = interp1(Time,SaO2,Pbetat,'nearest');
SpO2On = TempSpO2>10;

%figure(29); plot(Time,EpochsXHz,Pbetat,[NearSleep SpO2On]);

WakeTemp = EventsArXHz;
if ~WakeSleepInfo.scoredarousalsinwake
    WakeTemp(EpochsXHz==0)=1;
end
WakeTemp(EpochsXHz==8)=-1;

WakeSleep = interp1(Time,WakeTemp,Pbetat,'nearest'); %interp1 doesn't handle NaN

WakeSleep(WakeSleep==8)=NaN; %unknown sleep is NaN

if WakeSleepInfo.useonlySpO2on==1
    WakeSleep(SpO2On==0)=NaN;    
elseif WakeSleepInfo.useonlySpO2on==0
    WakeSleep(NearSleep==0|SpO2On==0)=NaN;
else
    %WakeSleep(NearSleep==0|SpO2On==0)=NaN;
end

%% Probability and Statistics
function [PbetaPrWake,PbetaPrWakeFilt,TempX,PrWake,PrWakeSEM,ydatamodel,Xmodel,AUC,OptimalThresholdYoudenJ,PrOptimalThresholdYoudenJ] = ProbabilityOfWake(PbetaLog,PbetaLogFilt,WakeSleepScored,WakeSleepInfo)
%WakeSleepInfo.Fwake = nanmean(WakeSleepScored);
gamma = WakeSleepInfo.Fwake/(1-WakeSleepInfo.Fwake);
TempX_ = prctile(PbetaLogFilt(~isnan(WakeSleepScored)),[0:10:90 99.5]); %powerlevels

PrWake = zeros(1,length(TempX_)-1);
PrWakeSEM = zeros(1,length(TempX_)-1);
PrWake_N = 0*PrWake;
TempX = 0*PrWake;
%TempXp75 = 0*PrWake;
%TempXp25 = 0*PrWake;
for i=1:(length(TempX_)-1)
    %     if i==length(TempX)
    %         data = WakeSleepScored(~isnan(WakeSleepScored)&PbetaLogFilt>=TempX(i));
    %     else
    data = WakeSleepScored(~isnan(WakeSleepScored)&PbetaLogFilt>=TempX_(i)&PbetaLogFilt<TempX_(i+1));
    
    %     end
    PrWake(i) = nanmean(data);
    
    %PrWake(i)=1/(1+Y)
    %Y=1/PrWake(i)-1; %Y=Nsleep/Nwake
    PrWake(i)=1/(1+(1/PrWake(i)-1)*gamma); %PrWake2(i)=1/(1+Y*gamma);
    PrWake_N(i) = length(~isnan(data));
    PrWakeSEM(i) = nanstd(data)/(((1-WakeSleepInfo.Foverlap)*PrWake_N(i))^.5);
    TempX(i) = nanmedian(PbetaLogFilt(~isnan(WakeSleepScored)&PbetaLogFilt>=TempX_(i)&PbetaLogFilt<TempX_(i+1)));
    %TempXp75(i) = prctile(PbetaLogFilt(~isnan(WakeSleepScored)&PbetaLogFilt>=TempX_(i)&PbetaLogFilt<TempX_(i+1)),75);
    %TempXp25(i) = prctile(PbetaLogFilt(~isnan(WakeSleepScored)&PbetaLogFilt>=TempX_(i)&PbetaLogFilt<TempX_(i+1)),25);
    
end

if mean(diff(PrWake))>0
    I = find(PrWake>0.5,1);
else
    I = find(PrWake<0.5,1);  
end
I = (I-1):(I+1);
I(I>length(TempX))=[];
I(I<1)=[];

if ~isempty(I)
    logpower50 = interp1(PrWake(I),TempX(I),0.5,'linear');
else
    logpower50 = TempX(end);
end

%formally fit a sigmoid to probability data
xdata=TempX; %logpower
ydata=PrWake;
lsqoptions=optimset('display','off');
if mean(diff(PrWake))>0
x0 = [logpower50 7];
upper = [logpower50+1 20];
lower = [logpower50-1 0.2];
fun = @(x,xdata)1./(1+exp(-x(2)*(xdata-x(1))));
else
x0 = [logpower50 -7];
upper = [logpower50+1 -20];
lower = [logpower50-1 -0.2];
fun = @(x,xdata)1./(1+exp(-x(2)*(xdata-x(1))));    
end
[Xmodel]=lsqcurvefit(fun,x0,xdata,ydata,lower,upper,lsqoptions);
sigmoidcoefficient_model = Xmodel(2);
logpower50_model = Xmodel(1);
ydatamodel = fun(Xmodel,xdata);

% Median Filter PAlpha+PBeta and plot (~15 s)
%medianfiltertime = 6;
    PbetaPrWake = fun(Xmodel,PbetaLog);
    PbetaPrWakeFilt = fun(Xmodel,PbetaLogFilt);

%Sensitivity
TempX2 = prctile(PbetaLogFilt(~isnan(WakeSleepScored)),[0.1 1:99 99.9]); %powerlevels
ConditionPositive = sum((WakeSleepScored(~isnan(WakeSleepScored))==1));
ConditionNegative = sum((WakeSleepScored(~isnan(WakeSleepScored))==0));
clear TestPositives TestNegatives TP TN FP PrWake2
if mean(diff(PrWake))>0
for i=1:length(TempX2)
    PrWake2(i)=nanmean(WakeSleepScored(PbetaLogFilt>TempX2(i)));
    TestPositives(i) = sum(PbetaLogFilt(~isnan(WakeSleepScored))>TempX2(i));
    TestNegatives(i) = sum(PbetaLogFilt(~isnan(WakeSleepScored))<=TempX2(i));
    TP(i) = sum((WakeSleepScored(~isnan(WakeSleepScored))==1)&(PbetaLogFilt(~isnan(WakeSleepScored))>TempX2(i))); %True Positives
    TN(i) = sum((WakeSleepScored(~isnan(WakeSleepScored))==0)&(PbetaLogFilt(~isnan(WakeSleepScored))<=TempX2(i))); %True Negatives
    FP(i) = sum((WakeSleepScored(~isnan(WakeSleepScored))==0)&(PbetaLogFilt(~isnan(WakeSleepScored))>TempX2(i))); %False Positives
end
else
for i=1:length(TempX2)
    PrWake2(i)=nanmean(WakeSleepScored(PbetaLogFilt>TempX2(i)));
    TestPositives(i) = sum(PbetaLogFilt(~isnan(WakeSleepScored))<=TempX2(i));
    TestNegatives(i) = sum(PbetaLogFilt(~isnan(WakeSleepScored))>TempX2(i));
    TP(i) = sum((WakeSleepScored(~isnan(WakeSleepScored))==1)&(PbetaLogFilt(~isnan(WakeSleepScored))<=TempX2(i))); %True Positives
    TN(i) = sum((WakeSleepScored(~isnan(WakeSleepScored))==0)&(PbetaLogFilt(~isnan(WakeSleepScored))>TempX2(i))); %True Negatives
    FP(i) = sum((WakeSleepScored(~isnan(WakeSleepScored))==0)&(PbetaLogFilt(~isnan(WakeSleepScored))<=TempX2(i))); %False Positives
end    
end
% 
% figure();
% plot(TempX2,PrWake2); ylabel('PrWake for Beta>X');

I=find(TempX2>logpower50_model,1);
WakeSleepInfo.OptimalThreshold = TempX2(I);
TPR = TP/ConditionPositive; %Sens
FPR = FP/ConditionNegative; %1-Spec
[YoudenJstatisticBest,YoudenJstatisticAlli] = max(TPR-FPR);
OptimalThresholdYoudenJ = TempX2(YoudenJstatisticAlli);
PrOptimalThresholdYoudenJ = fun(Xmodel,OptimalThresholdYoudenJ);
PPV = TP./TestPositives;
NPV = TN./TestNegatives;
Accuracy = (TP+TN)/(ConditionPositive+ConditionNegative);

TPatOP=TPR(I);
FPatOP=FPR(I);
AccuracyatOP=Accuracy(I);
TPatYoudenJ=TPR(YoudenJstatisticAlli);
FPatYoudenJ=FPR(YoudenJstatisticAlli);
AccuracyatYoudenJ=Accuracy(YoudenJstatisticAlli);

AUC = abs(trapz(FPR,TPR)); %abs is in case data are presented right to left.

%formally fit another sigmoid to probability data, forcing the optimal betapower to be the midpoint (50%).

    xdata=TempX; %logpower
    ydata=PrWake;
    lsqoptions=optimset('display','off');
    x0 = [7];
    upper = [20];
    lower = [0.2];
    fun = @(x,xdata)1./(1+exp(-x(1)*(xdata-OptimalThresholdYoudenJ)));
    
    [Xmodel2]=lsqcurvefit(fun,x0,xdata,ydata,lower,upper,lsqoptions);
if 0
    Xmodel=[OptimalThresholdYoudenJ Xmodel2(1)];
    ydatamodel = fun(Xmodel,xdata);

    PbetaPrWake = fun(Xmodel,PbetaLog);
    PbetaPrWakeFilt = fun(Xmodel,PbetaLogFilt);
end

%% filter EEG, Fpass=16Hz, Fstop=14Hz using fdatool
function y = filterEEG(x)

Hd = dsp.BiquadFilter( ...
        'Structure', 'Direct form II', ...
        'SOSMatrix', [1 -2 1 1 -1.39946994856687 0.928188937728236; 1 -2 1 1 ...
        -1.3062188506105 0.799707625432376; 1 -2 1 1 -1.22596816093582 ...
        0.689138268631091; 1 -2 1 1 -1.15699384015505 0.594105486788932; 1 -2 1 ...
        1 -1.09784030335654 0.512603775802443; 1 -2 1 1 -1.047287937652 ...
        0.442952753694169; 1 -2 1 1 -1.00432139291596 0.383753471610781; 1 -2 1 ...
        1 -0.968100845393655 0.333848820837472; 1 -2 1 1 -0.937937078727909 ...
        0.292289199450341; 1 -2 1 1 -0.91327051967812 0.258303627741331; 1 -2 1 ...
        1 -0.893654027815411 0.231276035869533; 1 -2 1 1 -0.878739105316625 ...
        0.210726263722812; 1 -2 1 1 -0.8682651697045 0.196295280905025; 1 -2 1 1 ...
        -0.862051565452499 0.187734180329483; 1 -1 0 1 -0.429996025243795 0], ...
        'ScaleValues', [0.831914721573777; 0.776481619010718; ...
        0.728776607391729; 0.687774831735995; 0.652611019789747; ...
        0.622560172836542; 0.597018716131686; 0.575487416557782; ...
        0.557556569544562; 0.542893536854863; 0.531232515921236; ...
        0.522366342259859; 0.516140112652381; 0.512446436445496; ...
        0.714998012621898; 1]);

y = step(Hd,x);
