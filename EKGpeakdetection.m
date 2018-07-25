function ECG_peak_i = EKGpeakdetection(EKG,Time,dt)

%% Find RR
fs = round(1/dt);
N=length(Time);
%Filter EKG to aid peak identification. 
filter_LFcutoff_butter = 5; %highpass, DC remove
filter_HFcutoff_butter = 20; %lowpass, smooth
filter_order = 2;
[B_butter,A_butter] = butter(filter_order,[filter_LFcutoff_butter filter_HFcutoff_butter]/(fs/2));
ECG_filtered = filtfilt(B_butter,A_butter,EKG);

ECG_flip_on=1;
FlipECG = ((mean(ECG_filtered)-prctile(ECG_filtered,0.5))>(prctile(ECG_filtered,99.5)-mean(ECG_filtered)))
if ECG_flip_on
if FlipECG
    ECG_filtered=-ECG_filtered;
    EKG=-EKG;
    'flipped EKG'
end
end

ECG_filtered_threshold_fraction=0.75;

count=1;
count_min=1;
ECG_filtered_threshold = ECG_filtered_threshold_fraction*[prctile(ECG_filtered,99)-prctile(ECG_filtered,1)]+prctile(ECG_filtered,1);
clear ECG_peak ECG_peak_i ECG_peak_t;
for i=1:N
    if (i>2)&&(i<N)
        if ((ECG_filtered(i)>ECG_filtered(i-1))&&(ECG_filtered(i)>ECG_filtered(i+1))||(ECG_filtered(i)>ECG_filtered(i-2))&&(ECG_filtered(i)>ECG_filtered(i+1))&&(ECG_filtered(i)==ECG_filtered(i-1)))&&(ECG_filtered(i)>ECG_filtered_threshold)
            ECG_peak(count)=ECG_filtered(i);
            ECG_peak_i(count)=i;
            ECG_peak_t(count)=Time(i);
            count=count+1;
        end       
    end
end

%Store data
ECG_peak_step1=ECG_peak;
ECG_peak_i_step1=ECG_peak_i;
ECG_peak_t_step1=ECG_peak_t;

ECG_peak(EKG(ECG_peak_i)<0)=[];
ECG_peak_t(EKG(ECG_peak_i)<0)=[];
ECG_peak_i(EKG(ECG_peak_i)<0)=[];
%% Remove array elements (RAE) 
    %Removes EKG peaks that are within X sec of another
    %The smaller of the two peaks are removed. 

%RAE Setup
ECG_peak_t_threshold_percentofmedian = 33; %specify
    
%RAE X-times
X=3;
for x=1:X
    ECG_peak_t_delta_median = prctile(diff(ECG_peak_t),90);
    ECG_peak_t_threshold_min = ECG_peak_t_delta_median*ECG_peak_t_threshold_percentofmedian/100; 
    i_=0;
    for i=1:length(ECG_peak)
        if ((i-i_)>1)&&((i-i_)<length(ECG_peak))
            if (((ECG_peak_t(i-i_))-(ECG_peak_t(i-i_-1)))<ECG_peak_t_threshold_min)&&((ECG_peak(i-i_))<(ECG_peak(i-i_-1))), %if true, remove element i-i_
                ECG_peak = [ECG_peak(1:(i-i_-1)) ECG_peak((i-i_+1):length(ECG_peak))];
                ECG_peak_i = [ECG_peak_i(1:(i-i_-1)) ECG_peak_i((i-i_+1):length(ECG_peak_i))];
                ECG_peak_t = [ECG_peak_t(1:(i-i_-1)) ECG_peak_t((i-i_+1):length(ECG_peak_t))];
                i_=i_+1;
            elseif (((ECG_peak_t(i-i_+1))-(ECG_peak_t(i-i_)))<ECG_peak_t_threshold_min)&&((ECG_peak(i-i_))<(ECG_peak(i-i_+1))), %if true, remove element i-i_
                ECG_peak = [ECG_peak(1:(i-i_-1)) ECG_peak((i-i_+1):length(ECG_peak))];
                ECG_peak_i = [ECG_peak_i(1:(i-i_-1)) ECG_peak_i((i-i_+1):length(ECG_peak_i))];
                ECG_peak_t = [ECG_peak_t(1:(i-i_-1)) ECG_peak_t((i-i_+1):length(ECG_peak_t))];
                i_=i_+1;                
            end
        end
    end
end

%% Move index to peak of unfiltered EKG
delta=2;
for i=1:length(ECG_peak_i)
    rangei = ECG_peak_i(i)-delta:ECG_peak_i(i)+delta;
    if rangei(1)<0||rangei(end)>length(EKG)
        continue
    end
    [~,temp] = max(EKG(rangei));
    ECG_peak_i(i)=ECG_peak_i(i)+temp-delta-1;
end
%% 
%remove first and last points.
ECG_peak=ECG_peak(2:length(ECG_peak)-1);
ECG_peak_i=ECG_peak_i(2:length(ECG_peak_i)-1);
ECG_peak_t=ECG_peak_t(2:length(ECG_peak_t)-1);

if 0
clear RR RR_t;
for i=1:length(ECG_peak)
    if i<length(ECG_peak)
        RR(i)=ECG_peak_t(i+1)-ECG_peak_t(i);
        RR_t(i)=ECG_peak_t(i+1);
        RR_i(i)=ECG_peak_i(i+1);
    end
end
end

%% Plot
plotfig=1;
if plotfig||Review_data
% Plot Peak detection
%figure('Name',['RR detect'],'color',[1 1 1]);
figure(5)
axx2(1)=subplot(2,1,1); plot(Time,EKG,'k',Time(ECG_peak_i),EKG(ECG_peak_i),'r.');
box('off');
xlabel('Time (s)');
ylabel('EKG');
axx2(2)=subplot(2,1,2); plot(Time,ECG_filtered,ECG_peak_t_step1,ECG_peak_step1,'r.',ECG_peak_t,ECG_peak,'k.',Time,ones(1,length(Time))*ECG_filtered_threshold,'r');
box('off');
xlabel('Time (s)');
ylabel('EKG Filt.');
linkaxes(axx2,'x');
end
%% Quadratic fit to peak to improve RR resolution from 1/fs
if 0
clear ECG_peak_quadfit_t ECG_peak_quadfit
for i=1:length(ECG_peak_i);
[ECG_peak_quadfit_t(i),ECG_peak_quadfit(i)]=PeakFitQuadratic(Time((ECG_peak_i(i)-1):(ECG_peak_i(i)+1)),EKG((ECG_peak_i(i)-1):(ECG_peak_i(i)+1)));
end

clear RR_quadfit RR_quadfit_t;
for i=1:length(ECG_peak_quadfit)
    if i<length(ECG_peak_quadfit)
        RR_quadfit(i)=ECG_peak_quadfit_t(i+1)-ECG_peak_quadfit_t(i);
        RR_quadfit_t(i)=ECG_peak_quadfit_t(i+1);
    end
end

%option:
use_quadfit=1;
if use_quadfit
    RR=RR_quadfit;
    RR_t=RR_quadfit_t;
end
end
%% Use local deviation from the mean >X ms to denote detection error.
if 0
for x=1:3
    deltaRRtotal=zeros(1,length(RR));
    deltaRRtotal2=zeros(1,length(RR));
    deltaRRflag=zeros(1,length(RR));
    for i=1:(length(RR)-1)
        deltaRRtotal(i) = (RR(i+1)+RR(i))/2-RR(i);
    end
    for i=2:(length(RR)-1)
        deltaRRtotal2(i) = abs(deltaRRtotal(i)-deltaRRtotal(i-1)-deltaRRtotal(i+1));
    end  
    %this method may only work for single erroneous beats:
    for i=2:(length(RR)-1)
        if deltaRRtotal2(i)>0.025 %25 ms threshold
            if (deltaRRtotal2(i)>deltaRRtotal2(i-1))&&(deltaRRtotal2(i)>deltaRRtotal2(i+1)) %abs value great than both adjacent values (could also check both are opposite directions?)
                deltaRRflag(i)=1;
                deltaRRflag(i+1)=1;
            end
        end
    end
    figure(1001)
    ax1001(1)=subplot(3,1,1), plot(RR_t,deltaRRtotal,RR_t,deltaRRtotal2);
    ax1001(2)=subplot(3,1,2), plot(RR_t,RR);
    ax1001(3)=subplot(3,1,3), plot(RR_t,deltaRRflag);
    linkaxes(ax1001,'x');
  
    for i=(length(RR)-1):-1:2
        if deltaRRflag(i)
        RR(i)=[];
        RR_i(i)=[];
        RR_t(i)=[]; 
        'single beat error removed';
        end
    end
end

plotfig=0;
if plotfig&&Review_data
figure(3), 
subplot(4,1,1), plot(Time,EKG,Time(ECG_peak_i),EKG(ECG_peak_i),'k.',ECG_peak_quadfit_t,ECG_peak_quadfit,'g.');
xlabel('Time (s)');
ylabel('EKG (mV)');
subplot(4,1,2), plot(Time,ECG_filtered,ECG_peak_t_step1,ECG_peak_step1,'r.',ECG_peak_t,ECG_peak,'k.',Time,ones(1,length(Time))*ECG_filtered_threshold,'r');
xlabel('Time (s)');
ylabel('EKG, peak detect, threshold (mV)');
subplot(4,1,3), plot(RR_t,RR,'k.',RR_quadfit_t,RR_quadfit,'g.');
xlabel('Time (s)');
ylabel('RR (1/s)');
end
end

%% Plots EKG Analysis
if 0
plotfig=1;
if plotfig||Review_data
% Plot Peak detection
figure('Name',['RR detect'],'color',[1 1 1]);
axx2(1)=subplot(3,1,1), plot(Time,EKG,ECG_peak_quadfit_t,ECG_peak_quadfit,'k.');
box('off');
xlabel('Time (s)');
ylabel('EKG');
axx2(2)=subplot(3,1,2), plot(Time,ECG_filtered,ECG_peak_t_step1,ECG_peak_step1,'r.',ECG_peak_t,ECG_peak,'k.',Time,ones(1,length(Time))*ECG_filtered_threshold,'r');
box('off');
xlabel('Time (s)');
ylabel('EKG Filt.');
axx2(3)=subplot(3,1,3), plot(RR_t,1000*RR,'.');
box('off');
xlabel('Time (s)');
ylabel('RR (ms)');
linkaxes(axx2,'x');
end

% Plot RR
plotfig=0;
if plotfig
figure('Name',[filename,' RR'],'color',[1 1 1]);
subplot(2,1,2), plot(RR_t,1000*RR,'.',time_rs,1000*RR_rs);
box('off');
xlabel('Time (s)');
ylabel('RR, resampled (ms)');
end
end
%% Mean
