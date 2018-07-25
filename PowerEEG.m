%% Power EEG in one or two bands
function [PbetaLog,PbetaLogFilt,Pbetat] = PowerEEG(WakeSleepInfo,Time,EEGselect1)
dt=Time(2)-Time(1);
w = hann(WakeSleepInfo.fft_length); w = w/rms(w);

T0 = WakeSleepInfo.fft_length*dt;
df = 1/T0;

betai = [ceil(WakeSleepInfo.beta(1)/df+1) ceil(WakeSleepInfo.beta(2)/df+1)];
if betai(1)==1
    betai(1)=2;
end
if betai(2)>floor(WakeSleepInfo.fft_length/2);
    betai(2)=floor(WakeSleepInfo.fft_length/2);
end
%beta_actual = [(betai(1)-1)*df-df/2 (betai(2)-1)*df+df/2];
%beta_actual_df = (betai(2)-betai(1)+1)*df;
if isfield(WakeSleepInfo,'alpha')
    alphai = [round(WakeSleepInfo.alpha(1)/df+1) floor(WakeSleepInfo.alpha(2)/df+1)];
end

%Setup window width
N=length(Time);
Nstep = round(WakeSleepInfo.fft_length*(1-WakeSleepInfo.Foverlap));
% Nwindows = floor((N/WakeSleepInfo.fft_length - 1)/(1-WakeSleepInfo.Foverlap) + 1);
Nwindows = floor((N - WakeSleepInfo.fft_length)/Nstep);
%Initialize
Pbeta = zeros(Nwindows,1);
Pbetat = zeros(Nwindows,1);
if isfield(WakeSleepInfo,'alpha')
    compIndex = [alphai(1):alphai(2) betai(1):betai(2)]; %bins to analyze power in
else
    compIndex = [betai(1):betai(2)]; %bins to analyze power in
end
for i=length(compIndex):-1:1 %code to remove duplicates in case beta band and alpha band overlap
    if sum(compIndex(i)==compIndex)>1
        compIndex(i)=[];
    end
end
F_compIndex = (compIndex-1)*df; %Frequency at each bin
for i=1:Nwindows
    I=(1:WakeSleepInfo.fft_length)+(i-1)*Nstep;
    X = fft((EEGselect1(I)))/WakeSleepInfo.fft_length*2; %faster without windowing
    %X = fft((EEGselect1(I)).*w)/WakeSleepInfo.fft_length*2; %faster without detrend and windowing
    %X = fft((EEGselect1(I)-mean(EEGselect1(I))).*w)/WakeSleepInfo.fft_length*2;
    %X = fft((detrend(EEGselect1(I))).*w)/WakeSleepInfo.fft_length*2;
    Pbeta(i,:) = sum(conj(X(compIndex)).*X(compIndex).*(F_compIndex'.^WakeSleepInfo.Fpower))*df;

    if 0
        figure(32)
        plot(0:df:(1/dt-df),abs(conj(X).*X)); xlim([0 24]);
    end
    
    Pbetat(i) = Time(I(end));
end

PbetaLog = log10(Pbeta);
Nmedianfilt = round(WakeSleepInfo.medianfiltertime/(WakeSleepInfo.fft_length*dt)/(1-WakeSleepInfo.Foverlap));
PbetaLogFilt = medfilt1(PbetaLog,Nmedianfilt); %median filter using data from +/- Nmedianfilt/2 elements