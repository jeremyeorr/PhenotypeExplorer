function [I,Vdot,VT,Ti,Te]=Vflowanalysis2(Vflow,time,dt,minimum_figs)
warning('off');
%% Respiratory Trace Analyse Using Volume as the primary method of breath ID
%[Pxy,F] = CPSD(X,Y,WINDOW,NOVERLAP,NFFT,Fs)
% Respiratory Trace Analyse
Nloops=1;

%transpose if needed:
if size(Vflow,2)<size(Vflow,1)
    Vflow=Vflow';
end
%transpose if needed:
if size(time,2)<size(time,1)
    time=time';
end

trendx=0*time;

vol1=cumsum(Vflow)*dt;
leak = mean(Vflow);
voldetrend1 = leak*(time-time(1));
vol2 = vol1-voldetrend1;
   
if ~minimum_figs
    figure(); 
    ax1(1)=subplot(Nloops+1,1,1); plot(time,vol1-voldetrend1);
    ax1(2)=subplot(Nloops+1,1,2); plot(time,vol2);
    linkaxes(ax1,'x');
end

if 1 %low pass filter original for apnea detection later
filter_HFcutoff_butter0 = 1;
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
Vflow = filtfilt(B_butter0,A_butter0,Vflow);
end

% figure(100);
% plot(vol2-mean(vol2))

%Setup gentle filter for use in analysis
if 0
filter_HFcutoff_butter0 = 6;
filter_order0 = 2;
[B_butterHcut,A_butterHcut] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
else
filter_HFcutoff_butter1 = 2;
filter_LFcutoff_butter1 = 1/60;
filter_order1 = 2;
[B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dt/2));    
end


if 0 %code to find first and last zero crossings...
    I1=find(detrend(Vflow)<0,1); I2=find(detrend(Vflow)>0);
    Index0 = I2(find(I2>I1,1));
    temp = detrend(Vflow);
    if size(Vflow,1)>size(Vflow,2)
        temp=temp';
    end
    temp=fliplr(temp);
    I1=find(temp>0,1); I2=find(temp<0);
    IndexN = length(Vflow)-I2(find(I2>I1,1));
end


filtervol=1;
if filtervol
vol2_filtered1 = filtfilt(B_butterHcut,A_butterHcut,detrend(vol2));
%vol2_filtered1 = filter(B_butterHcut,A_butterHcut,vol2);
Vflow_filtered1 = [diff(vol2_filtered1)/dt]; Vflow_filtered1=[Vflow_filtered1(1) Vflow_filtered1];
else
Vflow_filtered1 = filtfilt(B_butterHcut,A_butterHcut,Vflow-mean(Vflow));
vol2_filtered1 = cumsum(Vflow_filtered1)*dt;
end
voldetrend2 = vol2_filtered1 - vol2;


if ~minimum_figs
    ax1(2)=subplot(Nloops+1,1,2); plot(time,vol2_filtered1);
    linkaxes(ax1,'x');
end

%% loop -- estimate start insp values, then detrend vol, then re-estimate
%start insp times...

for XX=1:Nloops
%% detecting peaks

% accompanying filtered flow trace
%Vflow_filtered1 = [diff(vol2_filtered1)/dt]; Vflow_filtered1=[Vflow_filtered1(1) Vflow_filtered1];

[min_list,max_list] = peakdet(-vol2_filtered1,0.001*std(vol2_filtered1)); %in this configuration, max_list(i,1) > minlist(i,1).
BB_i_start=min_list(:,1);
BB_i_mid=max_list(:,1);
if BB_i_start(1)==1;
    BB_i_start(1)=[];
end
if BB_i_start(end)==length(vol2_filtered1);
    BB_i_start(end)=[];
end
if BB_i_start(1)>BB_i_mid(1)
    BB_i_mid(1)=[];
end
while BB_i_mid(end)>BB_i_start(end)
    BB_i_mid(end)=[];
end

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.-',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.-');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
end

%estimate normal tidal volume from the distance between upper and lower
%vols.

vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest
    
    
if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow,time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
end
volupperminusvollower = volupper-vollower;
VTmean = nanmean(volupper-vollower);

BB_i_end = BB_i_start(2:end);
BB_i_start = BB_i_start(1:end-1);

VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;
criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
Vdot = VT./Ttot;
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


%% Remove extreme mini-breaths 

VTi_thres = VTmean*0.04;
VTe_thres = VTmean*0.04;
Ti_thres = normalTtot/10;
Te_thres = normalTtot/10;

%If BOTH inspiration and expiration are either short or small, remove the breath
criteria=find(((VTi<VTi_thres)|(Ti<Ti_thres))&((VTe<VTe_thres)|(Te<Te_thres)));

if ~isempty(criteria)
firstbreathisbad=criteria(1)==1;
while firstbreathisbad
   criteria(1)=[];
   BB_i_start(1)=[];
   BB_i_mid(1)=[];
   BB_i_end(1)=[];   
   criteria=criteria-1;   
   if isempty(criteria)
       break
   end
   firstbreathisbad=criteria(1)==1;
end
   criteria2=criteria-1;
   BB_i_start(criteria)=[];
   BB_i_mid(criteria)=[];
   BB_i_end(criteria2)=[]; 
end

%% recalculate breath info:

for i=1:length(BB_i_start)
    [~,tempi]=max(vol2_filtered1(BB_i_start(i):BB_i_end(i)));
    BB_i_mid(i,1)=BB_i_start(i)+tempi-1;
end

vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest
  
if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;

criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
Vdot = VT./Ttot;
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));

%% Remove less extreme mini-breaths 

VTi_thres2 = VTmean*0.08;
Ttot_thres = normalTtot*0.5;

%If BOTH inspiration and expiration are either short or small, remove the breath
criteria=find(((VTi<VTi_thres2)|(Ttot<Ttot_thres)));

if ~isempty(criteria)
firstbreathisbad=criteria(1)==1;
while firstbreathisbad
   criteria(1)=[];
   BB_i_start(1)=[];
   BB_i_mid(1)=[];
   BB_i_end(1)=[];   
   criteria=criteria-1;   
   if isempty(criteria)
       break
   end
   firstbreathisbad=criteria(1)==1;
end
   criteria2=criteria-1;
   BB_i_start(criteria)=[];
   BB_i_mid(criteria)=[];
   BB_i_end(criteria2)=[]; 
end

%% recalculate breath info:

for i=1:length(BB_i_start)
    [~,tempi]=max(vol2_filtered1(BB_i_start(i):BB_i_end(i)));
    BB_i_mid(i,1)=BB_i_start(i)+tempi-1;
end

vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest
  
if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;

criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
Vdot = VT./Ttot;
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));

%% If we find a  short or small inspiration, presume it is part of the previous expiration:
criteria=find((Ti<Ti_thres)|(VTi<VTi_thres*0.67)|((Ti<Ti_thres*2)&(VTi<VTi_thres*2)));
%firstbreathisbad=sum(criteria==1);
if ~isempty(criteria)
firstbreathisbad=criteria(1)==1;
while firstbreathisbad
   criteria(1)=[];
   BB_i_start(1)=[];
   BB_i_mid(1)=[];
   BB_i_end(1)=[];   
   criteria=criteria-1;   
   if isempty(criteria)
       break
   end
   firstbreathisbad=criteria(1)==1;
end
   criteria2=criteria-1;
   BB_i_start(criteria)=[];
   BB_i_mid(criteria)=[];
   BB_i_end(criteria2)=[]; 
end

%% recalculate breath info:
for i=1:length(BB_i_start)
    [~,tempi]=max(vol2_filtered1(BB_i_start(i):BB_i_end(i)));
    BB_i_mid(i,1)=BB_i_start(i)+tempi-1;
end
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;

criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
Vdot = VT./Ttot;
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));



%% If we find a small(er) expiration above FRC presume the breath is part of one continuous inspiration:

%VolStartInsp = (VTi-VTe)/VTmean; 
FVTe = VTe./VTi;
criteria=find(((FVTe<0.2)&Te<normalTe/2)); %VolStartInsp>0.33 %was leaving double peaks on smaller breaths with incl exhalation

if ~isempty(criteria)
NN=length(BB_i_start);
lastbreathisbad=criteria(end)==NN;
while lastbreathisbad
        criteria(end)=[];
        BB_i_start(end)=[];
        BB_i_mid(end)=[];
        BB_i_end(end)=[];
        NN=NN-1;
       if isempty(criteria)
            break
       end
        lastbreathisbad=criteria(end)==NN;        
end
   criteria2=criteria+1;
   BB_i_start(criteria2)=[];
   BB_i_mid(criteria)=[];
   BB_i_end(criteria)=[];
end

%% recalculate breath info:
for i=1:length(BB_i_start)
    [~,tempi]=max(vol2_filtered1(BB_i_start(i):BB_i_end(i)));
    BB_i_mid(i,1)=BB_i_start(i)+tempi-1;
end
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;

criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
Vdot = VT./Ttot;
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));


% %% If there is small(er) inspiration after previous incomplete exhalation (mid-inspiratory pause), merge breath with previous expiration:
% 
% VolStartInsp = (VTi-VTe)/VTmean; VolStartInsp = [NaN VolStartInsp(1:end-1)];
% criteria=find((VolStartInsp>0.33&VTi<(VTi_thres*4)));
% 
% 
% if ~isempty(criteria)
% firstbreathisbad=criteria(1)==1;
% while firstbreathisbad
%    criteria(1)=[];
%    BB_i_start(1)=[];
%    BB_i_mid(1)=[];
%    BB_i_end(1)=[];   
%    criteria=criteria-1;   
%    if isempty(criteria)
%        break
%    end
%    firstbreathisbad=criteria(1)==1;
% end
%    criteria2=criteria-1;
%    BB_i_start(criteria)=[];
%    BB_i_mid(criteria)=[];
%    BB_i_end(criteria2)=[]; 
% end
% 
% %% recalculate breath info:
% vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
%     vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
%     vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
% volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
%     volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
%     volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest
% 
% if ~minimum_figs
%      
%     ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
%     ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow,time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
%     linkaxes(ax1,'x');
% end
% volupperminusvollower = volupper-vollower;
% VT = volupperminusvollower(BB_i_mid);
% VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
% VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
% Ti = (BB_i_mid-BB_i_start)'*dt;
% Te = (BB_i_end-BB_i_mid)'*dt;
% Ttot = Ti+Te;
% 
% criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
% Vdot = VT./Ttot;
% normalTtot=median(Ttot(criteriafornormal));
% normalTi=median(Ti(criteriafornormal));
% normalTe=median(Te(criteriafornormal));
% 
% 
% chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));


%% More
VE = VTe./Ttot;
VI = VTi./Ttot;



%% If we find a small expiration presume the breath is part of one continuous inspiration:
if 0
criteria=find(VTe<VTe_thres);

if ~isempty(criteria)
NN=length(BB_i_start);
lastbreathisbad=criteria(end)==NN;
while lastbreathisbad
        criteria(end)=[];
        BB_i_start(end)=[];
        BB_i_mid(end)=[];
        BB_i_end(end)=[];
        NN=NN-1;
       if isempty(criteria)
            break
       end
        lastbreathisbad=criteria(end)==NN;        
end
   criteria2=criteria+1;
   BB_i_start(criteria2)=[];
   BB_i_mid(criteria)=[];
   BB_i_end(criteria)=[];
end

% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;

criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
Vdot = VT./Ttot;
normalTtot=median(Ttot(criteriafornormal));
normalTi=median(Ti(criteriafornormal));
normalTe=median(Te(criteriafornormal));


chksum=sum(BB_i_start(2:end)~=BB_i_end(1:end-1));
end
%% Find probable apneic periods using original unfiltered flow
if 1&&length(BB_i_start)>5
    meanVflow = mean(Vflow(BB_i_start(1):BB_i_end(end)));
    slidew = normalTtot*2;
    if isnan(normalTtot)
        normalTtot=4;
    end
    dTi=25; dT=dt*dTi;
    slidewi = round(slidew/dt);
    Nwindows = ceil((length(Vflow)-slidewi)/dTi)+1;
    normalInspFlow = VTi/Ti;
    thres1 = normalInspFlow/3; %%%%%%%%%%%%%%%%%%%%%%%% low amplitude 
    thres2 = normalInspFlow/2; %%%%%%%%%%%%%%%%%%%%%%%% near flow = 0 (mean)
    Vflow_delta = zeros(Nwindows,1);
    Vflow_median = zeros(Nwindows,1);
    Vflow_apnea = 0*Vflow;
    for i=1:Nwindows
        li = 1+(i-1)*dTi;
        ri = li+slidewi;
        if ri>length(Vflow),ri=length(Vflow); end
        Vflow_delta(i)=prctile(Vflow(li:ri),100)-prctile(Vflow(li:ri),0);
        Vflow_median(i)=prctile(Vflow(li:ri),50);
        if Vflow_delta(i)<thres1&&(Vflow_median(i)-meanVflow)<thres2&&(Vflow_median(i)-meanVflow)>-thres2
            Vflow_apnea((li+round(0.10*slidewi)-1):(li+round(0.90*slidewi))) = 1;
        end
    end
    Vflow_apnea = Vflow_apnea(1:length(Vflow));
    
    trendx = Vflow_filtered1-Vflow;
    Vflow_apnea0 = Vflow_apnea;
%     if ~minimum_figs
%     ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,trendx,time,Vflow_apnea0,time,Vflow,time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
%     linkaxes(ax1,'x');
%     end
%     lengththreshold = round(normalTtot*0.25/dt); %if larger (e.g. 0.75) flat expiratory patterns get merged together into apnea (even if breathing is normal)
%     
%     Idiffapnea = find(diff(Vflow_apnea)==1);
%     Idiffnapnea = find(diff(Vflow_apnea)==-1);
%     if ~isempty(Idiffnapnea)&&~isempty(Idiffapnea)
%         if Idiffnapnea(1)<Idiffapnea(1), Idiffapnea = [1 Idiffapnea]; end
%         if Idiffnapnea(end)<Idiffapnea(end), Idiffnapnea = [Idiffnapnea length(Vflow_apnea)]; end
%         apnealengths = Idiffnapnea-Idiffapnea;
%         napnealengths = Idiffapnea(2:end)-Idiffnapnea(1:end-1);
%         
%         minnapnea=min(napnealengths);
%         minapnea=min(apnealengths);
%         while (~isempty(minnapnea)&&minnapnea<lengththreshold)||(~isempty(minapnea)&&minapnea<lengththreshold)
%             apneashorter=minapnea<minnapnea;
%             
%             if apneashorter
%                 [minapnea,tempi]=min(apnealengths); minnapnea=min(napnealengths); %min1=minapnea; mincriteria=minapnea;
%             else
%                 [minnapnea,tempi]=min(napnealengths); minapnea=min(apnealengths); %min1=minnapnea; mincriteria=minnapnea;
%             end
%             if minnapnea>=lengththreshold&&minapnea>=lengththreshold
%                 break
%             end
%                 
%             %[minnapnea,minapnea]
%             %while min1==mincriteria %remove all apneas/napneas of the same length
%                 if apneashorter %remove shortest apnealengths
%                     li = Idiffapnea(tempi)+1;
%                     ri = Idiffnapnea(tempi);
%                     chk = mean(Vflow_apnea(li:ri));
%                     Vflow_apnea(li:ri) = 0;
%                 else %remove shortest non-apnealengths
%                     li = Idiffnapnea(tempi)+1;
%                     ri = Idiffapnea(tempi+1);
%                     chk = mean(Vflow_apnea(li:ri));
%                     Vflow_apnea(li:ri) = 1;
%                 end
%                 if 1
%                 if ~minimum_figs
%                 ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,trendx,time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
%                 linkaxes(ax1,'x');
%                 end
%                 end
%                 Idiffapnea = find(diff(Vflow_apnea)==1);
%                 Idiffnapnea = find(diff(Vflow_apnea)==-1);
%                 if Idiffnapnea(1)<Idiffapnea(1), Idiffapnea = [1 Idiffapnea]; end
%                 if Idiffnapnea(end)<Idiffapnea(end), Idiffnapnea = [Idiffnapnea length(Vflow_apnea)]; end
%                 apnealengths = Idiffnapnea-Idiffapnea;
%                 napnealengths = Idiffapnea(2:end)-Idiffnapnea(1:end-1);
%                 
%                 minnapnea=min(napnealengths);
%                 minapnea=min(apnealengths);
%                 apneashorter=minapnea<minnapnea;
%                 if apneashorter
%                     [~,tempi]=min(apnealengths); mincriteria=minapnea;
%                 else
%                     [~,tempi]=min(napnealengths); mincriteria=minnapnea;
%                 end
% %             end
%         end
%     end
if ~minimum_figs
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,trendx,time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
else
    Vflow_apnea = time*0;
    Vflow_apnea0 = time*0;
end %if detect apneas

%% Inspirations occurring entirely 95% in apneas are considered apneas -- delete / merge with prior expiration
    if 0
    Apnea_B1 = 0*BB_i_start;

    for i=1:length(BB_i_start)
        if mean(Vflow_apnea0(BB_i_start(i):BB_i_mid(i)))>0.95
            Apnea_B1(i)=1;
        end
    end
    
    criteria = find(Apnea_B1==1);
    if ~isempty(criteria)
        firstbreathisbad=criteria(1)==1;
        while firstbreathisbad
            criteria(1)=[];
            BB_i_start(1)=[];
            BB_i_mid(1)=[];
            BB_i_end(1)=[];
            criteria=criteria-1;
            if isempty(criteria)
                break
            end
            firstbreathisbad=criteria(1)==1;
        end
        criteria2=criteria-1;
        BB_i_start(criteria)=[];
        BB_i_mid(criteria)=[];
        BB_i_end(criteria2)=[];
    end
    
% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,trendx,time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;
Vdot = VT./Ttot;
criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    end
%% In breaths with apneas - move the end-inspiration to the first zero crossing
%wtypically 
    %criteria longer than normal Ti, Fapnea>67%, existence of an earlier peak.
Fapnea = NaN*BB_i_start;
peak_i_earliest = NaN*BB_i_start;
Earlierpeak = NaN*BB_i_start;

CurrentPeakinApnea = Vflow_apnea0(BB_i_mid)';
peak_i_current = BB_i_mid-BB_i_start+1;
for i=1:length(BB_i_start)
    Fapnea(i,1)=mean(Vflow_apnea0(BB_i_start(i):BB_i_end(i)));
    tempvol = vol2_filtered1(BB_i_start(i):BB_i_end(i));
    [max_list,min_list] = peakdet((tempvol),0.2*std(tempvol)); 
    if isempty(max_list)
        peak_i_earliest(i,1)=peak_i_current(i);
        Earlierpeak(i,1) = 0;
        continue
    end
    peak_i_earliest(i,1)=max_list(1,1);
    
    Earlierpeak(i,1) = peak_i_earliest(i)<peak_i_current(i);
    
    %figure(100); plot(tempvol); hold('on'); 
end

criteria = Fapnea>0.67&Earlierpeak&Ttot'>normalTtot&CurrentPeakinApnea==1;

BB_i_mid(criteria)=BB_i_start(criteria)+peak_i_earliest(criteria)-1;

% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,trendx,time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;
Vdot = VT./Ttot;
criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
normalTtot=median(Ttot(criteriafornormal));


%% Detect start of inspiration at the end of apnea based on the time to 50% of peak flow. 
Fpeakflow = 0.5;
%find 5percent volume time
for i=1:length(BB_i_start)
    tempflow = Vflow_filtered1(BB_i_start(i):BB_i_mid(i));
    temppeakflow = max(tempflow);
    BB_Ti_X(i) = dt*(find(tempflow>temppeakflow*Fpeakflow,1)-1);
end

Ti_X_median = median(BB_Ti_X);
Ti_X_upper = prctile(BB_Ti_X,90);
Ti_X_75 = prctile(BB_Ti_X,75);

for i=1:length(BB_i_start)
    Fearlyinspisapnea(i)=mean(Vflow_apnea(BB_i_start(i):(BB_i_start(i)+round(BB_Ti_X(i)/dt))));
    if BB_Ti_X(i)>Ti_X_upper&&Fearlyinspisapnea(i)>0.5 %at least half of early insp is apnea
        BB_i_start(i)=BB_i_start(i)+round((BB_Ti_X(i)-Ti_X_upper)/dt);
        if i>1
            BB_i_end(i-1)=BB_i_start(i);
        end
    end
end
    
%% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,trendx,time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;
Vdot = VT./Ttot;
criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
normalTtot=median(Ttot(criteriafornormal));


%% Shallow breaths occurring predominantly in detected apneas with low ventilation is an apnea -- delete / merge with prior expiration
    
Apnea_B1 = 0*BB_i_start;
for i=1:length(BB_i_start)
    if mean(Vflow_apnea0(BB_i_start(i):BB_i_end(i)))>0.67&&Vdot(i)/(VTmean/normalTtot)<0.05
        Apnea_B1(i)=1;
    end
end

criteria = find(Apnea_B1==1);
if ~isempty(criteria)
    firstbreathisbad=criteria(1)==1;
    while firstbreathisbad
        criteria(1)=[];
        BB_i_start(1)=[];
        BB_i_mid(1)=[];
        BB_i_end(1)=[];
        criteria=criteria-1;
        if isempty(criteria)
            break
        end
        firstbreathisbad=criteria(1)==1;
    end
    criteria2=criteria-1;
    BB_i_start(criteria)=[];
    BB_i_mid(criteria)=[];
    BB_i_end(criteria2)=[];
end
    
%% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,trendx,time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;
Vdot = VT./Ttot;
criteriafornormal = VT>(0.2*VTmean)&Ttot<10;

%% Long expirations are apneas (later broken up into smaller pieces)
%But they should also contain detected apneas (Vflow_apnea0), otherwise slow large expirations get caught up
    FtotIsAnApnea = 1.0;
    VEaverage = VTmean/normalTtot;
    TeThresholdForApnea = (0.5*normalTe+FtotIsAnApnea*normalTtot);
    Apnea_B = Te'>TeThresholdForApnea;
    Iapnea = find(Apnea_B);
    for i=length(Iapnea):-1:1
        Te_total = dt*(BB_i_end(Iapnea(i))-BB_i_mid(Iapnea(i)));
%         figure(100)
%         subplot(211); plot(Vflow_filtered1(BB_i_mid(Iapnea(i)):BB_i_end(Iapnea(i))));
%         subplot(212); plot(vol2_filtered1(BB_i_mid(Iapnea(i)):BB_i_end(Iapnea(i)))); hold('on');
%         plot(Vflow_apnea0(BB_i_mid(Iapnea(i)):BB_i_end(Iapnea(i))));
%         
        clear temp
            temp = find(Vflow_apnea0(BB_i_mid(Iapnea(i)):BB_i_end(Iapnea(i))),1);
            Te_est1 = dt*(temp-1);
            if isempty(temp)
                Apnea_B(Iapnea(i))=0;
                continue
            end
            Te_est2 = normalTe; %previous version used this
            Te_est = max([Te_est1 Te_est2]); %will avoid the possiblity of apnea starting at end insp.
        
            if (Te_total-Te_est)>FtotIsAnApnea*normalTtot     %if the new period is now too short to be called apnea, do not proceed        
                BB_i_start_newapneabreath = BB_i_mid(Iapnea(i))+round(Te_est/dt);
                di=BB_i_end(Iapnea(i))-BB_i_start_newapneabreath;
                BB_i_mid_newapneabreath = BB_i_start_newapneabreath+round(normalTi/normalTtot*di);
                
                BB_i_start = [BB_i_start(1:Iapnea(i));BB_i_start_newapneabreath;BB_i_start(Iapnea(i)+1:end)];
                BB_i_end = [BB_i_end(1:Iapnea(i)-1);BB_i_start_newapneabreath;BB_i_end(Iapnea(i):end)];
                BB_i_mid = [BB_i_mid(1:Iapnea(i));BB_i_mid_newapneabreath;BB_i_mid(Iapnea(i)+1:end)];
                Apnea_B = [Apnea_B(1:Iapnea(i)-1);0;1;Apnea_B(Iapnea(i)+1:end)];
                
                if minimum_figs
                ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:',time(BB_i_start_newapneabreath),vol2_filtered1(BB_i_start_newapneabreath),'ro',time(BB_i_mid_newapneabreath),vol2_filtered1(BB_i_mid_newapneabreath),'ko');
                ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,trendx,time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
                linkaxes(ax1,'x');
                end
            else
                Apnea_B(Iapnea(i))=0;
            end
    end

%% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs
     
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow_filtered1,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;
Vdot = VT./Ttot;
criteriafornormal = VT>(0.2*VTmean)&Ttot<10;


%% Divide up apneas into Ttot-sized parcels
Iapnea = find(Apnea_B);
Nbreathsapnea=floor(Ttot(Iapnea)/normalTtot);
for i=length(Iapnea):-1:1
    if Nbreathsapnea(i)>0
    di=round((BB_i_end(Iapnea(i))-BB_i_start(Iapnea(i)))/Nbreathsapnea(i));
%     time(BB_i_start(Iapnea(i)))
    temp = (BB_i_start(Iapnea(i)):di:((BB_i_end(Iapnea(i)))+di/4))'; %di/5 just handles rounding issues
    BB_i_start = [BB_i_start(1:Iapnea(i));temp(2:(end-1));BB_i_start((Iapnea(i)+1):end)];
    BB_i_end = [BB_i_end(1:Iapnea(i)-1);temp(2:(end-1));BB_i_end(Iapnea(i):end)];
    BB_i_mid = [BB_i_mid(1:Iapnea(i)-1);temp(1:(end-1))+round((normalTi/normalTtot*di));BB_i_mid((Iapnea(i)+1):end)];
    Apnea_B = [Apnea_B(1:Iapnea(i));1+0*temp(2:(end-1));Apnea_B((Iapnea(i)+1):end)];
    end
end


%% recalculate breath info:
vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
volupper = interp1(time(BB_i_mid),vol2_filtered1(BB_i_mid),time,'linear');  
    volupper(1:BB_i_mid(1)-1)=vol2_filtered1(BB_i_mid(1)); %extrapolation using nearest
    volupper(BB_i_mid(end):end)=vol2_filtered1(BB_i_mid(end)); %extrapolation using nearest

if ~minimum_figs 
    ax1(XX)=subplot(Nloops+1,1,XX); plot(time,vol2_filtered1,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,volupper,'k:',time,vollower,'r:');
    ax1(Nloops+1)=subplot(Nloops+1,1,Nloops+1); plot(time,Vflow,time,Vflow_filtered1,time,Vflow_apnea0,time,Vflow_apnea,time(BB_i_start),Vflow_filtered1(BB_i_start),'r.',time(BB_i_mid),Vflow_filtered1(BB_i_mid),'k.');
    hold('on'); stairs(time(BB_i_start),Apnea_B*max(Vflow)/2,'b');  hold('off');
    linkaxes(ax1,'x');
end
volupperminusvollower = volupper-vollower;
VT = volupperminusvollower(BB_i_mid);
VTi = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_start);
VTe = vol2_filtered1(BB_i_mid) - vol2_filtered1(BB_i_end);
Ti = (BB_i_mid-BB_i_start)'*dt;
Te = (BB_i_end-BB_i_mid)'*dt;
Ttot = Ti+Te;
VE = VTe./Ttot;
VI = VTi./Ttot;
Vdot = VT./Ttot;

%% Additional information
% baselinedriftstd = std(vollower+voldetrend1);
% Endexpvol = vol2_filtered1(BB_i_start);
% Endexpvol_t = time(BB_i_start);
temp = polyfit(time,vollower,1);
leak = temp(1) + leak;
% if ~minimum_figs
% figure(99);
% plot(time,voldetrend1+vol2_filtered1,time(BB_i_start),voldetrend1(BB_i_start)+vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),voldetrend1(BB_i_mid)+vol2_filtered1(BB_i_mid),'k.',time,voldetrend1+volupper,'k:',time,voldetrend1+vollower,'r:');
% title(['leak=' num2str(leak) ', baselinedriftstd=' num2str(baselinedriftstd)]);
% set(gcf,'position',[50 100 500 500]);
% end

% vollower_record = - voldetrend2;
%baselinedriftstd = std(vollower_record);

% figure(900);
% plot(time,vol2_filtered1+vollower_record,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,vollower_record,'k:');
% set(gcf,'position',[50 100 500 500]);

%% detrend and re-run
if XX<Nloops
    %vol2_filtered1_backup = vol2_filtered1;
    %vol2_filtered1
    volarray = vol2_filtered1([BB_i_start;BB_i_end(end)]);
    %volarray = [volarray(1) volarray volarray(end)];
    %medfilt1
    %volarray2 = medfilt1(volarray,1)
    %volarray2 = volarray2(2:end-1);
    voltrend = interp1(time([BB_i_start;BB_i_end(end)]),volarray,time,'linear');        
    voltrend(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1));
    voltrend(BB_i_end(end):end)=vol2_filtered1(BB_i_end(end));
    %figure(102);plot(time,[vol2_filtered1;voltrend]);
    vol2_filtered1=vol2_filtered1-voltrend*0.33;
%     vollower_record = vollower_record + voltrend;
end

end
%% Force V positive, just in case.

VT(VT<0)=0;
Vdot(Vdot<0)=0;
VTi(VTi<0)=0;
VTe(VTe<0)=0;
VE(VE<0)=0;
VI(VI<0)=0;

%% Breaths with detected apneas have zero ventilation
VT(Apnea_B==1)=0;
Vdot(Apnea_B==1)=0;
VTi(Apnea_B==1)=0;
VTe(Apnea_B==1)=0;
VE(Apnea_B==1)=0;
VI(Apnea_B==1)=0;

%% Start/End breath indices in original indices
clear I;
I.starti = BB_i_start;
I.midi = BB_i_mid;
I.endi = BB_i_end;


function [x_peak,y_peak] = PeakFitQuadratic(x,y)

%x = [aD bD cD];
%y = [AD BD CD];

a = [y(3)*x(2)-y(3)*x(1)-x(3)*y(2)+x(3)*y(1)-y(1)*x(2)+x(1)*y(2)]/[(x(2)-x(1))*(-x(2)*x(3)+x(2)*x(1)+x(3)^2-x(1)*x(3))];
b = [y(2)-a*x(2)^2-y(1)+a*x(1)^2]/[x(2)-x(1)];
c = y(1)-a*x(1)^2-b*x(1);

%x1 = 2:0.001:2.02;
%y1 = a*x1.^2+b*x1+c;

x_peak = -b/(2*a);
y_peak = a*x_peak.^2+b*x_peak+c;

%figure(200), plot(x,y,'.',x1,y1,x_peak,y_peak,'o');
