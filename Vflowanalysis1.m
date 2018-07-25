function [I,VI,VTi,Ti,Te,VFlow_filtered1_aligned,VFlow_filtered2]=Vflowanalysis1(Vflow,time,dt,minimum_figs)

%% Respiratory Trace Analyse
%[Pxy,F] = CPSD(X,Y,WINDOW,NOVERLAP,NFFT,Fs)
% Respiratory Trace Analyse
detectapneas=1;
detrend_flow=1;

%transpose if needed:
if size(Vflow,2)<size(Vflow,1)
    Vflow=Vflow';
end

if detrend_flow
    Vflow=Vflow-mean(Vflow);
end

T0 = length(time)*dt;
minF=1/10;

[Pyy,F] = cpsd(Vflow,Vflow,length(Vflow),0,length(Vflow),1/dt);
Pyy=Pyy(F>minF);
F=F(F>minF);

if ~minimum_figs
    figure(1), plot(F,Pyy);
end
[Y,I] = max(Pyy); %%Use a parametric method to do this in the future.
BBpsd = 1/F(I);

%Filter gently for use in analysis
filter_LFcutoff_butter1 = 0.05;
filter_HFcutoff_butter1 = 10;
filter_order1 = 2;
[B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_HFcutoff_butter1]/(1/dt/2),'low');
[B_butterLcut,A_butterLcut] = butter(filter_order1,[filter_LFcutoff_butter1]/(1/dt/2),'high');
%[B_butter1,A_butter1] = butter(filter_order1,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dt/2));

VFlow_filteredHcut = filtfilt(B_butterHcut,A_butterHcut,Vflow);
VFlow_filtered1 = filtfilt(B_butterLcut,A_butterLcut,VFlow_filteredHcut);

fixedgeissuewithhighpass=0;
lefttrim_i=1;
righttrim_i=1;
if fixedgeissuewithhighpass
lowpassX=VFlow_filteredHcut-VFlow_filtered1;
ax2(2)=subplot(2,1,1); plot(time,Vflow,time,VFlow_filtered1,time,VFlow_filteredHcut,time,lowpassX);

for i=3:length(VFlow_filteredHcut)
    if VFlow_filteredHcut(i)>0&&VFlow_filteredHcut(i-1)<0
        break
    end
end
lefttrim_i=i;
for i=(length(VFlow_filteredHcut)-2):-1:1
    if VFlow_filteredHcut(i+1)>0&&VFlow_filteredHcut(i)<0
        break
    end
end
righttrim_i=i;

VFlow_filteredHcutTrim=VFlow_filteredHcut(lefttrim_i:righttrim_i);
timetrim=time(lefttrim_i:righttrim_i);
%figure(345);plot(timetrim,VFlow_filteredHcutTrim);

VFlow_filtered1trim = filtfilt(B_butterLcut,A_butterLcut,VFlow_filteredHcutTrim);
lowpassXtrim = VFlow_filteredHcutTrim-VFlow_filtered1trim;

lowpassXtrim_nearest = interp1(timetrim,lowpassXtrim,time,'nearest','extrap');
VFlow_filtered1_V2 = VFlow_filteredHcut - lowpassXtrim_nearest;
%figure(345);plot(timetrim,VFlow_filteredHcutTrim,timetrim,VFlow_filtered1trim,timetrim,lowpassXtrim,time,lowpassXtrim_nearest,time,VFlow_filtered1_V2);

VFlow_filtered1=VFlow_filtered1_V2;
end
VFlow_filtered1_aligned = [NaN*zeros(1,lefttrim_i-1) VFlow_filtered1 NaN*zeros(1,righttrim_i-1)];

ax2(2)=subplot(2,1,1); plot(time,Vflow,time,VFlow_filtered1);

%%

%% Filter strongly to aid peak identification
[Presp,F1] = pwelch(VFlow_filtered1-mean(VFlow_filtered1),length(Vflow),0,length(Vflow),1/dt);
df=1/T0;
guess_resp_freq=1/4;
Fi=floor((guess_resp_freq/2)/df+1); %25 sec breaths
Fii=ceil((guess_resp_freq*2)/df+1); %0.5 sec breaths
[temp,index_f]=max(Presp(Fi:Fii));
index_f2=index_f+Fi-1;
%peak fit quadratic:
range3 = (index_f2-1):(index_f2+1);
[fresp_peak,Presp_peak] = PeakFitQuadratic(F1(range3),Presp(range3));

%Filter strongly to aid peak identification, using LFcutoff as half the expected,
%and HF cutoff as twice the expected
filter_LFcutoff_butter2 = fresp_peak*.2;
filter_HFcutoff_butter2 = fresp_peak*1.67;

filter_order2 = 2;
[B_butter2,A_butter2] = butter(filter_order2,[filter_LFcutoff_butter2 filter_HFcutoff_butter2]/(1/dt/2));
VFlow_filtered2 = filtfilt(B_butter2,A_butter2,Vflow);

VFlow_filtered_threshold = 0;

if ~minimum_figs
    figure(200), plot(time,VFlow_filtered1,time,VFlow_filtered2)
end

%%detecting max flow values (i.e. expiration..)
clear VFlow_crossing BB_i_start VFlow_crossing_t;
count=1;
for i=1:length(Vflow),
    if (i>1)&&(i<length(Vflow))
        if (VFlow_filtered2(i)>VFlow_filtered2(i-1))&&(VFlow_filtered2(i)>VFlow_filtered2(i+1))&&(VFlow_filtered2(i)>VFlow_filtered_threshold)
            BB_i_start(count)=i;
            count=count+1;
        end
    end
end


if ~minimum_figs
    figure(3); subplot(3,1,1); plot(time,VFlow_filtered1,time,VFlow_filtered2,time(BB_i_start),VFlow_filtered2(BB_i_start),'.')
end
%Store data
VFlow_crossing_step1_i=BB_i_start;

% start from peak data and move leftwards to find zero crossing:
for count=1:length(BB_i_start)
    for k=1:999999
        if (VFlow_filtered2(BB_i_start(count))>=0)&&(VFlow_filtered2(BB_i_start(count)-1)<0)
            break;
        else
            if BB_i_start(count)==2 %hitting left hand limit
                break
            else
                BB_i_start(count)=BB_i_start(count)-1;%moving left...
            end
        end
    end
end
if BB_i_start(1)<3
    BB_i_start(1)=[];
end

i_=0;
for i=1:length(BB_i_start)
    if i>1
        if BB_i_start(i-i_)==BB_i_start(i-i_-1)
            BB_i_start(i-i_)=[];
            i_=i_+1;
        end
    end
end

if ~minimum_figs
    figure(3); subplot(3,1,2); plot(time,VFlow_filtered1,time(BB_i_start),VFlow_filtered1(BB_i_start),'.r')
end
%% Shift from start inspiration, head rightwards to find end inspiration
% 
% clear BB_i_mid
% BB_i_mid = BB_i_start+2;
% 
% for count=1:length(BB_i_start)
%     for k=1:999999
%         if (VFlow_filtered1(BB_i_mid(count))<0)&&(VFlow_filtered1(BB_i_mid(count)-1)>=0)||(BB_i_mid(count)>=length(VFlow_filtered1)) %last OR ensures you don't exceed the length of the trace
%             break;
%         else
%             if BB_i_mid(count)==2
%                 break
%             else
%                 BB_i_mid(count)=BB_i_mid(count)+1; %increase index
%             end
%         end
%     end
% end
% 
% i_=0;
% for i=1:length(BB_i_mid)
%     if i>1
%         if BB_i_mid(i-i_)==BB_i_mid(i-i_-1)
%             BB_i_mid(i-i_)=[];
%             i_=i_+1;
%         end
%     end
% end
clear BB_i_mid
for i=1:(length(BB_i_start)-1)
    [maxvalue,index]=max(cumsum(VFlow_filtered2(BB_i_start(i):(BB_i_start(i+1)-1))));
    BB_i_mid(i) = BB_i_start(i)+index-1;
end

if ~minimum_figs
    figure(3); subplot(3,1,3); plot(time,VFlow_filtered1,time(BB_i_start),VFlow_filtered1(BB_i_start),'.r',time(BB_i_mid),VFlow_filtered1(BB_i_mid),'.k')
end

%% Where is the zero crossing on the less filtered trace?
if 0
rangeleftright = 0.5;% look left and right by e.g. 0.5 s to find a better start of inspiration/expiration.
rangeleftrighti = round(0.5/dt);
A = VFlow_filtered2(BB_i_start)>VFlow_filtered1(BB_i_start);
for i=1:(length(BB_i_start)-1);
    lefti = (BB_i_start(i)-rangeleftrighti);
    righti = (BB_i_start(i)+rangeleftrighti);
    if lefti<1, lefti=1; end
    if righti>length(VFlow_filtered2), righti=length(VFlow_filtered2); end
    I_ = lefti:righti;
    
    if A(i)>0
        index=find(VFlow_filtered1(BB_i_start(i):righti)>0,1,'first');
        if ~isempty(index)
            BB_i_start(i) = BB_i_start(i)+index-1;
        end
    else
        index=find(VFlow_filtered1(lefti:BB_i_start(i))<0,1,'last');
        if ~isempty(index)
            BB_i_start(i) = BB_i_start(i)-lefti+index-1;
        end
    end
    if isempty(index)
    %where is the lowest value
    [~,index]=min(abs(VFlow_filtered1(I_)));
    BB_i_start(i) = I_(1)+index-1;
    end
    
    [maxvalue,index]=max(cumsum(VFlow_filtered1(BB_i_start(i):(BB_i_start(i+1)-1))));
    BB_i_mid(i) = BB_i_start(i)+index-1;
end

end


%% Breath-breath data

clear BB BB_t Ti Te BB_i_end VTi VTe;
for i=1:(length(BB_i_start)-1)
        BB_i_end(i)=BB_i_start(i+1);
        Ti(i)=time(BB_i_mid(i))-time(BB_i_start(i));
        Te(i)=time(BB_i_start(i+1))-time(BB_i_mid(i));%VFlow_crossing_t(i+1)-VFlow_crossingE_t(i);
        VTi(i)=dt*sum(VFlow_filtered1(BB_i_start(i):(BB_i_mid(i)-1))); % sum(VFlow_filtered1_pos(BB_i_start(i):(BB_i_start(i+1)-1)));
        VTe(i)=-dt*sum(VFlow_filtered1(BB_i_mid(i):(BB_i_end(i)-1)));
end
BB=Ti+Te;
BB_i_start(end)=[]; %may move this to the end -- might be easier to keep here
BB_t=time(BB_i_start);

if ~minimum_figs
    figure(2);
    ax2(2)=subplot(2,1,1); plot(time,VFlow_filtered1,time(BB_i_start),VFlow_filtered1(BB_i_start),'.r',time(BB_i_mid),VFlow_filtered1(BB_i_mid),'.k');
end


%% Remove extreme mini-breaths
if 1
VTi_thres = prctile(VTi,90)/50;
VTe_thres = prctile(VTe,90)/50;
Ti_thres = BBpsd/20;
Te_thres = BBpsd/20;

%If inspiration and expiration are both either short or small, remove the whole breath (rare)
criteria=find(((VTi<VTi_thres)|(Ti<Ti_thres))&((VTe<VTe_thres)|(Te<Te_thres)));
BB_i_end(criteria)=[];
BB_i_mid(criteria)=[];
BB_i_start(criteria)=[];
%recalculate breath info:
clear BB BB_t Ti Te VTi VTe;
for i=1:length(BB_i_start)
        Ti(i)=time(BB_i_mid(i))-time(BB_i_start(i));
        Te(i)=time(BB_i_end(i))-time(BB_i_mid(i));%VFlow_crossing_t(i+1)-VFlow_crossingE_t(i);
        VTi(i)=dt*sum(VFlow_filtered1(BB_i_start(i):(BB_i_mid(i)-1))); % sum(VFlow_filtered1_pos(BB_i_start(i):(BB_i_start(i+1)-1)));
        VTe(i)=-dt*sum(VFlow_filtered1(BB_i_mid(i):(BB_i_end(i)-1)));
end
BB=Ti+Te;
BB_t=time(BB_i_start);
VE = VTe./BB;
VI = VTi./BB;

%If we find a small/short inspiration, presume it is part of the previous expiration:
criteria=find((VTi<VTi_thres)|(Ti<Ti_thres));
while ~isempty(criteria)&&criteria(1)==1
    BB_i_start(1)=[];
    BB_i_mid(1)=[];
    BB_i_end(1)=[];
    criteria(1)=[];
    criteria=criteria-1;
    if isempty(criteria)
        break
    end
end
    BB_i_start(criteria)=[];
    BB_i_mid(criteria)=[];
    criteria2=criteria-1;
    BB_i_end(criteria2)=[];
%recalculate breath info:
clear BB BB_t Ti Te VTi VTe;
for i=1:length(BB_i_start)
        Ti(i)=time(BB_i_mid(i))-time(BB_i_start(i));
        Te(i)=time(BB_i_end(i))-time(BB_i_mid(i));%VFlow_crossing_t(i+1)-VFlow_crossingE_t(i);
        VTi(i)=dt*sum(VFlow_filtered1(BB_i_start(i):(BB_i_mid(i)-1))); % sum(VFlow_filtered1_pos(BB_i_start(i):(BB_i_start(i+1)-1)));
        VTe(i)=-dt*sum(VFlow_filtered1(BB_i_mid(i):(BB_i_end(i)-1)));
end
BB=Ti+Te;
BB_t=time(BB_i_start);
VE = VTe./BB;
VI = VTi./BB;

%If we find a small/short expiration, presume it is part of the previous (and next) inspiration:
criteria=find((VTe<VTe_thres)|(Te<Te_thres));
while ~isempty(criteria)&&(criteria(end)==length(BB_i_start))
    BB_i_start(end)=[];
    BB_i_mid(end)=[];
    BB_i_end(end)=[];
    criteria(end)=[];
    if isempty(criteria)
        break
    end  
end
BB_i_end(criteria)=[];
BB_i_mid(criteria)=[];
criteria2=criteria+1;
BB_i_start(criteria2)=[];

if ~minimum_figs
    figure(2);
    ax2(2)=subplot(2,1,1); plot(time,VFlow_filtered1,time(BB_i_start),VFlow_filtered1(BB_i_start),'.r',time(BB_i_mid),VFlow_filtered1(BB_i_mid),'.k');
end

%recalculate breath info:
clear BB BB_t Ti Te VTi VTe;
for i=1:length(BB_i_start)
        Ti(i)=time(BB_i_mid(i))-time(BB_i_start(i));
        Te(i)=time(BB_i_end(i))-time(BB_i_mid(i));%VFlow_crossing_t(i+1)-VFlow_crossingE_t(i);
        VTi(i)=dt*sum(VFlow_filtered1(BB_i_start(i):(BB_i_mid(i)-1))); % sum(VFlow_filtered1_pos(BB_i_start(i):(BB_i_start(i+1)-1)));
        VTe(i)=-dt*sum(VFlow_filtered1(BB_i_mid(i):(BB_i_end(i)-1)));
end
BB=Ti+Te;
BB_t=time(BB_i_start);
VE = VTe./BB;
VI = VTi./BB;

if ~minimum_figs
    figure(2);
    ax2(2)=subplot(2,1,1); plot(time,VFlow_filtered1,time(BB_i_start),VFlow_filtered1(BB_i_start),'.r',time(BB_i_mid),VFlow_filtered1(BB_i_mid),'.k');
end
end


%% Start/End breath indices in original indices
clear I;
I.starti = BB_i_start - lefttrim_i + 1;
I.midi = BB_i_mid - lefttrim_i + 1;
I.endi = BB_i_end - lefttrim_i + 1;

%% Long breaths are broken up into smaller pieces

Ttottemp=Ti+Te;
TionTtottemp=Ti./Ttottemp;
TionTtot=median(TionTtottemp(Ttottemp<6));
medianTtot = median(Ttottemp(Ttottemp<6));
TtotApneaThres = median(Ttottemp(Ttottemp<6));

T0=sum(Ttottemp);
VEaverage=sum(VI.*Ttottemp)/T0;

Iapnea=find((VI/VEaverage<0.2)&Ttottemp>TtotApneaThres);

Nbreathsapnea=round(Ttottemp(Iapnea)/medianTtot);
for i=length(Iapnea):-1:1
    di=round((I.endi(Iapnea(i))-I.starti(Iapnea(i)))/Nbreathsapnea(i));
    di_TiTtot=round((TionTtot*di));
    temp = I.starti(Iapnea(i)):di:((I.endi(Iapnea(i)))+di/5);
    I.starti = [I.starti(1:Iapnea(i)) temp(2:(end-1)) I.starti((Iapnea(i)+1):end)];
    I.endi = [I.endi(1:Iapnea(i)-1) temp(2:(end-1)) I.endi(Iapnea(i):end)];
    I.midi = [I.midi(1:Iapnea(i)) temp(2:(end-1))+di_TiTtot I.midi((Iapnea(i)+1):end)];
    VI = [VI(1:(Iapnea(i)-1)) VI(Iapnea(i))*ones(1,Nbreathsapnea(i)) VI((Iapnea(i)+1):end)];
    VTi = [VTi(1:(Iapnea(i)-1)) VTi(Iapnea(i))/Nbreathsapnea(i)*ones(1,Nbreathsapnea(i)) VTi((Iapnea(i)+1):end)];
end
Ti=(I.midi-I.starti)*dt; Te=(I.endi-I.starti)*dt; 


%% force ventilation non-zero
for i=1:length(VE)
    if VE(i)<0,
        VE(i)=0;
    end
    if VI(i)<0,
        VI(i)=0;
    end
end

%%
    figure(2);
    ax2(2)=subplot(2,1,1); plot(time,Vflow,time(I.starti),Vflow(I.starti),'.r',time(I.midi),Vflow(I.midi),'.k');



% end flow analysis


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
