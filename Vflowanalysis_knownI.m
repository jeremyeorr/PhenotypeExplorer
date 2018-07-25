function [I,VI,VTi,Ti,Te,leak]=Vflowanalysis_knownI(Vflow,time,dt,minimum_figs,I)

%% Respiratory Trace Analyse
%[Pxy,F] = CPSD(X,Y,WINDOW,NOVERLAP,NFFT,Fs)
% Respiratory Trace Analyse
detectapneas=0;

%transpose if needed:
if size(Vflow,2)<size(Vflow,1)
    Vflow=Vflow';
end
if size(time,2)<size(time,1)
    time=time';
end

    leak = mean(Vflow)
    Vflow = Vflow-leak;

%Filter gently for use in analysis
filter_HFcutoff_butter1 = 20;
filter_order1 = 2;
[B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_HFcutoff_butter1]/(1/dt/2),'low');

VFlow_filtered1 = filtfilt(B_butterHcut,A_butterHcut,Vflow);

figure(1)
ax2(1)=subplot(1,1,1); plot(time,Vflow,time,VFlow_filtered1);


%% Find max flow values (i.e. expiration..)

BB_i_start=0*I.starti;

for i=1:length(I.starti)
    [~,index]=max(VFlow_filtered1(I.starti(i):I.midi(i)));
    BB_i_start(i)=I.starti(i)-1+index;
end

if ~minimum_figs
    figure(3); subplot(3,1,1); plot(time,VFlow_filtered1,time(BB_i_start),VFlow_filtered1(BB_i_start),'.');
end
%Store data
VFlow_crossing_step1_i=BB_i_start;


% start from peak data and move leftwards to find zero crossing:
for i=1:length(I.starti)
    if i>1
        searchrange=I.midi(i-1):BB_i_start(i);
    else
        searchrange=I.starti(i):BB_i_start(i);
    end
        temp=find(VFlow_filtered1(searchrange)<=0,1,'last');
    if isempty(temp)
        temp=I.starti(i);
    else
        temp=temp+searchrange(1);
    end
    BB_i_start(i)=temp;
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
BB_i_start = [BB_i_start;I.endi(end)];

clear BB_i_mid
for i=1:(length(BB_i_start)-1)
    [maxvalue,index]=max(cumsum(VFlow_filtered1(BB_i_start(i):(BB_i_start(i+1)-1))));
    BB_i_mid(i) = BB_i_start(i)+index;
end

if ~minimum_figs
    figure(3);  subplot(3,1,3); plot(time,VFlow_filtered1,time(BB_i_start),VFlow_filtered1(BB_i_start),'.r',time(BB_i_mid),VFlow_filtered1(BB_i_mid),'.k')
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

%% If breaths are too small, use known I values instead

VTi_thres = mean(VTi)/20;
VTe_thres = mean(VTe)/20;

for i=1:length(BB_i_start)
    if VTe(i)<VTe_thres||VTi(i)<VTi_thres
        BB_i_start(i)=I.starti(i);
        BB_i_mid(i)=I.midi(i);
    end
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


%% Start/End breath indices in original indices
clear I;
I.starti = BB_i_start;
I.midi = BB_i_mid;
I.endi = BB_i_end;

%%
for i=1:length(VE) %force ventilation non-zero.
    if VE(i)<0,
        VE(i)=0;
    end
    if VI(i)<0,
        VI(i)=0;
    end
end
if ~minimum_figs
    figure(2);
    ax2(2)=subplot(2,1,1); plot(time,VFlow_filtered1,time(BB_i_start),VFlow_filtered1(BB_i_start),'.r',time(BB_i_mid),VFlow_filtered1(BB_i_mid),'.k');
end


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
