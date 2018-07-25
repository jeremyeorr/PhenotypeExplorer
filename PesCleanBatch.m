%read, decontaminate Pes and save as Pes_clean.values
for m=[10 3:9 11:12 14:29]
%%
m
filenameanddir=['J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\SpikeFiles\' subj{m}]
clear  Pes EKG PesToFlow FlowEdi Flow
load(filenameanddir, 'Pes', 'EKG','PesToFlow','FlowEdi','Flow');
dt = Flow.interval;
Time = (0:dt:(dt*(length(Flow.values)-1)))';

channels={'Pes','EKG','FlowEdi','Flow'};
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

Pes_clean.values = PesRemoveEKG(Pes.values,EKG.values,dt);
Pes_clean.interval = dt;

%%
%% Pes To Flow: Filtered Pes (gentle low pass)
filter_HFcutoff_butter0 = 8;
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');

%Pes(isnan(Pes))=nanmean(Pes);
Pes1_filtered = filtfilt(B_butter0,A_butter0,Pes_clean.values);

%% Estimate respiratory impedance (Pes)
rsquared_x=[];
parameters_x=[];
for x=1:size(PesToFlow.ranges,1)
    try
        %get data
        I=find(Time>PesToFlow.ranges(x,1)&Time<PesToFlow.ranges(x,2));
        
        xdata_filtered = Pes_clean.values(I);
        
        minimum_figs=1;
        [I_edi,~,~,~,~]=Vflowanalysis2(FlowEdi.values(I),Time(I),dt,minimum_figs);
        baselinePes = median(xdata_filtered(I_edi.starti));
        
        downsamplefactor=1;
        xdata.data = downsample(-xdata_filtered+baselinePes,downsamplefactor); %downsample(Pes(I),downsamplefactor);
        xdata.time = downsample(Time(I),downsamplefactor);
        xdata.dt = dt;
        ydata = downsample(Flow.values(I),downsamplefactor);
        
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
        
        figure(77);
        ax(1)=subplot(3,2,1);
        plot(xdata.time,ydata,'k');
        hold('on');
        plot(xdata.time,modelY,'r');
        hold('off');
        ax(2)=subplot(3,2,3);
        plot(xdata.time,cumsum(ydata-parameters(5))*dt,'k');
        hold('on');
        plot(xdata.time,cumsum(modelY-parameters(5))*dt,'r');
        hold('off');
        ax(3)=subplot(3,2,5);
        plot(xdata.time,xdata.data,'k');
        linkaxes(ax,'x');
        
        subplot(2,2,2);
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
parameters_x = [parameters_x;parameters_median];
%% Test: Use parameters from X on segment Y (Pes)
PesToFlow.Parameters = parameters_x;
PesToFlow.Rsquared=[];

x=size(parameters_x,1);
for i=1:size(parameters_x,1)-1
    try
        y=i;
        
        I=find(Time>PesToFlow.ranges(y,1)&Time<PesToFlow.ranges(y,2));
        
        xdata_filtered = Pes_clean.values(I);
        
        minimum_figs=1;
        [I_edi,~,~,~,~]=Vflowanalysis2(FlowEdi.values(I),Time(I),dt,minimum_figs);
        baselinePes = median(xdata_filtered(I_edi.starti));
        
        downsamplefactor=1;
        xdata.data = downsample(-xdata_filtered+baselinePes,downsamplefactor); %downsample(Pes(I),downsamplefactor);
        xdata.time = downsample(Time(I),downsamplefactor);
        xdata.dt = dt;
        
        ydata = downsample(Flow.values(I),downsamplefactor);
        
        modelY = PesToVflowModel(parameters_x(x,:),xdata); %is this use of 'x' right?
        
        SSE = sum(((modelY-mean(modelY))-(ydata-mean(ydata))).^2);
        SStot = sum((modelY-mean(modelY)).^2);
        
        figure(77)
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
parametersOSA=parameters_x(PesToFlow.useparameterset,:);
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

    save(filenameanddir,'FlowPes','PesToFlow','Pes_clean','-append');

%%
end
