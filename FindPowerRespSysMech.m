function [meanRsquared] = FindPowerRespSysMech(EdiToFlow,Powers,Edi,Flow,Time,dt)

for j=1:length(Powers)
parameters_x=[];
j
for x=1:size(EdiToFlow.ranges,1)
    
    %get data
    I=find(Time>EdiToFlow.ranges(x,1)&Time<EdiToFlow.ranges(x,2));
    
    filter_HFcutoff_butter0 = 8;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    tempEdi = Edi(I);
    tempEdi(tempEdi<0) = 0;
    xdata_filtered = filtfilt(B_butter0,A_butter0,tempEdi.^Powers(j));
    
    downsamplefactor=1;
    xdata.data = downsample(xdata_filtered,downsamplefactor); %downsample(Pes(I),downsamplefactor);
    xdata.time = downsample(Time(I),downsamplefactor);
    xdata.dt = dt;
    ydata = downsample(Flow(I),downsamplefactor);
    
    
    lsqoptions=optimset('display','off','MaxFunEvals',11);
    parameters=[6 10 0 0 0]; %guess: R E VL(1)-FRC alpha Vleak delX pwr %alpha = Younes exponent for length-tension, 0 is no effect    
    lower=[0.5 0.1 -1  -0.1 -60]; %R E VL(1)-FRC alpha Vleak delX pwr
    upper=[60 60 1 0.5 +60]; %R C VL(1)-FRC alpha Vleak delX  pwr
    
    [modelYguess] = EdiToVflowModel(parameters,xdata);
    
    for i=1:5
        [parameters,resnormC1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@EdiToVflowModel,parameters,xdata,ydata,lower,upper,lsqoptions);
    end
    
    [modelY] = EdiToVflowModel(parameters,xdata);
    
    rsquared_x(x) = sum((modelY-ydata).^2)/sum((ydata-mean(ydata)).^2);
    if 0
    figure(111)
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
    plot(cumsum(ydata-parameters(5))*dt,cumsum(modelY-parameters(5))*dt,'k');
    end
    parameters_x(x,:)=real(parameters);
    
end

parameters_median = median(parameters_x);
parameters_x = [parameters_x;parameters_median];
close all


%% Test: Use parameters from X on segment Y
%EdiToFlow.ParametersPower75 = parameters_x;
Rsquared=[];

x=size(parameters_x,1);
for i=1:size(parameters_x,1)-1
y=i;

I=find(Time>EdiToFlow.ranges(y,1)&Time<EdiToFlow.ranges(y,2));

filter_HFcutoff_butter0 = 8;
filter_order0 = 2;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');

tempEdi = Edi(I);
    tempEdi(tempEdi<0) = 0;
    xdata_filtered = filtfilt(B_butter0,A_butter0,tempEdi.^Powers(j));
    
downsamplefactor=1;
xdata.data = downsample(xdata_filtered,downsamplefactor); %downsample(Pes(I),downsamplefactor);
xdata.time = downsample(Time(I),downsamplefactor);
xdata.dt = dt;

ydata = downsample(Flow(I),downsamplefactor);

modelY = EdiToVflowModel(parameters_x(x,:),xdata);

SSE = sum(((modelY-mean(modelY))-(ydata-mean(ydata))).^2);
SStot = sum((modelY-mean(modelY)).^2);

if 0
figure(111)
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
plot(cumsum(ydata-parameters(5))*dt,cumsum(modelY-parameters(5))*dt,'k');
end
Rsquared(i) = 1-SSE/SStot;
end

meanRsquared(j)=mean(Rsquared);
close all
end


