function [VcritSEM,Vcritvalue,PcritSEM,Pcritvalue,PVSlope,Rsquared]=pcrit_model_run(Pmask,Vflow,exclude_zeroflowdata,plotpoints)
global fixedslope upperlimit

nandata = isnan(Pmask)|isnan(Vflow);
Pmask(nandata)=[];
Vflow(nandata)=[];
range=max(Pmask)-min(Pmask);

if exclude_zeroflowdata
Pmask(Vflow==0)=[];
Vflow(Vflow==0)=[];
end

%upper limit is under construction.
% if 1
%     upperlimit.on=1;
% else
%     upperlimit.on=0;
% end

if plotpoints==1
    plot(Pmask,Vflow,'.','markersize',16); hold('on'); box('off');
elseif plotpoints==2
    plot(Pmask,Vflow,'r.','markersize',16); hold('on'); box('off');
end

lsqoptions=optimset('display','off','maxiter',500,'tolx',10E-3);
Pcrit_lowerlimit=min(Pmask)-range*4;
Pcrit_upperlimit=max(Pmask);
lower=[0 Pcrit_lowerlimit]
upper=[max(Vflow)*10 Pcrit_upperlimit]
parameters=[1 Pcrit_lowerlimit];

if upperlimit.on==1
    lower(3)=min(Pmask);
    upper(3)=max(Pmask);
    parameters(3)=50;
end
    
if ~isempty(fixedslope)&&~isnan(fixedslope)
    lower(1)=[];
    upper(1)=[];
    parameters(1)=[];
end

modeloption=1; %not used just yet...
if ~isempty(fixedslope)&&~isnan(fixedslope)
    Ytemp=pcrit_model_fixedslope(parameters,Pmask,modeloption);
else
    Ytemp=pcrit_model(parameters,Pmask,modeloption);
end
%plot(Pmask,Ytemp,'r.');
clear Ytemp


for i=1:5
    if ~isempty(fixedslope)&&~isnan(fixedslope)
        [parameters,~,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,Pmask) pcrit_model_fixedslope(parameters,Pmask,modeloption),parameters,Pmask,Vflow,lower,upper,lsqoptions);
    else
        [parameters,~,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,Pmask) pcrit_model(parameters,Pmask,modeloption),parameters,Pmask,Vflow,lower,upper,lsqoptions);
    end
end
%have a good start point - but can we do better:
if ~isempty(fixedslope)&&~isnan(fixedslope)
lower=[parameters(1)-range/10];
upper=[parameters(1)+range/10];
else
lower=[parameters(1)/3 parameters(2)-range/10];
upper=[parameters(1)*2 parameters(2)+range/10];    
end

if upperlimit.on==1
    lower(end+1)=min(Pmask);
    upper(end+1)=max(Pmask);
end

for i=1:50
    parameters_start=randn*(upper-lower)+lower;
    if ~isempty(fixedslope)&&~isnan(fixedslope)
        [parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,Pmask) pcrit_model_fixedslope(parameters,Pmask,modeloption),parameters,Pmask,Vflow,lower,upper,lsqoptions);
    else
        [parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,Pmask) pcrit_model(parameters,Pmask,modeloption),parameters_start,Pmask,Vflow,lower,upper,lsqoptions);
    end
end
parametersi(isnan(Fres),:)=[]; Fres(isnan(Fres))=[];
[~,i]=min(Fres); parameters=parametersi(i,:);

%have a good better point - but can we do even better again:

if ~isempty(fixedslope)&&~isnan(fixedslope)
lower=[parameters(1)-range/50];
upper=[parameters(1)+range/50];
else
lower=[parameters(1)/1.5 parameters(2)-range/50];
upper=[parameters(1)*1.5 parameters(2)+range/50];
end
if upperlimit.on==1
    lower(end+1)=min(Pmask);
    upper(end+1)=max(Pmask);
end

for i=1:100
    parameters_start=randn*(upper-lower)+lower;
    if ~isempty(fixedslope)&&~isnan(fixedslope)
        [parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,Pmask) pcrit_model_fixedslope(parameters,Pmask,modeloption),parameters,Pmask,Vflow,lower,upper,lsqoptions);
    else
        [parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,Pmask) pcrit_model(parameters,Pmask,modeloption),parameters_start,Pmask,Vflow,lower,upper,lsqoptions);
    end
end
parametersi(isnan(Fres),:)=[]; Fres(isnan(Fres))=[];
[minFres,i]=min(Fres); parameters=parametersi(i,:);
%rsquared
Ftot=sum((Vflow-mean(Vflow)).^2);
Rsquared=1-minFres/Ftot;

Pmaskline=min(Pmask):0.01:max(Pmask);
%[modelVdotMax] = pcrit_model(parameters,Pmaskline,modeloption);
%confidence intervals
if ~isempty(fixedslope)&&~isnan(fixedslope)
    [modelVdotMax,delta] = nlpredci(@(parameters,Pmaskline) pcrit_model_fixedslope(parameters,Pmaskline,modeloption),Pmaskline,parameters,RESIDUAL,'Jacobian',JACOBIAN);
else
    [modelVdotMax,delta] = nlpredci(@(parameters,Pmaskline) pcrit_model(parameters,Pmaskline,modeloption),Pmaskline,parameters,RESIDUAL,'Jacobian',JACOBIAN);
end

upperSEM=modelVdotMax+delta/1.96;
lowerSEM=modelVdotMax-delta/1.96;

filly=[upperSEM fliplr(lowerSEM)];
fillx=[Pmaskline fliplr(Pmaskline)];

fill(fillx,filly,[0.8 0.2 0.2],'linestyle','none','facealpha',0.5); 
if plotpoints==0
hold('on'); box('off');
end

plot(Pmaskline,modelVdotMax,'r'); 

if plotpoints==1
    plot(Pmask,Vflow,'k.','markersize',18); 
elseif plotpoints==2
    plot(Pmask,Vflow,'r.','markersize',18); 
end

CI_parameters = nlparci(parameters,RESIDUAL,'jacobian',JACOBIAN);
Pcritvalue=parameters(end-upperlimit.on);
PcritSEM=(CI_parameters(end,end)-Pcritvalue)/1.96;
%ci = nlparci(beta,resid,'jacobian',J)

title(['Pcrit=' num2str(Pcritvalue,4) setstr(177) num2str(PcritSEM,3)]);

if ~isempty(fixedslope)&&~isnan(fixedslope)
    PVSlope=fixedslope;
else
    PVSlope=parameters(1);
end

Vcritvalue=-PVSlope*Pcritvalue;

Vcrit_lowerlimit=Vcritvalue-1;
Vcrit_upperlimit=Vcritvalue+1;

if upperlimit.on
upperbreakpoint=parameters(end)
else
upperbreakpoint=[];
end

parameters2=[PVSlope Vcritvalue upperbreakpoint];
modeloption=2; %not used just yet...
%pcrit_model(parameters2,Pmask,2)

if ~isempty(fixedslope)&&~isnan(fixedslope)
    pcrit_model_fixedslope(parameters2,Pmask,modeloption);
else
    pcrit_model(parameters2,Pmask,modeloption);
end

if ~isempty(fixedslope)&&~isnan(fixedslope)
lower=Vcrit_lowerlimit;
upper=Vcrit_upperlimit;
else
lower=[PVSlope/1.5 Vcrit_lowerlimit];
upper=[PVSlope*1.5 Vcrit_upperlimit];
end
if upperlimit.on==1
    lower(end+1)=min(Pmask);
    upper(end+1)=max(Pmask);
end

for i=1:50
    if ~isempty(fixedslope)&&~isnan(fixedslope)
        [parameters2,~,RESIDUAL2,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN2] = lsqcurvefit(@(parameters,Pmask) pcrit_model_fixedslope(parameters,Pmask,modeloption),parameters,Pmask,Vflow,lower,upper,lsqoptions);
    else
        [parameters2,~,RESIDUAL2,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN2] = lsqcurvefit(@(parameters,Pmask) pcrit_model(parameters,Pmask,modeloption),parameters,Pmask,Vflow,lower,upper,lsqoptions);
    end
end

if ~isempty(fixedslope)&&~isnan(fixedslope)
    [modelVdotMax2,delta2] = nlpredci(@(parameters2,Pmaskline) pcrit_model_fixedslope(parameters2,Pmaskline,modeloption),Pmaskline,parameters2,RESIDUAL2,'Jacobian',JACOBIAN2);
else
    [modelVdotMax2,delta2] = nlpredci(@(parameters2,Pmaskline) pcrit_model(parameters2,Pmaskline,modeloption),Pmaskline,parameters2,RESIDUAL2,'Jacobian',JACOBIAN2);
end
    
if 0 %check that models are identical
filly=[modelVdotMax2+delta2' fliplr(modelVdotMax2-delta2')]
fillx=[Pmaskline fliplr(Pmaskline)]
fill(fillx,filly,[0.2 0.8 0.2],'linestyle','none','facealpha',0.5); 
plot(Pmaskline,modelVdotMax2,'b:'); 
end

% if ~isempty(fixedslope)&&~isnan(fixedslope)
%     parameters2(1)=fixedslope;
% end

plot(0,Vcritvalue,'ro'); 

CI_parameters2 = nlparci(parameters2,RESIDUAL2,'jacobian',JACOBIAN2)

VcritSEM=abs((CI_parameters2(end-upperlimit.on,end)-parameters2(end-upperlimit.on))/1.96);

if Pcritvalue>0
I=find(abs(upperSEM)>0);
x=[0 Pmaskline(I)];
yupper=[Vcritvalue+VcritSEM upperSEM(I)];
ylower=[Vcritvalue-VcritSEM lowerSEM(I)];
Pmaskline2=0:0.01:Pmaskline(I);
upperSEMinterp=interp1(x,yupper,Pmaskline2,'spline')
lowerSEMinterp=interp1(x,ylower,Pmaskline2,'spline')

plot([0 Pmaskline(I)],[Vcritvalue modelVdotMax(I)],'r:'); 
plot(Pmaskline2,lowerSEMinterp,'r:'); 
plot(Pmaskline2,upperSEMinterp,'r:'); 
end

title(['Pcrit=' num2str(parameters(end-upperlimit.on),4) setstr(177) num2str(PcritSEM,3) '; ' 'Vcrit=' num2str(parameters2(end-upperlimit.on),4) setstr(177) num2str(VcritSEM,3) ] )
xlim([min([Pmask,0])-0.2,max(Pmask)])

herr=errorbar(0,Vcritvalue,VcritSEM,'r-');


function y = pcrit_model(x,xdata,modeloption)
% global fixedslope
% if ~isempty(fixedslope)&&~isnan(fixedslope)
%    x(1)=fixedslope;
% end
global upperlimit
if modeloption==1
y=x(1)*(xdata-x(2)); %x(2) is Pcrit
y(y<0)=0;
end
if modeloption==2
y=x(1)*(xdata)+x(2); %x(2) is Vcrit
y(y<0)=0;
end
if upperlimit.on
    if modeloption==1
        temp = x(1)*(x(3)-x(2));
    elseif modeloption==2
        temp = x(1)*(x(3)+x(2));
    end
    y(xdata>x(3))=temp; %upperlimit
end

function y = pcrit_model_fixedslope(x,xdata,modeloption)
global fixedslope upperlimit
if modeloption==1
y=fixedslope*(xdata-x(1)); %x(1) is Pcrit
y(y<0)=0;
end
if modeloption==2
y=fixedslope*(xdata)+x(1); %x(1) is Vcrit
y(y<0)=0;
end
if upperlimit.on
    if modeloption==1
        temp = fixedslope*(x(2)-x(1));
    elseif modeloption==2
        temp = fixedslope*(x(2)+x(1));
    end
    y(xdata>x(2))=temp; %upperlimit
end
