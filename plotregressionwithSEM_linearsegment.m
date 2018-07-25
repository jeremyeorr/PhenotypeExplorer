function [Rvalue,Pvalue,PVSlope,Pvalue2,VcritSEM,Vcritvalue,PcritSEM,Pcritvalue,parameters,Xmodel_sem]=plotregressionwithSEM_linearsegment(xdata,ydata,exclude_zeroflowdata,plotpoints,varargin)
global fixedslope upperlimit
%fixedslope=NaN;
%upperlimit.on=1;


if exclude_zeroflowdata
xdata(ydata==0)=[];
ydata(ydata==0)=[];
end



if length(varargin)>0&&strcmp(varargin(1),'lowerlimitonupperbreakpoint')
    lowerlimitonupperbreakpoint=varargin{2};
else
    lowerlimitonupperbreakpoint=NaN;
end

nandata = isnan(xdata)|isnan(ydata);
xdata(nandata)=[];
ydata(nandata)=[];
range=max(xdata)-min(xdata);

if plotpoints==1
    plot(xdata,ydata,'.','markersize',16); hold('on'); box('off');
elseif plotpoints==2
    plot(xdata,ydata,'r.','markersize',16); hold('on'); box('off');
end

lsqoptions=optimset('display','off','maxiter',500,'tolx',10E-3);
Pcrit_lowerlimit=min(xdata)-range*4;
Pcrit_upperlimit=max(xdata);
lower=[0 Pcrit_lowerlimit];
upper=[max(ydata)*10 Pcrit_upperlimit]; %max(ydata)*10
parameters=[0.1 Pcrit_lowerlimit];

if upperlimit.on==1
    if ~isnan(lowerlimitonupperbreakpoint)
        lower(3)=lowerlimitonupperbreakpoint;
    else
        lower(3)=min(xdata);
    end
    upper(3)=max(xdata);
    parameters(3)=max(xdata);
end
    
if ~isempty(fixedslope)&&~isnan(fixedslope)
    lower(1)=[];
    upper(1)=[];
    parameters(1)=[];
end

modeloption=1; %not used just yet...
if ~isempty(fixedslope)&&~isnan(fixedslope)
    Ytemp=pcrit_model_fixedslope(parameters,xdata,modeloption);
else
    Ytemp=pcrit_model(parameters,xdata,modeloption);
end
%plot(xdata,Ytemp,'r.');
clear Ytemp

for i=1:3 %overwrites results
    if ~isempty(fixedslope)&&~isnan(fixedslope)
        [parameters,Fres,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,xdata) pcrit_model_fixedslope(parameters,xdata,modeloption),parameters,xdata,ydata,lower,upper,lsqoptions);
    else
        [parameters,Fres,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,xdata) pcrit_model(parameters,xdata,modeloption),parameters,xdata,ydata,lower,upper,lsqoptions);
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
    if ~isnan(lowerlimitonupperbreakpoint)
        lower(end+1)=lowerlimitonupperbreakpoint;
    else
        lower(end+1)=min(xdata);
    end
    upper(end+1)=max(xdata);
end

parametersi(1,:)=parameters; 
for i=2:10
    parameters_start=randn*(upper-lower)+lower;
    if ~isempty(fixedslope)&&~isnan(fixedslope)
        [parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,xdata) pcrit_model_fixedslope(parameters,xdata,modeloption),parameters,xdata,ydata,lower,upper,lsqoptions);
    else
        [parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,xdata) pcrit_model(parameters,xdata,modeloption),parameters_start,xdata,ydata,lower,upper,lsqoptions);
    end
end
parametersi(isnan(Fres),:)=[]; Fres(isnan(Fres))=[];
[~,i]=min(Fres); parameters=parametersi(i,:);

%have a better point - but can we do even better again:

if ~isempty(fixedslope)&&~isnan(fixedslope)
lower=[parameters(1)-range/50];
upper=[parameters(1)+range/50];
else
lower=[parameters(1)/1.5 parameters(2)-range/50];
upper=[parameters(1)*1.5 parameters(2)+range/50];
end
if upperlimit.on==1
    if ~isnan(lowerlimitonupperbreakpoint)
        lower(end+1)=lowerlimitonupperbreakpoint;
    else
        lower(end+1)=min(xdata);
    end
    upper(end+1)=max(xdata);
end

tempI=size(parametersi,1);
for i=(tempI+1):(tempI+1)+10
    parameters_start=randn*(upper-lower)+lower;
    if ~isempty(fixedslope)&&~isnan(fixedslope)
        [parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,xdata) pcrit_model_fixedslope(parameters,xdata,modeloption),parameters,xdata,ydata,lower,upper,lsqoptions);
    else
        [parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,xdata) pcrit_model(parameters,xdata,modeloption),parameters_start,xdata,ydata,lower,upper,lsqoptions);
    end
end
parametersi(isnan(Fres),:)=[]; Fres(isnan(Fres))=[];
[minFres,i]=min(Fres); parameters=parametersi(i,:);
%rsquared
Ftot=sum((ydata-mean(ydata)).^2);
Rvalue=(1-minFres/Ftot)^0.5;


%get pvalue (adapted code from matlab's "NonLinearModel" function)
ssr = max(Ftot - minFres,0);
dfr = length(parameters) - 1;
dfe = length(xdata) - 1 - dfr;
f = (ssr./dfr)/(minFres/dfe);
Pvalue = fcdf(1./f,dfe,dfr); % upper tail


%get pvalue versus linear model
fun2 = @(x,xdata)(x(1)*xdata+x(2));
lower=[-Inf -Inf];
upper=[Inf Inf];
x02=polyfit(xdata,ydata,1);
[Xmodel2,Fres2,~,~,~,~,~]=lsqcurvefit(fun2,x02,xdata,ydata,lower,upper,lsqoptions);
ssr = max(Fres2 - minFres,0);
dfr = length(parameters) - length(x02);
dfe = length(xdata) - length(parameters);
f = (ssr./dfr)/(minFres/dfe);
Pvalue2 = fcdf(1./f,dfe,dfr); % upper tail


Xmodel_ci = nlparci(parameters,RESIDUAL,'jacobian',JACOBIAN);
Xmodel_sem = (Xmodel_ci(:,2)'-parameters)/1.96;


Pmaskline=min(xdata):0.01:max(xdata);
%[modelVdotMax] = pcrit_model(parameters,Pmaskline,modeloption);
%confidence intervals
if 0 %failed attempt at removing error: "SVD does not support sparse matrices"
    if upperlimit.on&&parameters(end)>=(max(xdata)-0.00001*range)
        parameters(end)=[];
        JACOBIAN(:,end)=[];
        upperlimit.on=0;
    end
end

if ~isempty(fixedslope)&&~isnan(fixedslope)
    predopt1='observation'; %'observation' 'curve'
    [modelVdotMax,delta] = nlpredci(@(parameters,Pmaskline) pcrit_model_fixedslope(parameters,Pmaskline,modeloption),Pmaskline,parameters,RESIDUAL,'Jacobian',JACOBIAN,'predopt',predopt1);
    predopt1='curve'; %'observation' 'curve'
    [modelVdotMax_,delta_] = nlpredci(@(parameters,Pmaskline) pcrit_model_fixedslope(parameters,Pmaskline,modeloption),Pmaskline,parameters,RESIDUAL,'Jacobian',JACOBIAN,'predopt',predopt1);
else
    predopt1='observation'; %'observation' 'curve'
    [modelVdotMax,delta] = nlpredci(@(parameters,Pmaskline) pcrit_model(parameters,Pmaskline,modeloption),Pmaskline,parameters,RESIDUAL,'Jacobian',JACOBIAN,'predopt',predopt1);
    predopt1='curve'; %'observation' 'curve'
    [modelVdotMax_,delta_] = nlpredci(@(parameters,Pmaskline) pcrit_model(parameters,Pmaskline,modeloption),Pmaskline,parameters,RESIDUAL,'Jacobian',JACOBIAN,'predopt',predopt1);    
end

% if 0&&upperlimit.on %kills the error bars on segment 3
% delta(Pmaskline>parameters(end))=0;
% end
if 1
    delta(modelVdotMax<=0)=0;
    delta_(modelVdotMax<=0)=0;
end

upperSEM=modelVdotMax+delta/1.96;
lowerSEM=modelVdotMax-delta/1.96;
filly=[upperSEM fliplr(lowerSEM)];
fillx=[Pmaskline fliplr(Pmaskline)];
fill(fillx,filly,[0.95 0.95 0.9875],'linestyle','none'); 

upperSEM_=modelVdotMax_+delta_/1.96;
lowerSEM_=modelVdotMax_-delta_/1.96;
filly=[upperSEM_ fliplr(lowerSEM_)];
fillx=[Pmaskline fliplr(Pmaskline)];
fill(fillx,filly,[0.8 0.8 0.95],'linestyle','none'); 


plot(Pmaskline,modelVdotMax,'k'); 
plot(xdata,ydata,'k.','markersize',16); 

CI_parameters = nlparci(parameters,RESIDUAL,'jacobian',JACOBIAN);
Pcritvalue=parameters(end-upperlimit.on);
PcritSEM=(CI_parameters(end,end)-Pcritvalue)/1.96;
%ci = nlparci(beta,resid,'jacobian',J)

%title(['Pcrit=' num2str(Pcritvalue,4) setstr(177) num2str(PcritSEM,3)]);

if ~isempty(fixedslope)&&~isnan(fixedslope)
    PVSlope=fixedslope;
else
    PVSlope=parameters(1);
end

Vcritvalue=-PVSlope*Pcritvalue;

Vcrit_lowerlimit=Vcritvalue-1;
Vcrit_upperlimit=Vcritvalue+1;

if upperlimit.on
upperbreakpoint=parameters(end);
else
upperbreakpoint=[];
end

parameters2=[PVSlope Vcritvalue upperbreakpoint];
modeloption=2; %not used just yet...
%pcrit_model(parameters2,xdata,2)

if ~isempty(fixedslope)&&~isnan(fixedslope)
    pcrit_model_fixedslope(parameters2,xdata,modeloption);
else
    pcrit_model(parameters2,xdata,modeloption);
end

if ~isempty(fixedslope)&&~isnan(fixedslope)
lower=Vcrit_lowerlimit;
upper=Vcrit_upperlimit;
else
lower=[PVSlope/1.5 Vcrit_lowerlimit];
upper=[PVSlope*1.5 Vcrit_upperlimit];
end
if upperlimit.on==1
    if ~isnan(lowerlimitonupperbreakpoint)
        lower(end+1)=lowerlimitonupperbreakpoint;
    else
        lower(end+1)=min(xdata);
    end
    upper(end+1)=max(xdata);
end

for i=1:1
    if ~isempty(fixedslope)&&~isnan(fixedslope)
        [parameters2,~,RESIDUAL2,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN2] = lsqcurvefit(@(parameters,xdata) pcrit_model_fixedslope(parameters,xdata,modeloption),parameters2,xdata,ydata,lower,upper,lsqoptions);
    else
        [parameters2,~,RESIDUAL2,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN2] = lsqcurvefit(@(parameters,xdata) pcrit_model(parameters,xdata,modeloption),parameters2,xdata,ydata,lower,upper,lsqoptions);
    end
end

try
if ~isempty(fixedslope)&&~isnan(fixedslope)
    [modelVdotMax2,delta2] = nlpredci(@(parameters2,Pmaskline) pcrit_model_fixedslope(parameters2,Pmaskline,modeloption),Pmaskline,parameters2,RESIDUAL2,'Jacobian',JACOBIAN2);
else
    [modelVdotMax2,delta2] = nlpredci(@(parameters2,Pmaskline) pcrit_model(parameters2,Pmaskline,modeloption),Pmaskline,parameters2,RESIDUAL2,'Jacobian',JACOBIAN2);
end
    
if 0 %check that models are identical
filly=[modelVdotMax2+delta2' fliplr(modelVdotMax2-delta2')]
fillx=[Pmaskline fliplr(Pmaskline)]
fill(fillx,filly,[0.8 0.8 0.95],'linestyle','none'); 
plot(Pmaskline,modelVdotMax2,'b:'); 
end

% if ~isempty(fixedslope)&&~isnan(fixedslope)
%     parameters2(1)=fixedslope;
% end

%plot(0,Vcritvalue,'ro'); 

catch me
      
end


CI_parameters2 = nlparci(parameters2,RESIDUAL2,'jacobian',JACOBIAN2);
VcritSEM=abs((CI_parameters2(end-upperlimit.on,end)-parameters2(end-upperlimit.on))/1.96);
if isnan(VcritSEM)
    Pmaskline2=0;
    [modelY,deltaY] = nlpredci(@(parameters,Pmaskline) pcrit_model_fixedslope(parameters,Pmaskline,modeloption),Pmaskline2,parameters,RESIDUAL,'Jacobian',JACOBIAN,'predopt',predopt1);
    VcritSEM=deltaY/1.96;
end

if Pcritvalue>0
I=find(abs(upperSEM)>0);
x=[0 Pmaskline(I)];
yupper=[Vcritvalue+VcritSEM upperSEM(I)];
ylower=[Vcritvalue-VcritSEM lowerSEM(I)];
Pmaskline2=0:0.01:Pmaskline(I);
upperSEMinterp=interp1(x,yupper,Pmaskline2,'spline')
lowerSEMinterp=interp1(x,ylower,Pmaskline2,'spline')

%plot([0 Pmaskline(I)],[Vcritvalue modelVdotMax(I)],'r:'); 
%plot(Pmaskline2,lowerSEMinterp,'r:'); 
%plot(Pmaskline2,upperSEMinterp,'r:'); 
end





title(['Pcrit=' num2str(parameters(end-upperlimit.on),4) setstr(177) num2str(PcritSEM,3) '; ' 'Vcrit=' num2str(parameters2(end-upperlimit.on),4) setstr(177) num2str(VcritSEM,3) ] )
%xlim([min([xdata,0])-0.2,max(xdata)])

%herr=errorbar(0,Vcritvalue,VcritSEM,'r-');

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
