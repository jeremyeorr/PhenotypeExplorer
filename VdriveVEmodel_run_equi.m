function [UAGSEM,UAGvalue,PcritSEM,Pcritvalue]=VdriveVEmodel_run_equi(Vdrive,VE,plotpoints,NEDoption,figurehandle)
global fixedslope fixedArthres 
%fixedslope=[];
figure(figurehandle);
if plotpoints
    %figurehandle=figure();
    %plot(Vdrive,VE,'.','markersize',16); 
    hold('on'); box('off');
end

%lsqoptions=optimset('display','off','maxiter',500,'tolx',10E-3);
 OPTIONS = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',200,'MaxFunEvals',200,'Algorithm','interior-point');
    
 
lower=[-1 -1.2]; %slope sleep curve; x=1,y=Vcrit, arthres
upper=[+2 1.2];
%parameters=[2 -0.5]; %starting parameters: the guess

modeloption=1;


% if 1
% %ymodeltest = pcrit_model([-0.2,1],[2:-0.01:1 1:0.01:2],modeloption,[2:-0.01:1 1:-0.005:0.5]);
% [restest,ymodeltest] = pcrit_model([-0.9,1],Vdrive,modeloption,VE);
% [restest,ymodeltest] = pcrit_model(parameters,Vdrive,modeloption,VE);
% figure(2)
% plot(Vdrive,[VE;ymodeltest],'.')
% end

parameters_start=zeros(20,2);
for i=1:30
    if i==1
        parameters_start(i,:)=[0 0];
    elseif i==2
        parameters_start(i,:)=[1 -1];
    elseif i==3
        parameters_start(i,:)=[1  1];
    elseif i==4
        parameters_start(i,:)=[-1 2];
    elseif i==5
        parameters_start(i,:)=[0.5 0.5];
    elseif i==6
        parameters_start(i,:)=[0.1 -1];
    elseif i==7
        parameters_start(i,:)=[0.99 0.99];    
    else
        parameters_start(i,:)=rand(1,2).*(upper-lower)+lower;
    end
    %[parametersi(i,:),Fres(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(@(parameters,Vdrive) pcrit_model(parameters,Vdrive,modeloption,VE),parameters_start,Vdrive,VE,lower,upper,lsqoptions);
    [parametersi(i,:),Fres(i),~,~] = fmincon(@(parameters) pcrit_model(parameters,Vdrive,modeloption,VE,NEDoption),parameters_start(i,:),[],[],[],[],lower,upper,[],OPTIONS); %[],[] = @(Parameters) model_nonlinear_constraints() ...
end
% 
% figure(200);
% plot(parameters_start(:,1),parameters_start(:,2),'.')

parametersi(isnan(Fres),:)=[]; Fres(isnan(Fres))=[];
[f,i]=min(Fres); parameters=parametersi(i,:);
disp(['SSE=',num2str(f)]);

noJacobian=0;
if ~noJacobian
    clear modelY_ modelY1 modelY2 f_ f1 f2
    dX=0.001;
    for i=1:length(parameters)
        parametersX1=parameters;
        parametersX2=parameters;
        parametersX1(i)=parameters(i)+dX;
        parametersX2(i)=parameters(i)-dX;
        [f_(i),modelY_]=pcrit_model(parameters,Vdrive,modeloption,VE,NEDoption);
        [f1(i),modelY1(:,i)]=pcrit_model(parametersX1,Vdrive,modeloption,VE,NEDoption);
        [f2(i),modelY2(:,i)]=pcrit_model(parametersX2,Vdrive,modeloption,VE,NEDoption);
        JACOBIAN = (modelY2-modelY1)/dX;
        RESIDUAL = (modelY_-VE)';
    end
end


Pmaskline=min(Vdrive):0.01:max(Vdrive);
%[modelVdotMax] = pcrit_model(parameters,Pmaskline,modeloption);
%confidence intervals

if ~noJacobian
[modelVdotMax,delta] = nlpredci(@(parameters,Pmaskline) pcrit_model2(parameters,Pmaskline,modeloption,VE),Pmaskline,parameters,RESIDUAL,'Jacobian',JACOBIAN);
if size(modelVdotMax,1)~=1 %to do -- work out why modelVdotMax is sometimes sideways...
    modelVdotMax=modelVdotMax';
end
if size(delta,1)~=1 %to do -- work out why modelVdotMax is sometimes sideways...
    delta=delta';
end
upperSEM=modelVdotMax+delta/1.96;
lowerSEM=modelVdotMax-delta/1.96;
if size(upperSEM,2)==1 %to do -- work out why modelVdotMax is sometimes sideways...
upperSEM=upperSEM';
lowerSEM=lowerSEM';
end

filly=[upperSEM fliplr(lowerSEM)];
fillx=[Pmaskline fliplr(Pmaskline)];
fill(fillx,filly,[0.8 0.2 0.2],'linestyle','none','facealpha',0.5);
end

if ~plotpoints
    hold('on'); box('off');
end

if noJacobian
    [f,modelVdotMax]=pcrit_model(parameters,Pmaskline,modeloption,NaN*Pmaskline,NEDoption);
end

plot(Pmaskline,modelVdotMax,'r','linewidth',2);

% if plotpoints
%     plot(Vdrive,VE,'k.','markersize',18);
% end

Pcritvalue=parameters(2);
UAGvalue=parameters(1);

if ~noJacobian
CI_parameters = nlparci(parameters,RESIDUAL,'jacobian',JACOBIAN);
PcritSEM=(CI_parameters(2,2)-Pcritvalue)/1.96;
UAGSEM=(CI_parameters(1,2)-UAGvalue)/1.96;
herr=errorbar(1,Pcritvalue,PcritSEM,'-','color',[0.7 0 0.7]);
else
  PcritSEM=NaN;
  UAGSEM=NaN;
end
% ArThresvalue=parameters(end);
% ArThresSEM=(CI_parameters(end,2)-ArThresvalue)/1.96;




% if Pcritvalue>0
%     I=find(abs(upperSEM)>0);
%     x=[0 Pmaskline(I)];
%     yupper=[Vcritvalue+VcritSEM upperSEM(I)];
%     ylower=[Vcritvalue-VcritSEM lowerSEM(I)];
%     Pmaskline2=0:0.01:Pmaskline(I);
%     upperSEMinterp=interp1(x,yupper,Pmaskline2,'spline');
%     lowerSEMinterp=interp1(x,ylower,Pmaskline2,'spline');
%     
%     plot([0 Pmaskline(I)],[Vcritvalue modelVdotMax(I)],'r:');
%     plot(Pmaskline2,lowerSEMinterp,'r:');
%     plot(Pmaskline2,upperSEMinterp,'r:');
% end

title(['Vpassive=' num2str(Pcritvalue,3) setstr(177) num2str(PcritSEM,3) '; ' 'UAG=' num2str(UAGvalue,2) setstr(177) num2str(UAGSEM,2) ] )
if 0
xlim([min([Vdrive,0])-0.2,max(Vdrive)]);
ylim([min([Vdrive,0])-0.2,max(Vdrive)]);
else
xlim([-1 5]);
ylim([-1 5]);
end





function [f,y] = pcrit_model(x,xdata,modeloption,ydata,NEDoption)
global fixedslope
if ~isempty(fixedslope)&&~isnan(fixedslope)
        y=fixedslope*(xdata-1)+x(1); %x(1) is Vcrit y=m(x-1)+c
        y(y>xdata)=xdata(y>xdata);
        y(y<0)=0;
else
        y=x(1)*(xdata-1)+x(2); %x(2) is Pcrit
        y(y>xdata)=xdata(y>xdata);
        y(y<0)=0;
end

minx=min(xdata); maxx=max(xdata);
xrange=minx:0.005:maxx;
if ~isempty(fixedslope)&&~isnan(fixedslope)
        yrange=fixedslope*(xrange-1)+x(1); %x(1) is Vcrit y=m(x-1)+c
        yrange(yrange>xrange)=xrange(yrange>xrange);
        yrange(yrange<0)=0;
else
        yrange=x(1)*(xrange-1)+x(2); %x(2) is Pcrit
        yrange(yrange>xrange)=xrange(yrange>xrange);
        yrange(yrange<0)=0;
end
if 0
   figure(99);plot(xrange,yrange,xdata,ydata,'.');
end

errorXYsq=zeros(1,length(xdata));
errorXYsqUnity=zeros(1,length(xdata));
for i=1:length(xdata)
   errorXYsq(i)=min((xdata(i)-xrange).^2+(ydata(i)-yrange).^2); %squared euclidian error
end

% errorsqY=(y-ydata).^2;
if NEDoption==1
    thres1=0;
elseif NEDoption==2
    thres1=0.33;
else
    thres1=-100;
end
    
if 1
    if x(1)>thres1 %penalise as normal
        error=errorXYsq;
    else %penalize data using distance to y=x if closer to unity than to NED model fit.
        for i=1:length(xdata)
            errorXYsqUnity(i)=min((xdata(i)-xrange).^2+(ydata(i)-xrange).^2); %squared euclidian distance to line of unity
        end
        error=min([errorXYsq;errorXYsqUnity]);
    end
else
    error=errorXYsq;
end

error(error>prctile(error,95))=[];
f = sum(error);
% 
% function [c,ceq] = cons()
% c=[];
% ceq = [];

function y = pcrit_model2(x,xdata,modeloption,ydata)
global fixedslope
if ~isempty(fixedslope)&&~isnan(fixedslope)
        y=fixedslope*(xdata-1)+x(1); %x(1) is Vcrit y=m(x-1)+c
        y(y>xdata)=xdata(y>xdata);
        y(y<0)=0;
else
        y=x(1)*(xdata-1)+x(2); %x(2) is Pcrit
        y(y>xdata)=xdata(y>xdata);
        y(y<0)=0;
end

