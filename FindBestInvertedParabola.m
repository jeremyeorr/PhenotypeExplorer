function [c,indices,y] = FindBestInvertedParabola(FlowData,Thres) 

%Thres = 75;
start=[max(FlowData)*2]; % A reasonable estimate of a starting physiological value for each parameter
lower=[0]; % Minimum conveivable physiological value for each parameter
upper=[max(FlowData)*20]; % Maximum conveivable physiological value for each parameter

    OPTIONS = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',4000,'MaxFunEvals',150,'Algorithm','interior-point');
    Parameters=start;
    [c,Error,~,~] = fmincon(@(Parameters) TheInvertedParabola(Parameters,FlowData,Thres),Parameters,[],[],[],[],lower,upper,@(Parameters) model_nonlinear_constraints(),OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]
    
    Ploton = 1;
    [~,y,Rsquared]=TheInvertedParabola(c,FlowData,Thres);
    F_all = sum(FlowData)/sum(y);
    mid50i = round(length(FlowData)*0.25):round(length(FlowData)*0.75);
    F_mid50 = sum(FlowData(mid50i))/sum(y(mid50i));
    F_max = 1-max(y-FlowData')/c;
    indices = [F_all F_mid50 F_max Rsquared];
end
%%
function [c,ceq] = model_nonlinear_constraints()
c=[];
ceq = [];
end