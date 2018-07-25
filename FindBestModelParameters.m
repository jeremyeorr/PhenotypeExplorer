function [Error,Vchem,Varousal,LoopGainParameters,BestParameters,BestSSres,FitQuality,i_1,i_end] = FindBestModelParameters(Data,VraOn,Veupnea,polyfitorder) 
% Written by Philip I Terrill and Scott A Sands
% This MATLAB function accompanies the manuscript entitled:
% "Quantifying the Ventilatory Control Contribution to Obstructive Sleep Apnea"
% **************************************************************************
% **************************************************************************
% This function is used to identify the set of parameters (LG_0, Tau,
% delay, Vdr[1]) that results in the best fit between the first order
% ventilation model and the observed ventilation. The MATLAB function
% fmincon.m is used to search for the optimal parameter set. 
% >>>>>>>>>>>>>>>>>>>>>>>>Input parameters<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% ***Data: This matrix contains information about the patients observed
% ventilation
%             - Column 1: Breath-by-Breath Minute Ventilation -- Must already be "mean subtracted" i.e. eupnea=0. 
%             - Column 2: Obstructive respiratory event series. Elements=1 for an un-obstructed breath; Elements=0 for an obstructed breath.
%             - Column 3: Arousal Scoring series: defines whether there is an arousal present in each breath (arousal present, AR[n]=1; No arousal present, AR[n]=0);
%             - Column 4: The time that each breath in the series occurs. 
%             - Column 5: Scored Central Apnea array series: defines whether each breath is scored as a central apnea. elements=1 for breaths in a scored central apnea; Elements=0 for breaths not classified as a central apnea.
%             - Column 6: T_tot series (currently unused)
%             - Column 7: T_tot series for the previous breath. 
% ***VraOn: This binary variable allows the inclusion or exclusion of the
% ventilatory response to arousal parameter.  
% ***Veupnea: This is the eupneic ventilation level. This is typically normalised to 1. 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%
% >>>>>>>>>>>>>>>>>>>>>>>>Output parameters<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% - Error: 1*N Array. This is the breath-by-breath different between the
% model predicted ventilation and the observed ventilation
% - Vchem: 1*N Array. The ventilatory drive due to chemical drive
% - Varousal: 1*N Array. The ventilatory drive due to arousal
% - LoopGainParameters: 1*8 Array of the derived loop gain parameters: 
%           - Element 1: Steady State Loop Gain 
%           - Element 2: Time Constant
%           - Element 3: Delay
%           - Element 4: Loop gain at the natural frequency
%           - Element 5: The natural cycling period
%           - Element 6: Loop gain at 1 cycle per minute
%           - Element 7: Loop gain at 2 cycles per minute
%           - Element 8: Additional ventilatory drive during breaths
%           occuring during an arrousal. 
% - BestParametersForEachDelay: 1*5 Array of the fitted ventilatory
% control model parameters:
%           - Element 1: Steady State Loop Gain 
%           - Element 2: Time Constant
%           - Element 3: Delay
%           - Element 4: Error in observed ventilattion to calculate Vchem_1
%           - Element 5: Additional ventilatory drive during breaths
%           occuring during an arrousal. 
% - BestSSres: 1*1 Array. Mean Square Error (SSres) of the model fit to observed
% ventilation data.

% - FitQuality: 1*4 Array. Summary Measures of the quality of the fit between
% the model and observed ventilation
%           - Element 1: Mean Square Error (SSres) of the model fit to observed
%           ventilation data.
%           - Element 2: R_squared for the correlation between the observed
%           ventilation and model predicted ventilation.
% - i_1: The index of the first breath (n=1) in the ventilation series for which Vchem and Vdrive are estimated
% - i_end: The index of the last breath in the ventilation series for which
% the model was applied
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%
%**************************************************************************
%**************************************************************************

Minute_Ventilation=Data(:,1); % The Minute Ventilation array is a breath-by-breath minute ventilation series which is normalised 
ObstructiveEvents=Data(:,2); % The ObstructiveEvents array defines whether each breath is obstructed. Penalize[n]=1 for an un-obstructed breath; Penalize[n]=0 for an obstructed breath.
TimeB=Data(:,4); % The time of each breath, from the first breath in the series. 
Ttot_ser=Data(:,6); % T_tot Series
Ttot_median=median(Ttot_ser); % The median t-tot based on the T_tot series [diff(TimeB)]. 


%**************************************************************************
% Set the constraints for the optimisation procedure:
%     [LG_0         tau         Err1[0]          VRA        Delay    ]
%**************************************************************************
start=[ 4           3          0                0           NaN   ]; % A reasonable estimate of a starting physiological value for each parameter
lower=[ 0.1         2           -Veupnea*3      0           1     ]; % Minimum conveivable physiological value for each parameter
upper=[ 30          120         Veupnea*3       Veupnea*3   5     ]; % Maximum conveivable physiological value for each parameter
%**************************************************************************
%**************************************************************************
    
%**************************************************************************
% The binary variable VraOn allows the ventilatory response to arousal
% parameter of the model to be selected - i.e., as per supplementary methods:
%    Vdrive=Vchem+Varousal                                 (Equation 1)
%    Vdrive[n]= alpha*Vchem[n?1]+beta*VE'[n]+gamma*Ar[n]   (Equation S5)
% if VraOn=0, gamma is set to zero, and contribution to ventilatory drive
% is assumed to be entirely due to chemical drive. 
%**************************************************************************
if VraOn==0
    start(4)=0;
    lower(4)=0;
    upper(4)=0;
end
%**************************************************************************
%**************************************************************************


%**************************************************************************
% To ensure model fitting is applied to identical sets of observed
% ventilation regardless of the delay time, we need to identify the first
% breath in the series for which there is adequate previous breaths in the
% series to construct a delayed ventilation series
%**************************************************************************
%delay_secs_Max=upper(5)*Ttot_median; % The maximum delay constraint in seconds
delay_secs_Max=5*Ttot_median; % The maximum delay constraint in seconds
i_min=4; 
while delay_secs_Max>TimeB(i_min)
    i_min=i_min+1;
end
i_1=10; % The index of the first breath for which the model will be fit.
%i_1=i_min; % The index of the first breath for which the model will be fit.
i_end=length(Minute_Ventilation);
%**************************************************************************
%**************************************************************************

%**************************************************************************
% Count the number of obstructive events within the analysis window. If
% there is not at least one obstructive event commencing within the window,
% pass all output parameters as NaN's
%**************************************************************************
N_events=0;
for k=i_1:length(ObstructiveEvents)
    if ObstructiveEvents(k)==0&&ObstructiveEvents(k-1)==1
        N_events=N_events+1;
    end
end
if  N_events<1&&0
    LoopGainParameters(1:8)=NaN;
    FitQuality(1:4)=NaN;
    Error=NaN*Minute_Ventilation;
    Vchem=NaN*Minute_Ventilation;
    Varousal=NaN*Minute_Ventilation;
    BestSSres=NaN;
    BestParameters=NaN;
    i_1=NaN;
    i_end=NaN;
    return
end
%**************************************************************************
%**************************************************************************

delay_breaths=lower(5):upper(5); % Use the constraints to create a vector of candidate delay values
start(5)=[];
lower(5)=[];
upper(5)=[];
for dd=1:length(delay_breaths)

    %**********************************************************************
    % Generate the delayed VE:
    %**********************************************************************
    delay_secs=delay_breaths(dd)*Ttot_median;
    delayed_VE{dd}=zeros(length(Minute_Ventilation),1);
    for k=i_1:length(Minute_Ventilation)
        d_i=Ttot_ser(k-1);
        k_back=1;
        while delay_secs>d_i
            k_back=k_back+1;
            d_i=sum(Ttot_ser(k-k_back:k-1));
        end
        delayed_VE{dd}(k)=Minute_Ventilation(k-k_back)+(1-(delay_secs-sum(Ttot_ser(k-k_back+1:k-1)))/Ttot_ser(k-k_back))*(Minute_Ventilation(k-k_back+1)-Minute_Ventilation(k-k_back)); 
        delayed_VE{dd}(1:i_1-1)=-1;
    end
    %**********************************************************************
    %**********************************************************************

    %**************************************************************************
    % In order to identify the optimal parameter set (of LG_0, tau, Vdr[1],
    % VRA) the inbuilt MATLAB function fmincon.m is utilised with an interior
    % point algorithm.
    % 
    % To maximize the robustness of the fit (avoiding local minima), a set of
    % five starting conditions (initial parameter values) were applied. The
    % optimisation algorithm is applied for up to 400 iterations. 
    %
    % The starting set that yielded the the lowest final SSres is identified. The
    % final end parameter set for this optimal starting point is now used for
    % up to 1000 further iterations to identify the final optimal paramter set:
    %**************************************************************************

    % Trial 5 initial conditions:
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    OPTIONS = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',10,'MaxFunEvals',10,'Algorithm','interior-point');
    
    XX=100;
    Parameter_Mat=zeros(length(upper),XX);
    SSres_Mat=zeros(1,XX);
    for i=1:XX
    Parameters = rand(1,length(upper)).*(upper-lower)+lower;
    [Parameter_Mat(:,i),SSres_Mat(i),~,~] = fmincon(@(Parameters) TheModel(Parameters,Data,i_1,delayed_VE{dd},polyfitorder),Parameters,[],[],[],[],lower,upper,@(Parameters) model_nonlinear_constraints(),OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % Identify the best of the 5 initial starting points:
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    [~,ind_k]=min(SSres_Mat);
    Parameters=Parameter_Mat(:,ind_k);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    % Use the end-point of the best starting condition as the start point for further iterations:
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    OPTIONS = optimset('Display','off','TolX',1e-6,'TolFun',1e-6,'MaxIter',50,'MaxFunEvals',50,'Algorithm','interior-point');
    [BestParametersForEachDelay(:,dd),Optimal_SSres(dd),~,~] = fmincon(@(Parameters) TheModel(Parameters,Data,i_1,delayed_VE{dd},polyfitorder),Parameters,[],[],[],[],lower,upper,@(Parameters) model_nonlinear_constraints(),OPTIONS);
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %**************************************************************************
    %**************************************************************************

end
[~, ind]=min(Optimal_SSres);
BestParameters=BestParametersForEachDelay(:,ind);
[BestSSres,Error,Rsq,Vchem,Varousal] = TheModel(BestParameters,Data,i_1,delayed_VE{ind},polyfitorder);
FitQuality=[BestSSres Rsq ind_k]; % This array compiles the parameters describing the quality of the model fit to observed ventilation data
BestParameters(5)=delay_breaths(ind); % The optimal delay (breaths)

%**************************************************************************
% Following the optimisation procedure the optimal values of LG_0, tau, Vdr[1],
% gamma and delay are now known. Equation S4 can now be used to calculate
% the dynamic loop gain parameters LG1, LG2, LGn and natural cycling
% frequency:
% LG(f)=-LG_0*e^(delay*i*2*pi*f)/(1+i*2*pi*f*tau)
%**************************************************************************
LG0=BestParameters(1); % The steady state loop gain
tau=BestParameters(2); % The estimated time-constant, Tau
gamma=BestParameters(4)/Veupnea; % The ventilatory response to arousal (VRA)
delay=BestParameters(5)*Ttot_median; % Convert the optimal delay from number of breaths back to a duration in seconds my multiplying by the patients median T_tot

Ff=[0.002:0.00001:1.09999 1.1.^(1:100)];
phase=pi-atan(2*pi*tau*Ff)-2*pi*delay*Ff;
fn=interp1(phase,Ff,0,'spline'); % The natural cycling frequency
Tn=1/fn; % The natural cycling period
LGn=abs(LG0*exp(-1i*2*pi*fn*delay)/(1+1i*2*pi*fn*tau)/(1)); % The loop gain at the natural frequency
LG1=abs(LG0*exp(-1i*2*pi*(1/60)*delay)/(1+1i*2*pi*(1/60)*tau)/(1)); % The loop gain at 1 cycles per minute
LG2=abs(LG0*exp(-1i*2*pi*(1/30)*delay)/(1+1i*2*pi*(1/30)*tau)/(1)); % The loop gain at 2 cycles per minute
%**************************************************************************
%**************************************************************************

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
LoopGainParameters=[LG0 tau delay LGn Tn LG1 LG2 gamma];
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%Vchem=Vchem+Veupnea;

end



%%
function [c,ceq] = model_nonlinear_constraints()
c=[];
ceq = [];
end