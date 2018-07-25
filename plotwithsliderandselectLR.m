function plotwithsliderandselectLR(tminmax)
global P h1 valuelast xlimlast sliderstep xvalues yvalues range ax2
sliderstep = [0.025 1];

%range = 2; % visible x range
P = [tminmax(1)+range/2 tminmax(2)-range/2]; % min and max mean(xlim Positions)
axis tight; xlim([tminmax(1) tminmax(1)+range])
%plot(x,y); 

set(gca,'ylimmode','auto');
set(gca,'ylimmode','manual');
h1=uicontrol('style','slider','units','normalized','position',[0 0 1 0.05],'KeyPressFcn',@key_pressed_fcn,...
    'callback',@slider1);
valuelast=get(h1,'Value');
xlimlast=get(gca,'xlim');
set(h1,'sliderstep',sliderstep*range/(P(2)-P(1)))
uicontrol(...
    'Units','normalized',...
    'Position',[0.01,0.99-0.025,0.15,0.025],...
    'String','Select',...
    'KeyPressFcn', @key_pressed_fcn,...
    'Callback',@selectdata);
uicontrol(...
    'Units','normalized',...
    'Position',[0.17,0.99-0.025,0.15,0.025],...
    'String','Optimize',...
    'KeyPressFcn', @key_pressed_fcn,...
    'Callback',@optimize);
uicontrol(...
    'Units','normalized',...
    'Position',[0.17+0.16,0.99-0.025,0.15,0.025],...
    'String','ZoomInX',...
    'KeyPressFcn', @key_pressed_fcn,...
    'Callback',@zoominX);
uicontrol(...
    'Units','normalized',...
    'Position',[0.17+2*0.16,0.99-0.025,0.15,0.025],...
    'String','ZoomOutX',...
    'KeyPressFcn', @key_pressed_fcn,...
    'Callback',@zoomoutX);
uicontrol(...
    'Units','normalized',...
    'Position',[0.17+3*0.16,0.99-0.025,0.15,0.025],...
    'String','ClrLastSelect',...
    'KeyPressFcn', @key_pressed_fcn,...
    'Callback',@ClrLastSelect);
set(gcf,'KeyPressFcn', @key_pressed_fcn);
end

function slider1(hObj,event)
global P valuelast xlimlast sliderstep range
% called when slider was used
global ax2
for i=1:length(ax2)
    set(ax2(i),'ylimmode','manual');
end
    valtemp=get(hObj,'value');
    rangeold = range;
    attemptedwindowstepnorm=(valtemp-valuelast)*(P(2)-P(1))/rangeold;
    
    xlimnew=get(gca,'xlim');
    range=xlimnew(2)-xlimnew(1);
    set(hObj,'sliderstep',sliderstep*range/(P(2)-P(1)));
    if sum(xlimnew==xlimlast)==0 %xrange moved not via slider
        mxl = xlimnew(1)+range/2;
        val = (mxl-P(1))/(P(2)-P(1)); %first move slider to appropriate location
        set(hObj,'value',val);
        switch round(1000*attemptedwindowstepnorm)/1000 %minor rounding issues
            case sliderstep(2) %large step forward
            val=val+range*sliderstep(2)/(P(2)-P(1));
            case sliderstep(1) %small step forward
            val=val+range*sliderstep(1)/(P(2)-P(1));
            case -sliderstep(2) %large step backward
            val=val-range*sliderstep(2)/(P(2)-P(1));    
            case -sliderstep(1) %small step backward
            val=val-range*sliderstep(1)/(P(2)-P(1));    
        end
        set(hObj,'value',val);
    else %move x limits simply according to slider
        val = get(hObj,'value');
    end
    mxl = val*(P(2)-P(1))+P(1); % mean xlim, denormalize
    xlim([mxl-range/2 mxl+range/2]);
    xlimlast=get(gca,'xlim');
    valuelast=val;
end

function selectdata(varargin)
global xvalues yvalues
[tempx,tempy] = ginput(2);
xvalues = [xvalues;tempx'];
yvalues = [yvalues;tempy'];
outputlive = [xvalues yvalues]
hold(gca,'on');
plot([xvalues(end,1) xvalues(end,2) xvalues(end,2) xvalues(end,1) xvalues(end,1)],[yvalues(end,1) yvalues(end,1) yvalues(end,2) yvalues(end,2) yvalues(end,1)],'color',[1 0 0]);
hold(gca,'off');
end %inputCallback

function ClrLastSelect(varargin)
global xvalues yvalues
hold(gca,'on');
if xvalues(end,1)~=xvalues(end,2)||yvalues(end,1)~=yvalues(end,2)
    plot([xvalues(end,1) xvalues(end,2) xvalues(end,2) xvalues(end,1) xvalues(end,1)],[yvalues(end,1) yvalues(end,1) yvalues(end,2) yvalues(end,2) yvalues(end,1)],'color',[1 0.8 0.8]);
else
    plot([xvalues(end,1)],[yvalues(end,1)],'.','color',[1 0.8 0.8]);
end    
hold(gca,'off');
xvalues(end,:)=[];
yvalues(end,:)=[];
outputlive = [xvalues yvalues]
end %inputCallback


function optimize(varargin)
global ax2
for i=1:length(ax2)
    set(ax2(i),'ylimmode','auto');
end
% figure(gcf);
% set(gca,'ylimmode','manual');
end %inputCallback

function zoominX(varargin)
global range h1 sliderstep P
xlimnew=get(gca,'xlim');
xlimdelta=xlimnew(2)-xlimnew(1);
xlimnew(2)=xlimnew(1)+xlimdelta/2;
set(gca,'xlim',xlimnew);
  
    range=xlimnew(2)-xlimnew(1);
    set(h1,'sliderstep',sliderstep*range/(P(2)-P(1)));
    mxl = xlimnew(1)+range/2;
    val = (mxl-P(1))/(P(2)-P(1)); %first move slider to appropriate location
    set(h1,'value',val);
    valuelast=val;
end %inputCallback

function zoomoutX(varargin)
global range h1 sliderstep P valuelast
xlimnew=get(gca,'xlim');
xlimdelta=xlimnew(2)-xlimnew(1);
xlimnew(2)=xlimnew(1)+xlimdelta*2;
set(gca,'xlim',xlimnew);

    range=xlimnew(2)-xlimnew(1);
    set(h1,'sliderstep',sliderstep*range/(P(2)-P(1)));
    mxl = xlimnew(1)+range/2;
    val = (mxl-P(1))/(P(2)-P(1)); %first move slider to appropriate location
    set(h1,'value',val);
    valuelast=val;
end %inputCallback

function key_pressed_fcn(varargin)
global h1 range sliderstep P xvalues yvalues
switch varargin{2}.Key
    case 'c'
        zoomoutX();
    case 'x'
        zoominX();
    case 's'
        selectdata();
    case 'z'
        optimize();
    case 'rightarrow'
        %'rightarrow'
        val=get(h1,'value');
        val=val+range*sliderstep(1)/(P(2)-P(1));
        set(h1,'value',val);
        slider1(h1);
    case 'leftarrow'
        %'rightarrow'
        val=get(h1,'value');
        val=val-range*sliderstep(1)/(P(2)-P(1));
        set(h1,'value',val);
        slider1(h1);    
    case 'pageup'
        %'rightarrow'
        val=get(h1,'value');
        val=val-range*sliderstep(2)/(P(2)-P(1));
        set(h1,'value',val);
        slider1(h1);   
    case 'pagedown'
        %'rightarrow'
        val=get(h1,'value');
        val=val+range*sliderstep(2)/(P(2)-P(1));
        set(h1,'value',val);
        slider1(h1);   
    case 'q' %pageup
        %'rightarrow'
        val=get(h1,'value');
        val=val-range*sliderstep(2)/(P(2)-P(1));
        set(h1,'value',val);
        slider1(h1);   
    case 'w' %pagedown
        %'rightarrow'
        val=get(h1,'value');
        val=val+range*sliderstep(2)/(P(2)-P(1));
        set(h1,'value',val);
        slider1(h1);
    case 'n' %jump to next arousal
        It=evalin('base','It'); %isempty(eval(signallist{i}))
        t_wincenter = get(gca,'xlim') + range/2;
        t_wincenter(2)=[];
        I=find(It>(t_wincenter+1),1); %add 1 sec for certain right shift
        if ~isempty(I)
        try 
            val=(It(I)-P(1))/(P(2)-P(1)); 
            set(h1,'value',val);
            slider1(h1);
        catch me, 
        end
        end
    case 'b' %jump to previous arousal
        It=evalin('base','It'); %isempty(eval(signallist{i}))
        t_wincenter = get(gca,'xlim') + range/2;
        t_wincenter(2)=[];
        I=find(It<=(t_wincenter-1),1,'last');
        if ~isempty(I)
        try 
            val=(It(I)-P(1))/(P(2)-P(1)); 
            set(h1,'value',val);
            slider1(h1);
        catch me, 
        end
        end
    case 'e'
        [tempx,tempy] = ginput(1);
        xvalues = [xvalues;[tempx tempx]];
        yvalues = [yvalues;[tempy tempy]];
        outputlive = [xvalues yvalues]
        hold(gca,'on');
        plot([xvalues(end,1)],[yvalues(end,1)],'.','color',[1 0 0]);
        hold(gca,'off');
end
end