function PlotAndSelectData(signallist,Time,clearvalues)
global ax2 xvalues yvalues range
try clf(2), catch me, end

dsf=5;

if clearvalues
    xvalues=[];
    yvalues=[];
end

figure(2);
set(gcf,'color',[1 1 1]);
signalnothere=[];
Nsignals=length(signallist);
n=1;
for i=1:Nsignals
    try
        signal{n}=evalin('base',signallist{i}); %isempty(eval(signallist{i}))
        n=n+1;
    catch me
        disp(me.message);
        signallist{i}=[];
    end
end

X=length(signal);

ax2=[];
for i=1:length(signal)
    ax2(i)=subplot(X,1,i); plot(downsample(Time,dsf),downsample(signal{i},dsf)); ylabel(signallist{i});
    postemp = get(gca,'Position');
    posnew = [postemp(1)/2 postemp(2)-0.5*(postemp(4)*0.3) postemp(3)+2*(postemp(1)/2) postemp(4)+postemp(4)*0.3];
    set(gca,'Position',posnew);
end
linkaxes(ax2,'x');

for i=1:length(ax2)-1
    set(ax2(i),'Xtick',[],'Xcolor',[1 1 1]);
end
for i=1:length(ax2)
    set(ax2(i),'tickdir','out','box','off','fontname','arial narrow');
    TickLength=get(ax2(i),'TickLength');
    set(ax2(i),'TickLength',TickLength/3);
end

plotwithsliderandselectLR([Time(1) Time(end)]);

