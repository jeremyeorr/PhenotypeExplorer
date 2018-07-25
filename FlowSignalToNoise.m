function noisewav = FlowSignalToNoise(timewav,respwav,showfigures)

Fs = 1./(timewav(2)-timewav(1));

runStoNdata=1;
scalenoise=4; %default 1
if runStoNdata
    %signal to noise per window
    secslide2=10;
    windur2 = 120; %note no savings in speed by keeping nfft = 2^X (actually lost speed).
    % nfft is usually set as next larger power of 2, as:
    % nfft = 2^nextpow2(N); % next larger power of 2
    nwin=floor((timewav(end)-timewav(1)-windur2)/secslide2 + 1);
    StoN_x=zeros(1,nwin);
    P1_x=zeros(1,nwin);
    Time_x=zeros(1,nwin);
    for i=1:nwin%winNum=winnumrange
        lefti=1+(i-1)*secslide2*Fs;
        righti=lefti+windur2*Fs-1;
        Flow=respwav(lefti:righti);
        Time=timewav(lefti:righti);
        Flowfft = fft(Flow-mean(Flow))*2/length(Flow);
        Flowpsd = Flowfft.*conj(Flowfft)*2/length(Flow);
        if i==1
            df = 1/(Time(end)-Time(1));
            F = 0:df:(length(Time)-1)*df; Fmax = 20;
            Fmaxi = round(Fmax/df+1);
            ranges = [0.1 1 1 10];
            rangesi = round(ranges/df)-1;
        end
        F((Fmaxi+1):end)=[]; Flowpsd((Fmaxi+1):end)=[]; Flowpsd(1)=0;
        P1 = sum(Flowpsd(rangesi(1):rangesi(2))); %note not divided by bandwidth
        P2 = sum(Flowpsd(rangesi(3):rangesi(4))); %note not divided by bandwidth
        StoN = P1/P2;
        StoN_x(i)=StoN;
        P1_x(i) = P1;
        Time_x(i) = Time(1); 
    end
    P1scale = 10^mean(log10(P1_x(StoN_x>10)));
    criteria = P1_x/P1scale.*StoN_x;
    noisewav=0*timewav;
    for i=1:nwin
        if criteria(i)<(1/scalenoise)
            lefti=1+(i-1)*secslide2*Fs;
            righti=lefti+windur2*Fs-1;
            noisewav(lefti:righti)=1;
        end
    end
    for i=1:nwin
        if criteria(i)<(0.5/scalenoise)
            lefti=1+(i-1)*secslide2*Fs;
            righti=lefti+windur2*Fs-1;
            noisewav(lefti:righti)=2;
        end
    end
    for i=1:nwin
        if criteria(i)<(0.25/scalenoise)||isnan(criteria(i))
            lefti=1+(i-1)*secslide2*Fs;
            righti=lefti+windur2*Fs-1;
            noisewav(lefti:righti)=3;
        end
    end  
    StoNData.Fnoiseovernight = [sum(noisewav>=1) sum(noisewav>=2) sum(noisewav>=3)]/length(noisewav);
    StoNData.StoNovernight = [sum(StoN_x<20) sum(StoN_x<10) sum(StoN_x<5)]/length(noisewav);
    if showfigures %StoN figure
        figure(101); 
        set(gcf,'color',[1 1 1]);
        ax30(1)=subplot(4,1,3); stairs(Time_x,P1_x/P1scale);
        set(gca,'yscale','log');
        hold('on');
        stairs(Time_x,StoN_x,'k'); box('off');
        stairs(Time_x,P1_x/P1scale.*StoN_x,'r'); box('off');
        stairs(timewav,10.^noisewav,'g'); box('off'); hold('off');
        legend('Power','StoN','PxStoN','lowQ');
        ax30(2)=subplot(4,1,4); plot(timewav,respwav); ylabel('Flow'); box('off');
        %ax30(3)=subplot(3,1,3); plot([1:nwin],NoiseData); ylabel('Noise'); box('off');
        linkaxes(ax30(1:2),'x');
    end
end