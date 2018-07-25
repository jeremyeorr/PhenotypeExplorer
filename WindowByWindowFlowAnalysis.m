function Data = WindowByWindowFlowAnalysis(Time,Flow,noisewav,CPAPoff)
    plotfigs=0;
    secslide=120;
    Fs = 1/(Time(2)-Time(1));
    WindowDuration=300;
    skipCPAPon = 1;
    numwind=floor((length(Flow)-WindowDuration*Fs)/(secslide*Fs)+1); % calculates the number window (of given size, and of given slide) required to encapsulate range
    Data = [];
    pcompletestep=1;
    pcompletenext=pcompletestep; 
    for winNum=1:numwind
        if 1
            if floor(winNum/numwind*100)>=pcompletenext
                disp([num2str(floor(winNum/numwind*100)) '% complete']);
                pcompletenext = pcompletenext + pcompletestep;
            end
        end
        li=(winNum-1)*secslide*Fs+1;
        ri=(winNum-1)*Fs*secslide+WindowDuration*Fs;
        I=(li:ri);
        Fnoise2=sum(noisewav(I)>=2)/length(noisewav(I));
        if Fnoise2>0.1
            continue
        end
        if skipCPAPon
            FCPAPoff=sum(CPAPoff(I))/length(CPAPoff(I));
            if FCPAPoff<0.99
                continue
            end
        end
    [~,~,BB_i_start,BB_i_mid,BB_i_end,BB_t,VI,VE,Ttot,leak,IEratio,VT] = Vflowanalysis3(Time(I),Flow(I),0,1-plotfigs);
    Data_ = [winNum+0*BB_i_start,BB_i_start,BB_i_mid,BB_i_end,BB_t',VI',VE',Ttot',VT'];
    Data = [Data;Data_];
    end
    
    %Data_backup = Data;