function Pes_clean = PesRemoveEKG(Pes,EKG,dt)
%load('J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\SpikeFiles\1341.mat', 'Pes', 'EKG');
%Pes.values(end:end+1)=Pes.values(end);
%Pesclean = PesRemoveEKG(Pes.values,EKG.values,Pes.interval)
%addpath('E:\Work\MatlabFunctionDatabase')
polyorder=3;
Time = (0:dt:(dt*(length(Pes)-1)))';
Niterations=3;
Nbeatspersegment = 200;

ECG_peak_i = EKGpeakdetection(EKG,Time,dt);
        medianRR = median(diff(ECG_peak_i))*dt;
        
        medianRRi = round(medianRR/dt);
    lefttemplatetime = 0.33*medianRR; righttemplatetime = 1.67*medianRR;
    leftdi1=round(1/dt*lefttemplatetime); rightdi1=round(1/dt*righttemplatetime);
    F1 = (leftdi1+rightdi1)/medianRRi;
    contaminationmagnitudethreshold = 0.1;
    
    
% clear temp contaminationmagnitude contaminationmagnitudepost
% processPes=1;
% loadprocessedPesifexists = 0; saveprocessedPes = 0; foundamatch=0;
% if loadprocessedPesifexists
%         foundamatch=0;
%         temp='Pes_clean';
%         for j=1:length(w)
%             foundamatch=strcmp(w(j).name,temp);
%             if foundamatch
%                 eval([temp '=filehandle.' temp ';']);
%                 break
%             end
%         end
%     if foundamatch
%         [~,template] = crosscontaminated(Pes,ECG_peak_i,leftdi,rightdi,1);
%         contaminationmagnitude=max(abs(template));
%         [~,template] = crosscontaminated(Pes_clean,ECG_peak_i,leftdi,rightdi,1);
%         contaminationmagnitudepost=max(abs(template));
%     end
% end
%%
% if (~foundamatch||loadprocessedPesifexists==0)&&processPes
try clf(111), catch me, end

Pes_clean=Pes;
plotalready=0;
    dividenightinto = floor(length(ECG_peak_i)/Nbeatspersegment);
    
for n=1:dividenightinto
    %%
    beatleft = Nbeatspersegment*(n-1)+1; beatright = beatleft + Nbeatspersegment + 1; %last 2 beats overlap
    
    if beatright>length(ECG_peak_i), beatright=length(ECG_peak_i); end
        ileft = ECG_peak_i(beatleft);
        iright = ECG_peak_i(beatright)+1;
    
    if n==1, ileft=1; end
    if n==dividenightinto, iright=length(EKG); end
    
    ECG_peak_i_segment = ECG_peak_i(beatleft:beatright) - ileft + 1;
    
    Pes_clean_segment=Pes_clean(ileft:iright);
    try
    leftdi=leftdi1; rightdi=rightdi1;
    
    for repeat=1:Niterations %repeats template subtraction with new template
        % EKG artifact removal
        [~,template] = crosscontaminated(Pes_clean_segment,ECG_peak_i_segment,leftdi,rightdi,1,polyorder,0);
        
        contaminationmagnitude=max(abs(template));
        
        %disp(['Decontaminating ' num2str(sum(contaminationmagnitude>contaminationmagnitudethreshold)) ' Pes signal']);
        if contaminationmagnitude>contaminationmagnitudethreshold
            Pes_clean_segment=crosscontaminated(Pes_clean_segment,ECG_peak_i_segment,leftdi,rightdi,0,polyorder,0);
        end
        
        % after
        %[~,template] = crosscontaminated(Pes_clean_segment,ECG_peak_i_segment,leftdi,rightdi,1);
        %contaminationmagnitudepost(i)=max(abs(template));
        if repeat<Niterations
            leftdi=round(leftdi*1.1); rightdi=round(rightdi*1.1);
            leftdi=leftdi+leftdi; rightdi=rightdi-leftdi;
        end
    end
    catch me
        disp(me.message)
    end
    Pes_clean(ileft:iright)=Pes_clean_segment;
    if 1
        %if ~plotalready
        figure(111)
        if ~plotalready
            axx(1)=subplot(2,1,1); plot(Time,EKG);
                axx(1).FontSmoothing = 'off';
            axx(2)=subplot(2,1,2); plot(Time,Pes,'k'); hold('on');
                axx(2).FontSmoothing = 'off';
            plotalready=1;
            linkaxes(axx,'x');
            xlim([Time(ileft)-5 Time(iright)+5]);
        end
        axx(2)=subplot(2,1,2); plot(Time(ileft:iright),Pes_clean(ileft:iright),'r');
        xlim([Time(ileft)-5 Time(iright)+5]);
    n
    end
end
%     if saveprocessedPes
%         currentdir = cd;
%         savestring = ['Pes_clean'];
%         temp = '-append';
%         cd(directory);
%         eval(['save ' filename ' ' savestring ' ' temp]);
%         cd(currentdir);
%     end
% else %no decontamination needed, just copy Pes into Pes_clean
%     Pes_clean=Pes;
% end

%%
if 1
    try clf(111), catch me, end
    figure(111)
    axx(1)=subplot(2,1,1); plot(Time,EKG);
    axx(2)=subplot(2,1,2); plot(Time,Pes,'r',Time,Pes_clean,'b');
    linkaxes(axx,'x');
end

