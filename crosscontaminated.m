function [decontaminated,mediantemplate] = crosscontaminated(contaminated,index,leftdi,rightdi,templateonly,polyorder,disppercentcomplete)
plotany=0;
decontaminated = contaminated;
%%
mag=1; %degree of noise to subtract, 1=full, less is more conservative
%polyorder=5;
warning('off');

templatedata = NaN*zeros(length(index),leftdi+1+rightdi);
if 1
    w = hann(leftdi+1+rightdi);
else
    w = tukeywin(leftdi+1+rightdi,0.33);
end
    for i=1:length(index)
        try
            lefti = index(i)-leftdi;
            righti = index(i)+rightdi;    
            templatedata(i,:)=contaminated(lefti:righti)-mean(contaminated(lefti:righti));
        catch me 
            disp me
        end
    end
    mediantemplate = nanmedian(templatedata)-nanmedian(nanmedian(templatedata));

    if disppercentcomplete
    try clf(99); catch me, end
    figure(99); 
    plot(mediantemplate,'r'); hold('on');
    plot(mediantemplate.*w','k'); title('template'); hold('off');
     if templateonly
         decontaminated=[];
         return
     end
    end
%%  

mediantemplate = mediantemplate.*w';
clear M
    count=1;
    t_start = clock;
    ploteveryN=round(length(index)/2); 
    if plotany
        figure(54);
    end
    ratios = [0 0.2 0.5 1 1.5 2 4]; b1 = -1:1;
    %errorterm = NaN*zeros(length(index),length(ratios));
    tic
    for i=1:length(index)
        walls=0;
        for b=1:length(b1)
            clear errorterm
            for j=1:length(ratios)
                lefti = index(i)-leftdi+b1(b);
                righti = index(i)+rightdi+b1(b);  
                if lefti<1||righti>length(contaminated)
                    walls=1;
                    break
                end
                temp = contaminated(lefti:righti);

                errortemp = temp-ratios(j)*mediantemplate';
                
                if polyorder==1&&1 %linear detrend, fastest
                    temp2=detrend(errortemp,'linear'); %detrend residual "corrected" trace
                else %polynomial detrend
                    if 0 %slowest
                    polyp = polyfit((1:length(errortemp))',errortemp,polyorder);
                    temp2 = errortemp - polyval(polyp,(1:length(errortemp))');
                    else
                        y=errortemp;
                        x=(1:length(errortemp))';
                        if ~exist('M','var') %doesn't change each step so no need to rewrite
                            M = repmat(x,1,polyorder+1); M = bsxfun(@power,M,polyorder:-1:0); %credit to Matt Tearle on 23 Feb 2011, see http://www.mathworks.com/matlabcentral/answers/1836-multiple-use-of-polyfit-could-i-get-it-faster
                        end
                        polyp = M\y;
                        polytrend = polyval(polyp,(1:length(y))');
                        temp2 = y - polytrend;
                    end
                
                end
                
                if 1
                    errorterm(j)=std(temp2);
                else
                    errorterm(j)=std(diff(temp2));
                end
            end %k
            if walls
                break
            end
            [error(b),mini]=min(errorterm);
            if mini>1&&mini<length(ratios)
                [ratio(b),error(b)]=PeakFitQuadratic(ratios(mini-1:mini+1),errorterm(mini-1:mini+1)); %[error,ratio]
            else
                ratio(b)=ratios(mini); %ratio
            end
            if mod(i,ploteveryN)==1&&plotany

                subplot(2,1,2);semilogx(ratios,errorterm,'.',ratio(b),error(b),'ko');
                hold('on');
            end
        end
        if walls
            continue
        end
        [minerror,bestbi] = min(error); bestb = b1(bestbi); bestratio = ratio(bestbi);
            lefti = index(i)-leftdi+bestb;
            righti = index(i)+rightdi+bestb;   
            temp = contaminated(lefti:righti);
            temp2=(temp-mag*bestratio*(mediantemplate)');
            decontaminated(lefti:righti) = temp2;
            
            
            if disppercentcomplete&&ceil(i/length(index)*100)>count %display percentage complete
                disp([num2str(ceil(i/length(index)*100)) '% complete']);
                count=count+1;
                disp([num2str(etime(clock,t_start)/count*100 - etime(clock,t_start)) 's remaining']);
            end
            
            if plotany
            if mod(i,ploteveryN)==1
                
            if 1 %plot trend
            y=temp2;
            x=(1:length(y))';
            if ~exist('M') %doesn't change so no need to rewrite
                M = repmat(x,1,polyorder+1); M = bsxfun(@power,M,polyorder:-1:0); %credit to Matt Tearle on 23 Feb 2011, see http://www.mathworks.com/matlabcentral/answers/1836-multiple-use-of-polyfit-could-i-get-it-faster
            end
            polyp = M\y;
            polytrend = [NaN*(1:10)';polyval(polyp,(1:length(y))');NaN*(1:10)'];
            end
            
            figure(54);
            subplot(2,1,1); 
                plot(contaminated(lefti-10:righti+10),'r'); hold('on');
                plot(decontaminated(lefti-10:righti+10)); 
                plot(polytrend,'g'); hold('off');
                xlabel('Sample Number'); ylabel('Signal'); legend('Contaminated','Decontaminated','Trend');
            subplot(2,1,2);semilogx(bestratio,minerror,'ro'); 
            xlabel('Fraction of template subtracted'); ylabel('Residual S.D.');
                hold('off');   
                pause(0.01)
            end
            end
    end