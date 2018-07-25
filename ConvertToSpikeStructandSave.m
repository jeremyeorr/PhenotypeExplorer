function ConvertToSpikeStructandSave(varlist,data,dt,directory,filename)

%         for i=1:length(EEGsignals)
%             varlist{i} = [EEGsignals{i} '_clean'];
%             data{i}=eval(varlist{i});
%         end
        
        
        currentdir = cd;
        savestring=[];
        for n=1:length(varlist)
            eval([varlist{n} '.values=data{n};']);
            eval([varlist{n} '.interval=dt;']);
            if n>1
                temp = ' ';
            else
                temp = '';
            end
            savestring = [savestring temp varlist{n}];
        end
        temp = '-append';
        cd(directory);
%         eval(['save([directory filename],' savestr ',temp' ');']);
        eval(['save ' filename ' ' savestring ' ' temp]);
        cd(currentdir);

        