function [NL MDATE dur]=getMeanMinutePGNL(fname, fband_idx)




curdir=cd();

%cd('C:\Users\charlotte\Desktop\From Luke\BinaryFiles');
fnames=dir(fullfile(fname,'Noise_*pgdf'))

NL=[];
MDATE=[];

for ii=1:length(fnames)
    load_name=[fname '\' fnames(ii).name];
    [noises] = loadBandNoiseFile(load_name);
    

    if length(noises)>0 
        NL_temp=[];
        Jdate_temp=[];
        for jj=1:length(noises)
            NL_temp(jj)=log10(mean(10.^noises(1).noise(:,fband_idx)));
            Jdate_temp(jj)=noises(jj).date;
        end
        NL=[NL NL_temp];
        MDATE=[MDATE Jdate_temp];
        clear NL_temp Jdate_temp
    end
    
    
    
end
 


NL(NL==0)=[];
MDATE(MDATE==0)=[];
cd(curdir);

end