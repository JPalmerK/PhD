function [NL Jdate dur]=getPGNL(fname, fband_idx)




curdir=cd();

%cd('C:\Users\charlotte\Desktop\From Luke\BinaryFiles');
fnames=dir(fullfile(fname,'Noise_*pgdf'))

NL=[];
Jdate=[];

for ii=1:length(fnames)
    load_name=[fname '\' fnames(ii).name];
    [noises] = loadBandNoiseFile(load_name);
    

    if length(noises)>0 
        for jj=1:length(noises)
            if length(noises(jj).noise)> fband_idx
            %Median of the 10 minute band
            NL_1(jj)=median(noises(jj).noise(:,fband_idx)); 
            Jdate_1(jj)=noises(jj).date;
            NL=[NL NL_1]; 
            Jdate=[Jdate Jdate_1];
            end
        end
          % Get the duration in seconds of the ambient noise measure
    %dur=(Jdate_1(2)-Jdate_1(1))*24*60*60;
    end
    
    
end
 


NL(NL==0)=[];
Jdate(Jdate==0)=[];
cd(curdir);

end