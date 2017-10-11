function [NL Jdate dur]=getMedianHourlyPGNL(fname, fband_idx)

% Just as getPGNL, gets the median hourly pamguard noise measures in band
% of intrest

curdir=cd();

cd('C:\Users\charlotte\Desktop\From Luke\BinaryFiles');
fnames=dir(fullfile(fname,'Noise_Monitor*pgdf'));

NL=[];
Jdate=[];

for ii=1:length(fnames)
    load_name=[fname '\' fnames(ii).name];
    [noises] = loadBandNoiseFile(load_name);

    if ~isempty(noises)
        for jj=1:length(noises)
            NL_1(jj)=median(noises(jj).noise(:,fband_idx));
            Jdate_1(jj)=noises(jj).date;
            NL=[NL NL_1]; 
            Jdate=[Jdate Jdate_1];
        end
          % Get the duration in seconds of the ambient noise measure
    end
    
    if length(Jdate_1)<2
        dur=0;
    else
    dur=(Jdate_1(2)-Jdate_1(1))*24*60*60;
    end
end
 

% For each Hour get the median
% Determine how many hours
Jdate_hrs=floor(Jdate*24);
N_hrs=ceil((max(Jdate_hrs)-min(Jdate_hrs)));
Hr_idx=unique(Jdate_hrs);
MedNL=zeros(N_hrs,1);
for jj=1:length(Hr_idx)
    
    MedNL(jj)=median(NL(Jdate_hrs==Hr_idx(jj)));   
    
end


NL=MedNL;
NL(NL==0)=[];
Jdate=unique(Jdate_hrs)/24;
cd(curdir);

end