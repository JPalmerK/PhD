
close all; clear all; clc
cd('W:\KJP PHD\3-Detection Function\Propagation Model\Code')
Parent_loc='J:\GliderData for CH3 NL comp';

SMFBands=readtable('W:\KJP PHD\SM2M Processing\PAMGuardCenterFrequencyNLBands.csv');
Gliderfband=readtable('J:\GliderData for CH3 NL comp\PAMGuardCenterFrequencyNLBands.csv');

% Get index of Glider fband nearest the highest SM band
GliderBandIdx=[];
GliderBandIdx(2)=nearest2(Gliderfband.CenterFreq,SMFBands.CenterFreq(end));
GliderBandIdx(3)=GliderBandIdx(2)+1;
GliderBandIdx(1)=GliderBandIdx(2)-1;

bw=Gliderfband.CenterFreq(GliderBandIdx)*1.222-Gliderfband.CenterFreq(GliderBandIdx)*.891;


Parent_loc='J:\GliderData for CH3 NL comp\';
filedays=dir(fullfile([Parent_loc], '20*'));

NL_out1=[];
NL_out2=[];
NL_out3=[];

Jdate_out1=[];
Jdate_out2=[];
Jdate_out3=[];

    for jj=1:length(filedays)-1
        
        fdate=[Parent_loc filedays(jj).name];
        [NL_dn Jdate1]=getPGNL(fdate,GliderBandIdx(1));
        [NL_sm Jdate2]=getPGNL(fdate,GliderBandIdx(2));
        [NL_up Jdate3]=getPGNL(fdate,GliderBandIdx(3));
        
        NL_out1=[ NL_out1 NL_dn];
        NL_out2=[ NL_out2 NL_sm];
        NL_out3=[ NL_out3 NL_up];
        
        Jdate_out1=[Jdate_out1 Jdate1];
        Jdate_out2=[Jdate_out2 Jdate2];
        Jdate_out3=[Jdate_out3 Jdate3];

    end


    NL_lower=NL_out1(1:round(length(NL_out1)/500):end);
    Jdate_lower=Jdate_out1(1:round(length(NL_out1)/500):end);
    
    NL_sm=NL_out2(1:round(length(NL_out2)/500):end);
    Jdate_sm=Jdate_out2(1:round(length(NL_out2)/500):end);
    
    NL_glider=NL_out3(1:round(length(NL_out3)/500):end);
    Jdate_glider=Jdate_out3(1:round(length(NL_out3)/500):end);
    
 
    
    [~,I]=sort(Jdate_glider);
    
    NL_lower=NL_lower(I)-10*log10(bw(1));
    Jdate_lower=Jdate_lower(I);
    
    NL_sm=NL_sm(I)-10*log10(bw(2));
    Jdate_sm=Jdate_sm(I);
    
    NL_glider=NL_glider(I)-10*log10(bw(3));
    Jdate_glider=Jdate_glider(I);
    
    
    figure
    set(findall(gcf,'type','text'),'FontSize',12,'FontName','TimesNewRoman')
    subplot(2,2,1)
    boxplot([NL_lower(20:end)', NL_sm(20:end)', NL_glider(20:end)'], 'Labels',{'32.2 kHz','40.6 kHz','51.3kHz'}, 'Notch', 'on')
    title('Noise Levels')
    ylabel('dB_{ASL}')
    
    subplot(2,2,2)
    hist(NL_glider(20:end)-NL_sm(20:end))
    title('Histogram of NL Differences')
    xlabel('Difference dB_{ASL}')
    
    subplot(2,2,3)
    plot(Jdate_lower(20:end), smooth(NL_lower(20:end),20),'c', Jdate_glider(20:end), smooth(NL_glider(20:end),20),'b', Jdate_sm(20:end), smooth(NL_sm(20:end), 20), 'r')
    datetick
    title('Noise Level')
    ylabel('dB_{ASL}')
    legend({'32.2 kHz','40.6 kHz','51.3kHz'})
    
    subplot(2,2,4)
    plot(Jdate_glider(20:end), smooth((NL_glider(20:end)-NL_sm(20:end))), 'k')
    title('40.6 kHz- 51.3kHz')
    ylabel('NL Diff dB_{ASL}')
    datetick
    
    pname={'Glider NL Comparison'};
    print(cell2mat(pname),'-dpng','-r250')

