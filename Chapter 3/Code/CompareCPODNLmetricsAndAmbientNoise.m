close all; clear all; clc;
cd('W:\KJP PHD\3-Detection Function\Propagation Model\Code')
floc='W:\KJP PHD\3-Detection Function\Propagation Model\BathData\';
bath_floc='Z:\PropModelBathData\';
set(0, 'DefaulttextInterpreter', 'none')
% Setup noise level directory
NL_Loc='Z:\2013SM_output\2013_SM_BinaryData';


%% Deployment lat/long
Dep_loc=table(['Arb_10'; 'Arb_15';'Arb_05';...
    'Cro_10';	'Cro_15';	'Cro_05';...
    'Cru_10';	'Cru_15';	'Cru_05';...
    'Fra_10';	'Fra_15';	'Fra_05';...
    'Hel_10';	'Hel_15';	'Hel_05';...
    'Lat_10';	'Lat_15';	'Lat_05';...
    'SpB_10';	'SpB_15';	'SpB_05';...
    'Stb_10';	'Stb_15';	'Stb_05';...
    'StA_10';	'StA_15';	'StA_05';...
    'Sto_10';	'Sto_15';	'Sto_05'],...
[56.49963; 56.45947; 56.55413; 57.68911; 57.70653; 57.67517;...
    57.38005; 57.37742; 57.37994; 57.77086; 57.84925; 57.71136;...
    58.00499; 57.97563; 58.05338; 58.22928; 58.18665; 58.26933;...
    57.74148; 57.78529; 57.68989; 55.96354; 56.03334; 55.9292;...
    56.25781; 56.29021; 56.26576; 56.95927; 56.98069; 56.94688],...
   [-2.38019;	-2.29861;	-2.48347;	-3.88155;	-3.81024;	-3.98799;...
    -1.73735;	-1.61812;	-1.82851;	-2.13963;	-2.08936;	-2.12987;...
    -3.61127;	-3.5361;	-3.71525;	-3.2065;	-3.13512;	-3.31806;...
    -3.03874;	-3.05898;	-3.06196;	-2.16201;	-2.07546;	-2.17725;...
    -2.49871;	-2.43311;	-2.571429;	-2.11361;	-2.02194;	-2.17688],...
    zeros(30,1), zeros(30,1), 'VariableNames',{'Loc', 'Lat', 'Lon', 'UTMX', 'UTMY'});

Dep_loc.Var6(14,1) = 7;
Dep_loc.Properties.VariableNames{6} = 'SMNumber';

FBands=readtable('W:\KJP PHD\SM2M Processing\PAMGuardCenterFrequencyNLBands.csv');
df=(FBands.CenterFreq*2^(1/6))-(FBands.CenterFreq*2^(-1/6));
FBands=[FBands table(df, 'VariableNames', {'df'})];
fband_idx=nearest2(FBands.CenterFreq, 41000);
%% Do the analysis for all years 
years=[2013 2014];
nyears= length(years);



    Year=years(1);
    
    %% Setup noise level directory 2013

        Parent_loc='Z:\2013SM_output\2013_SM_BinaryData';
        SM_dir=dir('Z:\2013SM_output\2013_SM_BinaryData\SM2M*');
        SM_nums=find(vertcat(SM_dir.isdir));
        SM_dir(find(strcmpi('.',cellstr({SM_dir.name})) |...
        strcmpi('..',cellstr({SM_dir.name}))))=[]; 

    % 2013 Deployment SM2M numbers
        Dep_loc.SMNumber(12) = 4;
        Dep_loc.SMNumber(30) = 12;
        Dep_loc.SMNumber(24) = 3;
        Dep_loc.SMNumber(25) = 9;
        Dep_loc.SMNumber(19) = 10;
        Dep_loc.SMNumber(5) = 15;
        Dep_loc.SMNumber(18) = 13;
        Dep_loc.SMNumber(9) = 6;
        Dep_loc.SMNumber(1) = 5;
  


    [Dep_loc.UTMX, Dep_loc.UTMY] = deg2utm(Dep_loc.Lat,Dep_loc.Lon);
    %load nltxt using the import function, read date as a text
    cpod_time=datenum(nltxt.ChunkEnd, 'dd/mm/yyyy HH:MM');
    nltxt=[nltxt, table(cpod_time, 'VariableNames', {'MatlabDate'})];
%% 

for ii=8 %ii=1:length(SM_nums)

    NL_out=[];
    Jdate_out=[];
    filedays=dir(fullfile([NL_Loc '\' SM_dir(SM_nums(ii)).name], '20*'));
    % List of days for which there are binary files
    
    
    
    for jj=1:length(filedays)-1
        
        fdate=[NL_Loc '\' SM_dir(ii).name '\' filedays(jj).name];
        [NL Jdate]=getMeanMinutePGNL(fdate, fband_idx);
        NL_out=[ NL_out NL];
        Jdate_out=[Jdate_out Jdate];
    end
     
    % Trim the CPOD data
    min(Jdate_out)
    
    CPOD_dates=STOminWNoise.min
    CpodNose=STOminWNoise.
    Jdate_trim=STOminWNoise.min(find(STOminWNoise.min>=min(Jdate_out))); 
    NL_trim=NL(idx1:idx2)
    

%556 1996
end

%% Plot the raw data
% figure
% subplot(2,1,1)
% plot(Jdate_out, smooth(NL_out,60))
% xlim([datenum('25/07/2013', 'dd/mm/yyyy') datenum('02/12/2013', 'dd/mm/yyyy')])
% datetick('x','dd-mmm','keeplimits', 'keepticks')
% ylabel('1/3RD Octave Energy Pa re 1upa')
% 
% subplot(2,1,2)
% plot(cpod_time, smooth(nltxt.MMM, 60))
% xlim([datenum('25/07/2013', 'dd/mm/yyyy') datenum('02/12/2013', 'dd/mm/yyyy')])
% datetick('x','dd-mmm','keeplimits', 'keepticks')
% ylabel('Random CPOD noise Metric')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross Correlation %

% Zero Pad the SM data

SMtable=table(Jdate_out', NL_out', cellstr(datestr(Jdate_out, 'dd/mm/yyyy HH:MM')),...
    'VariableNames', {'MatlabDate', 'SM_NL', 'ChunkEnd'});

C=outerjoin(nltxt, SMtable,'Keys','ChunkEnd');
C = sortrows(C,'MatlabDate_nltxt','ascend');

figure
subplot(2,1,1)
plot(C.MatlabDate_nltxt, C.SM_NL)
xlim([datenum('10/06/2013', 'dd/mm/yyyy') datenum('12/12/2013', 'dd/mm/yyyy')])
datetick('x','dd-mmm','keeplimits', 'keepticks')
ylabel('1/3RD Octave Energy Pa re 1upa')

subplot(2,1,2)
plot(C.MatlabDate_nltxt, C.MMM)
xlim([datenum('10/06/2013', 'dd/mm/yyyy') datenum('12/12/2013', 'dd/mm/yyyy')])
datetick('x','dd-mmm','keeplimits', 'keepticks')
ylabel('Mean MultiPath Minema')



TempC=C;
idx=find(TempC.SM_NL<10);
NewDat=datasample(TempC.SM_NL(find(TempC.SM_NL>0)),...
    length(find(TempC.SM_NL<10)));
TempC.SM_NL(idx)=NewDat;

idx=find(TempC.MMM==0);
NewDat=datasample(TempC.MMM(find(TempC.MMM>0)),...
    length(find(TempC.MMM==00)));
TempC.MMM(idx)=NewDat;


[acor,lag] = xcorr(C.SM_NL,C.MMM);
figure
subplot(3,1,1)
plot(C.MatlabDate_nltxt,C.MMM)
datetick('x','dd-mmm','keeplimits', 'keepticks')
ylabel('Mean Multipath Minima')


subplot(3,1,2)
plot(C.MatlabDate_nltxt,C.SM_NL)
datetick('x','dd-mmm','keeplimits', 'keepticks')
ylabel('SM Measured Noise Level')

subplot(3,1,3)
plot(lag,acor)
ylabel('Autocorrelation MMM & SM noise levels')


% 
% [~,I] = max(abs(acor));
% lagDiff = lag(I)
% timeDiff = lagDiff
% a3 = gca;
% a3.XTick = sort([-3000:1000:3000 lagDiff]);
%% Try only looking at high noise levels

 C(find(C.SM_NL<103),:)=[];
 temtable1=C;
 temtable2=C;
 temtable3=C;
 
temtable1(find(temtable1.meanClkPasq<1),:)=[];
temtable2(find(temtable2.MMM<1),:)=[];
temtable3(find(temtable3.SumPasqm<1),:)=[];

figure
subplot(3,1,1)
scatter(temtable1.SM_NL, temtable1.meanClkPasq)
ylabel('Mean  ClkPasq')

subplot(3,1,2)
scatter(temtable2.SM_NL, temtable2.MMM)
ylabel('Mean MultiPath Minima')

subplot(3,1,3)
scatter(temtable3.SM_NL, (temtable3.SumPasqm))
ylabel('SumPasqm')
xlabel('SM Noise Level')


