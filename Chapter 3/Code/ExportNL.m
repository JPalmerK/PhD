
% This code exports calls and exports the median hourly noise levels for
% the SM2M's in 2013. 

% To transform 1/3 octave NL's into subtract 10*log10(BW) % Meeting with
% Doug April 15 2016. Note the factor of 10

% The code calls the getHrlyMedNL file 


close all; clear all; clc;
cd('W:\KJP PHD\3-Detection Function\Propagation Model\Code')
floc='W:\KJP PHD\3-Detection Function\Propagation Model\BathData\';
bath_floc='Z:\PropModelBathData\';
set(0, 'DefaulttextInterpreter', 'none')
% Setup noise level directory
Parent_loc='Z:\2013SM_output\2013_SM_BinaryData';

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
%% Get the f bands for the 1/3 octave filter of choice
FBands=readtable('W:\KJP PHD\SM2M Processing\PAMGuardCenterFrequencyNLBands.csv');


% The frequency (hz) 1/3 octave band to run the bellhop model (5 or 40khz)
freq=45000;
% Calculate the delta frequency for each band 1/3 octave band
df=(FBands.CenterFreq*2^(1/6))-(FBands.CenterFreq*2^(-1/6));
FBands=[FBands table(df, 'VariableNames', {'df'})];

fband_idx=nearest2(FBands.CenterFreq, freq);

% delta frequency in the 40108.92 hz band
dfband=FBands.CenterFreq(fband_idx)

%% Noise levels and Combination 

SM_dir=dir('Z:\2013SM_output\2013_SM_BinaryData\SM2M*');
SM_nums=find(vertcat(SM_dir.isdir));
SM_dir(find(strcmpi('.',cellstr({SM_dir.name})) |...
    strcmpi('..',cellstr({SM_dir.name}))))=[];

for ii=2:length(SM_nums)
  
tic
NL_out=[];
Jdate_out=[];
filedays=dir(fullfile([Parent_loc '\' SM_dir(SM_nums(ii)).name], '20*'));
% List of days for which there are binary files

    for jj=1:length(filedays)-1
        
        fdate=[Parent_loc '\' SM_dir(ii).name '\' filedays(jj).name];
%        [NL Jdate]=getPGNL(fdate,fband_idx);
        % Get the median hourly noise level
        [NL Jdate]=getMedianHourlyPGNL(fdate,fband_idx)
        NL_out=[ NL_out; NL];
        Jdate_out=[Jdate_out Jdate];

    end
    
    % Convert 1/3 band level to mean Sound Pressure Density Spectrum Level
    NL_SPDS=NL_out-(10*log10(FBands.df(fband_idx)));
    
    figure
    hist(NL_SPDS,50)
    toc
    
    NL1=table(NL_SPDS, Jdate_out', cellstr(datestr(Jdate_out)), 'VariableNames',....
        {'NLasl', 'MatlabDate', 'Date'});
    NL1.SM= cellstr(repmat(SM_dir(SM_nums(ii)).name, [height(NL1),1]));
    
    NL1.DepLoc=cellstr(repmat(Dep_loc.Loc(find(Dep_loc.SMNumber==str2num(SM_dir(ii).name(end-1:end))),:), [height(NL1),1]));
    
    
    csv_loc='W:\KJP PHD\3-Detection Function\Propagation Model\Hourly NL 41khz thrdOctBand\';
    csv_nam=[SM_dir(SM_nums(ii)).name '.csv'];
    
    %csvwrite([csv_loc csv_nam], NL_SPDS)
   
    writetable(NL1, [csv_loc csv_nam], 'Delimiter',',','WriteRowNames',false)
   %writetable(Pdet, csv_nam, 'Delimiter',',','WriteRowNames',true);
    
   
           
    end

