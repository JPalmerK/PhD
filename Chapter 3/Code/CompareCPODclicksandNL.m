% This code compares the median and standard deviation of the SM noise
% values to the number of impulses detected by the adjacent CPOD to
% establish a relationship for the time periods and locations for which
% there were not attached SM units

close all; clear all; clc;
cd('W:\KJP PHD\3-Detection Function\Propagation Model\Code')

% CP1 file directory
CP1_Loc='W:\KJP PHD\CPOD Processing\2013_CP1\CP1 Noise Levels'

% Setup noise level directory
NL_Loc='Z:\2013SM_output\2013_SM_BinaryData';

% Unit numbers available
SM2M_nums=dir(fullfile(NL_Loc, 'SM2M_*'));
SM_dir=dir('Z:\2013SM_output\2013_SM_BinaryData')
SM_nums=find(vertcat(SM_dir.isdir));
SM_dir(find(strcmpi('.',cellstr({SM_dir.name})) | strcmpi('..',cellstr({SM_dir.name}))))=[];


% Load F1 and F2 bands
FBands=readtable('W:\KJP PHD\SM2M Processing\PAMGuardCenterFrequencyNLBands.csv');
FBands.f1=FBands.CenterFreq*2^(-1/6);
FBands.f2=FBands.CenterFreq*2^(1/6);

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
%% Load the NL values

% Will change this for all fbands
fband_idx=height(FBands);


for ii=1:length(SM_nums)
    
    NL_out=[];
    Jdate_out=[];
    filedays=dir(fullfile([NL_Loc '\' SM_dir(SM_nums(ii)).name], '20*'));
    % List of days for which there are binary files
    
    
    
    for jj=1:length(filedays)-1
        
        fdate=[NL_Loc '\' SM_dir(ii).name '\' filedays(jj).name];
        [NL Jdate]=getPGNL(fdate, fband_idx);
        NL_out=[ NL_out NL];
        Jdate_out=[Jdate_out Jdate];
    end
    
        % Get median hourly Noise level
    [med_NL Jdate_hrs stdev_NL]=getHrlyMedNL(Jdate_out, NL_out);
    SMhR=
    
    % Load the corrisponding CP1 file
    Dep_id=Dep_loc.Loc(find(Dep_loc.SMNumber==str2num(SM_dir(SM_nums(ii)).name(end-1:end))),:);
    
    % Add column for matlab Hour
    CP1name=[CP1_Loc '\' Dep_id '.txt'];
    CP1=readtable(CP1name, 'Delimiter', '\t');
    CP1.matlabdate=datenum(CP1.Minute,'dd/mm/yyyy HH:MM')+CP1.microsec/86400000000;
    
    % Subset the data to only frequency values within the band and only
    % time values for which there were also SM recordings
    
    CP1_sub=CP1(CP1.matlabdate>min(Jdate_hrs) & CP1.matlabdate<max(Jdate_hrs),:);
    
    CP1_freqs=CP1_sub(CP1_sub.kHz>=round(FBands.f1(end)/1000) &...
        CP1_sub.kHz<=round(FBands.f2(end)/1000),:);
    
    
    vals=hist(CP1_freqs.matlabdate, Jdate_hrs-(30/1440));
    
   [~, CP1NL_order]=sort(vals);
   CP1NL_orderA=CP1NL_order(CP1NL_order~=length(CP1NL_order));
   
   [~, NL_order]=sort(med_NL);
   
   scatter([1:length(vals)], log10(vals(CP1NL_order)))
   figure
   scatter([1:length(vals)-1],med_NL(CP1NL_orderA))
    
    
end


