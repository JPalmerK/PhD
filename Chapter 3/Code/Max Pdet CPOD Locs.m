

close all; clear all; clc;
cd('W:\KJP PHD\3-Detection Function\Propagation Model\Code')
floc='W:\KJP PHD\3-Detection Function\Propagation Model\BathData\';
bath_floc='Z:\PropModelBathData\';
set(0, 'DefaulttextInterpreter', 'none')


%% Set up the Deployment Table
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
Dep_loc.Properties.VariableNames{6} = 'SMNumber2013';

Dep_loc.SMNumber2013(12) = 4;
Dep_loc.SMNumber2013(30) = 12;
Dep_loc.SMNumber2013(24) = 3;
Dep_loc.SMNumber2013(25) = 9;
Dep_loc.SMNumber2013(19) = 10;
Dep_loc.SMNumber2013(5) = 15;
Dep_loc.SMNumber2013(18) = 13;
Dep_loc.SMNumber2013(9) = 6;
Dep_loc.SMNumber2013(1) = 5;


Dep_loc=[Dep_loc array2table(zeros(30,1), 'VariableNames',{'SMNumber2014'})];
Dep_loc.SMNumber2014(1) = 7;
Dep_loc.SMNumber2014(5) = 13;
Dep_loc.SMNumber2014(9) = 3;
Dep_loc.SMNumber2014(12) = 12;
Dep_loc.SMNumber2014(14) = 5;
Dep_loc.SMNumber2014(18) = 4;
Dep_loc.SMNumber2014(19) = 16;
Dep_loc.SMNumber2014(24) = 15;
Dep_loc.SMNumber2014(25) = 1;
Dep_loc.SMNumber2014(30) = 6;



[Dep_loc.UTMX, Dep_loc.UTMY] = deg2utm(Dep_loc.Lat,Dep_loc.Lon);
[Dep_loc.UTMX, Dep_loc.UTMY] = deg2utm(Dep_loc.Lat,Dep_loc.Lon);


MaxPdet2013=zeros(height(Dep_loc),1);
Dep_loc=[Dep_loc table(MaxPdet2013),...
    table(MaxPdet2013, 'VariableNames', {'MaxPdet2014'}),...
    table(MaxPdet2013, 'VariableNames', {'MaxPdet2015'})];

%% Get the Frequency Bands %%
FBands=readtable('W:\KJP PHD\SM2M Processing\PAMGuardCenterFrequencyNLBands.csv');
df=(FBands.CenterFreq*2^(1/6))-(FBands.CenterFreq*2^(-1/6));
FBands=[FBands table(df, 'VariableNames', {'df'})];
fband_idx=nearest2(FBands.CenterFreq, 41000);

csv_loc='W:\KJP PHD\3-Detection Function\Propagation Model\Area_vs_NL\';
csv_nam=[csv_loc 'Pdet_41kHzBand_5km.txt'];
Pdet=readtable(csv_nam, 'Delimiter',',', 'ReadRowNames', 1);

% Noise level values associated with the Pdet table in dB ASL
NL_vals=[0:.1:70];

%% Grab the median Pdet for a given noise level and plug it into the DepLoc table

Years=[2013 2014];


for yy=1:length(Years)

    
    % Setup noise level directory
    Parent_loc=['Z:\', num2str(Years(yy)), 'SM_output\', num2str(Years(yy)), '_SM_BinaryData'];
    SM_dir=dir(['Z:\', num2str(Years(yy)), 'SM_output\', num2str(Years(yy)), '_SM_BinaryData\SM2M*']);
    SM_nums=find(vertcat(SM_dir.isdir));
    SM_dir(find(strcmpi('.',cellstr({SM_dir.name})) |...
    strcmpi('..',cellstr({SM_dir.name}))))=[];

    for ii=1:length(SM_nums)

        NL_out=[];
        Jdate_out=[];
        Pdet_out=[];
        filedays=dir(fullfile([Parent_loc '\' SM_dir(SM_nums(ii)).name], '20*'));
            % List of days for which there are binary files
        for jj=1:length(filedays)-1
            fdate=[Parent_loc '\' SM_dir(ii).name '\' filedays(jj).name];
            %[NL Jdate]=getPGNL(fdate,fband_idx);
            % Get the median hourly noise level
            [NL Jdate]=getMedianHourlyPGNL(fdate,fband_idx);
            NL_out=[ NL_out; NL];
            Jdate_out=[Jdate_out Jdate];
        end

            % Convert NL Out to NL_ASL
            NL_out=NL_out-10*log10(FBands.df(fband_idx));

            % Median Noise Level
            NL_med=median(NL_out);

            % Minimum Noise Level....
            NL_min=min(NL_out);

            % Index of adjacent CPOD units
            
            loc_idx=find(table2array(Dep_loc(:,(5+yy)))== str2num(SM_dir(ii).name(6:7)));
            
            loc_name=table2array(Dep_loc(loc_idx,1));
            mm=cellstr(Dep_loc.Loc);
            CPODidx=find(strncmpi(mm, loc_name, 3));

            % Noise level index nearest to the median noise level
            NL_idx=nearest2(NL_vals, NL_min);

            % Get the detection probability for given noise level
            Loc_names=mm(CPODidx);
            Med_Pdets=Pdet(Loc_names,NL_idx);

            % Plug the Pdet Values back in
            Dep_loc(CPODidx,7+yy)=Med_Pdets;



    end
end
%% Export the Detection table %%

kk={'Abr_10', 'Abr_15', 'Abr_05'};
Dep_loc.Loc(1:3,:)=char(kk');
writetable(Dep_loc, 'W:\KJP PHD\4-Bayesian Habitat Use\Pdet at Time\MaxPdet.csv');




