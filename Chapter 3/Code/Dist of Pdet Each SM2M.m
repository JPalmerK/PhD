% This code relates the probability of detecting a call at a given time to
% the ambient noise level.

% Reads in 1/3rd octave noise levels from 41 kHz band and assumes same
% amount of energy in the octave bands above. 

% Estimates total noise level in the energy band of the echolcation clicks

% Compares that noise level to the detection probability and area monitored
% calcualted using the Area vs 40khz NL code

close all; clear all; clc;
cd('W:\KJP PHD\3-Detection Function\Propagation Model\Code')
floc='W:\KJP PHD\3-Detection Function\Propagation Model\BathData\';
bath_floc='Z:\PropModelBathData\';
set(0, 'DefaulttextInterpreter', 'none')

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
%% Load the CSV files with the NL as a fucntion of PDET


    csv_loc='W:\KJP PHD\3-Detection Function\Propagation Model\Area_vs_NL\';

    CSV_names=dir([csv_loc '*Tao*']);
    
    % Probability of detecting something within 5km
    Area_1=readtable([csv_loc, CSV_names(1).name], 'Delimiter',',', 'ReadRowNames', 1);
    Area_2=readtable([csv_loc, CSV_names(2).name], 'Delimiter',',', 'ReadRowNames', 1);
    Area_3=readtable([csv_loc, CSV_names(3).name], 'Delimiter',',', 'ReadRowNames', 1);
    
    Range_1=readtable([csv_loc, CSV_names(4).name], 'Delimiter',',', 'ReadRowNames', 1);
    Range_2=readtable([csv_loc, CSV_names(5).name], 'Delimiter',',', 'ReadRowNames', 1);
    Range_3=readtable([csv_loc, CSV_names(6).name], 'Delimiter',',', 'ReadRowNames', 1);
    
    Pdet_1=readtable([csv_loc, CSV_names(7).name], 'Delimiter',',', 'ReadRowNames', 1);
    Pdet_2=readtable([csv_loc, CSV_names(8).name], 'Delimiter',',', 'ReadRowNames', 1);
    Pdet_3=readtable([csv_loc, CSV_names(9).name], 'Delimiter',',', 'ReadRowNames', 1);
    

    
    NL_vals=[70:.1:150]; % from Area vs 40khz NL  30/11/2016
    
    %% Make Some Figures
    
    % Detection Probability
    figure
    subplot(3,1,1)
    plot(NL_vals(2:end),table2array(Pdet_1(:,2:end))')
    title(['SNR Thresh = ', (CSV_names(7).name(19:20)), ' dB, \tau_e = ',...
        (CSV_names(7).name(25:30)), ' sec'],'Interpreter', 'TeX')
    ylabel('P_d_e_t', 'Interpreter', 'TeX')
    subplot(3,1,2)
    plot(NL_vals(2:end),table2array(Pdet_2(:,2:end))')
    title(['SNR Thresh = ', (CSV_names(8).name(19:20)), ' dB, \tau_e = ',...
        (CSV_names(8).name(25:29)), ' sec'],'Interpreter', 'TeX')
    ylabel('P_d_e_t', 'Interpreter', 'TeX')
    subplot(3,1,3)
    plot(NL_vals(2:end),table2array(Pdet_3(:,2:end))')
    title(['SNR Thresh = ', (CSV_names(9).name(19:20)), ' dB, \tau_e = ',...
        (CSV_names(9).name(25:28)), ' sec'],'Interpreter', 'TeX')
    xlabel('Noise Level (dB re: 1\muPa)','Interpreter', 'TeX')
    ylabel('P_d_e_t', 'Interpreter', 'TeX')
    
    % Area Monitored
    figure
    subplot(3,1,1)
    plot(NL_vals(2:end),table2array(Area_1(:,2:end))./1000000')
    title(['SNR Thresh = ', (CSV_names(7).name(19:20)), ' dB, \tau_e = ',...
        (CSV_names(7).name(25:30)), ' sec'],'Interpreter', 'TeX')
    ylabel('Area km^2', 'Interpreter', 'TeX')
    
    subplot(3,1,2)
    plot(NL_vals(2:end),table2array(Area_2(:,2:end))./1000000')
    title(['SNR Thresh = ', (CSV_names(8).name(19:20)), ' dB, \tau_e = ',...
        (CSV_names(8).name(25:29)), ' sec'],'Interpreter', 'TeX')
   ylabel('Area km^2', 'Interpreter', 'TeX')
   
    subplot(3,1,3)
    plot(NL_vals(2:end),table2array(Area_3(:,2:end))./1000000')
    title(['SNR Thresh = ', (CSV_names(9).name(19:20)), ' dB, \tau_e = ',...
        (CSV_names(9).name(25:28)), ' sec'],'Interpreter', 'TeX')
    xlabel('Noise Level (dB re: 1\muPa)','Interpreter', 'TeX')
    ylabel('Area km^2', 'Interpreter', 'TeX')
    

%% Do the analysis for all years 
years=[2013 2014 ];

%2015 different NL values :(
nyears= length(years);

for ll=1:nyears

    Year=years(ll);
    
    %% Setup noise level directory 2013
    if Year==2013
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
    end

    if Year==2014
        % Setup noise level directory 2014
        Parent_loc='Z:\2014SM_output\2014_SM_BinaryData';
        SM_dir=dir('Z:\2014SM_output\2014_SM_BinaryData\SM2M*');
        SM_nums=find(vertcat(SM_dir.isdir));
        SM_dir(find(strcmpi('.',cellstr({SM_dir.name})) |...
        strcmpi('..',cellstr({SM_dir.name}))))=[];

        % 2014 Deployment SM2M numbers
        Dep_loc.SMNumber(1) = 7;
        Dep_loc.SMNumber(5) = 13;
        Dep_loc.SMNumber(9) = 3;
        Dep_loc.SMNumber(12) = 12;
        Dep_loc.SMNumber(14) = 5;
        Dep_loc.SMNumber(18) = 4;
        Dep_loc.SMNumber(19) = 16;
        Dep_loc.SMNumber(24) = 15;
        Dep_loc.SMNumber(25) = 1;
        Dep_loc.SMNumber(30) = 6;
    end
    
    if Year==2015
        Parent_loc='Z:\2015SM_output';
       
       fnames_2015=ls(Parent_loc);
        fnames_2015(find(strcmpi('.',cellstr(fnames_2015)) |...
        strcmpi('..',cellstr(fnames_2015))),:)=[];
    
    
    end


    [Dep_loc.UTMX, Dep_loc.UTMY] = deg2utm(Dep_loc.Lat,Dep_loc.Lon);




%% 
   
    for ii=1:length(SM_nums)

        NL_out=[];
        Jdate_out=[];
        Pdet_out1=[];
        Pdet_out2=[];
        Pdet_out3=[];
        
        Area_out1=[];
        Area_out2=[];
        Area_out3=[];
        
        MaxRange_out1=[];
        MaxRange_out2=[];
        MaxRange_out3=[];
        
        % If the year isn't 2015 extract the median hourly noise levels
        if Year ~= 2015
            filedays=dir(fullfile([Parent_loc '\' SM_dir(SM_nums(ii)).name], '20*'));
                % List of days for which there are binary files
            for jj=1:length(filedays)-1
                fdate=[Parent_loc '\' SM_dir(ii).name '\' filedays(jj).name];
                % Get the median hourly noise level
                [NL Jdate]=getMedianHourlyPGNL(fdate,fband_idx);
                NL_out=[ NL_out; NL];
                Jdate_out=[Jdate_out Jdate];
            end

        else
            fname=fullfile([Parent_loc '\' fnames_2015(ii,:)])
            NL_dat=readtable(fname, 'Delimiter',',','ReadVariableNames',false);
            NL_dat = table2array(array2table(table2array(NL_dat).'));
            Jdate_out=NL_dat(2:end,1);
            NL_out=NL_dat(2:end,2);
        end
        

        % Estimate total energy in the 20-120kHz noise band (after meeting
        % with Doug 07/12/2016)
        NL_out=NL_out-10*log10(df(end))+10*log10(100000);
        
        if Year ~= 2015
            loc_name=table2array(Dep_loc(find(Dep_loc.SMNumber== str2num(SM_dir(ii).name(6:7))),1)) 
        else
            loc_name= fname(32:37)
        end
        % Reference the NL to the Pdet value
        for kk=1:length(NL_out)
            NL_idx=nearest2(NL_vals, NL_out(kk));
            Pdet_out1(kk)=Pdet_1{loc_name, NL_idx+1};
            Pdet_out2(kk)=Pdet_2{loc_name, NL_idx+1};
            Pdet_out3(kk)=Pdet_3{loc_name, NL_idx+1};
            
            
            Area_out1(kk)=Area_1{loc_name,NL_idx+1};
            Area_out2(kk)=Area_2{loc_name,NL_idx+2};
            Area_out3(kk)=Area_3{loc_name,NL_idx+3};
            
            MaxRange_out1(kk)=Range_1{loc_name,NL_idx+1};
            MaxRange_out2(kk)=Range_2{loc_name,NL_idx+2};
            MaxRange_out3(kk)=Range_3{loc_name,NL_idx+3};
        end
       
        
        
        % Create a table for the output
             Pdetsite=table(Jdate_out', NL_out, repmat(cellstr(loc_name),length(NL_out),1),...
                                            Pdet_out1', Pdet_out2', Pdet_out3',...
                                            Area_out1', Area_out2', Area_out3',...
                                            MaxRange_out1',MaxRange_out2',MaxRange_out3',...
                                            repmat(Year,length(NL_out),1),...
                       'VariableNames', {'MatlabDate', 'MedianNoiseLevel','DepLoc',...
                                        'MedianPdet1', 'MedianPdet2', 'MedianPdet3',...
                                        'MedianArea1', 'MedianArea2', 'MedianArea3',...
                                        'MedianRange1', 'MedianRange2','MedianRange3',...
                                         'Year'});
        
        
        % Make and write the Pdet table
        csv_name=['W:\KJP PHD\4-Bayesian Habitat Use\Pdet at Time\' num2str(Year) '\' loc_name '_NLandPdetVals.csv'];
        writetable(Pdetsite, csv_name, 'Delimiter', ',');
        
    
    end
end