%% Listening Range as a function of SNR %%

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
%% Set up some Constants and Empty Variables %% 

% Make up some Noise Levels
NL=[20:1:200];

% The frequency (hz) to run the bellhop model (5 or 40khz)
freq=[50000];


%% Read in the bathymetric Data %%


% Click SL
SL_click=213; %dB from Zimmer pg 89
Min_SNR=1;


for ff=1:length(freq)
    % place to store the Area
    area_monitored=zeros(length(NL), 30);

    for ii=1:30
      

        % Load the bathymetry grid;
        fname=Dep_loc.Loc(ii,:);

        % Use this if grids are unprocessed or single
        [bath_grid, A]=getbathgrid(bath_floc, [fname '.asc']);
        % index of the hydrophone in MATLAB space
        HydXYZ=[nearest2(A.Lat, Dep_loc.UTMX(ii,:)) nearest2(A.Lon, Dep_loc.UTMY(ii,:))];
        HydXYZ(3)=bath_grid(HydXYZ(2), HydXYZ(1));

        
        TL_fname=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TL_' fname '_'...
        num2str(freq(ff)/1000) 'kHz_20rad'];
        load(TL_fname)
        GridX=A.meta{9,2};
        GridY=A.meta{10,2};

        % Get the Area Monitored
  
        for jj=1:length(NL)      
            % Add rose plot
            [area_monitored(jj,ii), ]=SNRMonitored(TL_grids, Min_SNR, SL_click, NL(jj));       
         
        end
        

    end
    
   
    
    Area=array2table(area_monitored, 'VariableNames', cellstr(Dep_loc.Loc)');
    csv_loc='W:\KJP PHD\3-Detection Function\Propagation Model\Area_vs_NL\';
    csv_nam=['Nlvs' num2str(freq(ff)/1000) 'kHz'];
    writetable(Area, [csv_loc csv_nam])
    
end

figure
plot(area_monitored)
title([num2str(freq(ff)/1000) 'kHz'])
xlabel('Noise Level dB re: 1uPa')
ylabel('Square Km Monitored')
legend(cellstr(Dep_loc.Loc)')
