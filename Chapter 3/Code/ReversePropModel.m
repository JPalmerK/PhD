% This code checks that the values for the prop model (using bellhop)
% roughly match those of the reverse. E.g. This checks that the acoustic
% recriprocity principle is ok in this situation


% Update 18/04/2016
% Calculate the error rate out to a maximum of 1, 3, and 5km

%% Section 1-Initilizations

close all; clear all; clc;
cd('W:\KJP PHD\3-Detection Function\Propagation Model\Code')
floc='W:\KJP PHD\3-Detection Function\Propagation Model\BathData\';
bath_floc='Z:\PropModelBathData\';
set(0, 'DefaulttextInterpreter', 'none')


% Maximum range (km) to run the bellhop model
rmax=7;
%number of radial lines to draw
nmapx=20;

% Assume a Sourece Level (dB)
SL=213; %RMS from zimmer page 89
% SL=SL-log10( 9.2878e+03) % subtract 

% Noise Level (dB)- Assume for now, measure later
% Update 18/04/2016 median noise level seems to be about 77 dB after
% running PAMGuard on source locations

NL=77;
% Minimum SNR for dB monitored
Min_SNR=120;


% Deployment lat/long
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

[Dep_loc.UTMX, Dep_loc.UTMY] = deg2utm(Dep_loc.Lat,Dep_loc.Lon);

Area_monitored=zeros(height(Dep_loc),1);
Area_monitored_rev=zeros(height(Dep_loc),1);

%% 

% The frequency (hz) to run the bellhop model
freq=40000;

for KK=1:30
    fname=Dep_loc.Loc(KK,:);
    fin=bath_floc;

    % Use this if grids are unprocessed or single
    [bath_grid, A]=getbathgrid(fin, [fname '.asc']);
    GridX=A.meta{9,2};
    GridY=A.meta{10,2};
    
    % index of the hydrophone in MATLAB space
    HydXYZ=[nearest2(A.Lat, Dep_loc.UTMX(KK,:)) nearest2(A.Lon, Dep_loc.UTMY(KK,:))];
    HydXYZ(3)=bath_grid(HydXYZ(2), HydXYZ(1));

    
    % Get or set the name of the propmodel
     TL_fname=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TL_' ...
         Dep_loc.Loc(KK,:) '_' num2str(freq/1000) 'kHz_' num2str(nmapx) 'rad' ];
     
    % Create the reverse TL grid
    if exist([TL_fname '.mat'])>=1
         load(TL_fname)
    else
        lat=Dep_loc.Lat(KK);
        lon=Dep_loc.Lon(KK);
        [TL_grids]=ESCBellhop_transmission1(nmapx, bath_grid, GridX, GridY,...
         HydXYZ, lat, lon, freq, rmax);
        TL_grids=rmfield(TL_grids,{'pa'});
        save(TL_fname, 'TL_grids');
    end
     
    
    % Set the reverse TL grid name
    TL_fname_rev=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TLrev_' Dep_loc.Loc(KK,:) '_'...
    num2str(freq/1000) 'kHz' num2str(nmapx) 'rad' ];    
    
    % Create the reverse TL grid
    if exist([TL_fname_rev '.mat'])==0
    
        lat=Dep_loc.Lat(KK);
        lon=Dep_loc.Lon(KK);
        [TL_grids_rev]=ESCReverseBellhop_transmission(TL_grids,nmapx, bath_grid, GridX, GridY,...
         HydXYZ, lat, lon, freq, rmax);
        save(TL_fname_rev, 'TL_grids_rev');
    else
        load(TL_fname_rev);
    end
      
    
    clear  TL_grids_rev
    
    
  KK
end




%% Produce Mean Propagation Models For Each Depth
Area_monitored=[];

Whistle_SL=[158 169]; %Janik 2000
SL=220; % For clicks
NL=91; % Mode noise level in the 41khz band from fband table
RL_thresh=122; %ptp from Dahne  from Dahne paper

R_max=5000; % This is the maximum detection range over which to calculate the differences
% between 

SNR_thresh=10;
tau=5000*10^-6;

map_angles=[1:360/nmapx:360];

freqs=[30000 35000 40000 45000];

for jj=1:length(freqs)
    freq=freqs(jj);
    tot_area=[];
    
    for ii=1:30
    % Get or set the name of the propmodel
     TL_fname=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TL_' ...
         Dep_loc.Loc(ii,:) '_' num2str(freq/1000) 'kHz_' num2str(nmapx) 'rad' ];
        load(TL_fname)
    
    
    % Set the reverse TL grid name
    TL_fname_rev=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TLrev_'...
        Dep_loc.Loc(ii,:) '_' num2str(freq/1000) 'kHz' num2str(nmapx) 'rad' ];        
    
    load(TL_fname_rev)
    


        
%     % Adjust the TL grids if the maximum pressure is different
%     for ll=1:length(TL_grids)
%         P1max=max(max(TL_grids(ll).pa(:,2:end)));
%         P2max=max(max(TL_grids_rev(ll).pa(:,2:end)));
%         Pdiff=P1max-P2max;
% %         TL_grids(ll).pa(:,2:end)+Pdiff;
%         
%     end
        [Area_monitored, R]=ENRMonitored(TL_grids, SNR_thresh, RL_thresh, SL, NL, tau);
        [Area_monitored_rev R]=ENRMonitored(TL_grids_rev, SNR_thresh, RL_thresh, SL, NL, tau);
        
        tot_area(ii,:)=[Area_monitored Area_monitored_rev];
        disp(ii)

    end


% Get the final area monitored difference 
tot_area(:,3)=diff(tot_area,1,2)./tot_area(:,1);
col_name=['PrctErrBellhop' num2str(freq/1000)];
%%
Dep_loc=[Dep_loc table(tot_area(:,3), 'VariableNames', {col_name})];

Dep_loc{:,end}=abs(Dep_loc{:,end});
end
