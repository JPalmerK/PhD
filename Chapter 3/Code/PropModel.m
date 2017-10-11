%% Propagation Modeling for Detection Functions %%
% This script is designed to look at the variation in detection functions
% between different deployment locations
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


[Dep_loc.UTMX, Dep_loc.UTMY] = deg2utm(Dep_loc.Lat,Dep_loc.Lon);




%% Read in the bathymetric Data %%
% The frequency (hz) to run the bellhop model
freq=5000;
% Maximum range (km) to run the bellhop model
rmax=7;
%number of radial lines to draw
nmapx=20;

for KK=1:30
 % This is the deployment location which may differ from the Bath Grid

fname=Dep_loc.Loc(KK,:);
fin=bath_floc;

% Use this if grids are unprocessed or single
[bath_grid, A]=getbathgrid(fin, [fname '.asc']);


%%
%%%% data were downloaded from http://digimap.edina.ac.uk/%%%
% generally 1 arc second with 5000 me surrounding


% index of the hydrophone in MATLAB space
HydXYZ=[nearest2(A.Lat, Dep_loc.UTMX(KK,:)) nearest2(A.Lon, Dep_loc.UTMY(KK,:))];
HydXYZ(3)=bath_grid(HydXYZ(2), HydXYZ(1));
% % 
% figure
% imagesc(bath_grid);
% set(gca,'YDir','normal')
% hold on
% title([Dep_loc.Loc(KK,:)])
% scatter(HydXYZ(1), HydXYZ(2), 'w', 'filled')
% ax=gca;
% set(ax, 'XtickLabel', round(linspace(0,abs(max(A.Lat)-min(A.Lat))/1000, length(get(ax,'XTick')))));
% set(ax, 'YtickLabel', round(linspace(0,abs(max(A.Lon)-min(A.Lon))/1000, length(get(ax,'YTick')))));
% xlabel('Distance (km)')
% ylabel('Distance (km)')

%% Set up the for the sound propagation model

% grid spacing of the bathymetry map
GridX=A.meta{9,2};
GridY=A.meta{10,2};



lat=Dep_loc.Lat(KK);
lon=Dep_loc.Lon(KK);
[TL_grids]=ESCBellhop_transmission1(nmapx, bath_grid, GridX, GridY,...
     HydXYZ, lat, lon, freq, rmax)
 
 % Save the resulting TL grid
 
  TL_fname=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TL_' Dep_loc.Loc(KK,:) '_'...
    num2str(freq/1000) 'kHz_20rad' ];
save(TL_fname, 'TL_grids');
KK
end

%% Produce Mean Propagation Models For Each Depth
Area_monitored=[];
db_down=100; % Decide db down value
Whistle_SL=[158 169]; %Janik 2000
Click_SL=[218 228];
freq=[40000 10000];
rmax=[5 7];
NL=28;
Min_SNR=3;
% Maximum range (km) to run the bellhop model
tot_area=[];

map_angles=[1:360/nmapx:360];
% db_down=fliplr([50:10:150]);
db_down=90;
for kk=1:2
    SL=Whistle_SL(kk);
    freq=5000;
    for ii=1:30
        TL_fname=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TL_' Dep_loc.Loc(ii,:) '_'...
        num2str(freq/1000) 'kHz' ];
        load(TL_fname)

        fname=Dep_loc.Loc(ii,:);

        % Use this if grids are unprocessed or single
        [bath_grid, A]=getbathgrid(bath_floc, [fname '.asc']);
        GridX=A.meta{9,2};
        GridY=A.meta{10,2};
        % index of the hydrophone in MATLAB space
        HydXYZ=[nearest2(A.Lat, Dep_loc.UTMX(ii,:)) nearest2(A.Lon, Dep_loc.UTMY(ii,:))];
        HydXYZ(3)=bath_grid(HydXYZ(2), HydXYZ(1));
        % 
    %     figure
    %     imagesc(bath_grid);
    %     set(gca,'YDir','normal')
    %     hold on
    %     title([Dep_loc.Loc(ii,:) '_' num2str(freq/100) 'kHz'])
    %     scatter(HydXYZ(1), HydXYZ(2), 'w', 'filled')
    %     ax=gca;
    %     set(ax, 'XtickLabel', round(linspace(0,abs(max(A.Lon)-min(A.Lon))/1000, length(get(ax,'XTick')))));
    %     set(ax, 'YtickLabel', round(linspace(0,abs(max(A.Lat)-min(A.Lat))/1000, length(get(ax,'YTick')))));
    %     xlabel('Distance (km)')
    %     ylabel('Distance (km)')

%         for jj=1:length(db_down)
            [Area_monitored R]=SNRMonitored(TL_grids, Min_SNR, SL, NL);
%           [Area_monitored(jj,ii), R]=DbMonitored(TL_grids, db_down(jj));
            % Plot the Rose Diagram
%             col_idx=(length(jet)/length(db_down)*jj); 
    %          Plot_Rose(R, map_angles, GridX, GridY, HydXYZ, col_idx);
%         end
        tot_area(ii,kk)=Area_monitored;

    end
end










