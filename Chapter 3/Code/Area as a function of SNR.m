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
%% Read in the bathymetric Data %%
% The frequency (hz) to run the bellhop model (5 or 40khz)
freq=40000;
% Maximum range (km) to run the bellhop model
rmax=5;
% number of radial lines to draw
nmapx=20;
% map angles

map_angles=[1:360/nmapx:360];

% Click SL
SL_click=213; %dB from Zimmer pg 89
Min_SNR=1;

% Setup noise level directory
Parent_loc='Z:\2013SM_output\2013_SM_BinaryData';

% Unit numbers available
SM2M_nums=dir(fullfile(Parent_loc, 'SM2M_*'));

% Load F1 and F2 bands
FBands=readtable('W:\KJP PHD\SM2M Processing\PAMGuardCenterFrequencyNLBands.csv');

fband_idx=nearest2(FBands.CenterFreq, freq);

%%

SM_dir=dir('Z:\2013SM_output\2013_SM_BinaryData')
SM_nums=find(vertcat(SM_dir.isdir));
SM_dir(find(strcmpi('.',cellstr({SM_dir.name})) | strcmpi('..',cellstr({SM_dir.name}))))=[];

for ii=1:length(SM_nums)
tic
NL_out=[];
Jdate_out=[];
filedays=dir(fullfile([Parent_loc '\' SM_dir(SM_nums(ii)).name], '20*'));
% List of days for which there are binary files

    for jj=1:length(filedays)-1
        
        fdate=[Parent_loc '\' SM_dir(ii).name '\' filedays(jj).name];
        [NL Jdate]=getPGNL(fdate,fband_idx);
        NL_out=[ NL_out NL];
        Jdate_out=[Jdate_out Jdate];
    end
    
    % Get median hourly Noise level
    [med_NL Jdate_hrs]=getHrlyMedNL(Jdate_out, NL_out)

    
    
    
    % Load the bathymetry grid;
    SM_id=str2num(SM_dir(ii).name(end-1:end));
    fidx=find(Dep_loc.SMNumber==SM_id);
    fname=Dep_loc.Loc(fidx,:);
    
    % Use this if grids are unprocessed or single
    [bath_grid, A]=getbathgrid(bath_floc, [fname '.asc']);
    % index of the hydrophone in MATLAB space
    HydXYZ=[nearest2(A.Lat, Dep_loc.UTMX(fidx,:)) nearest2(A.Lon, Dep_loc.UTMY(fidx,:))];
    HydXYZ(3)=bath_grid(HydXYZ(2), HydXYZ(1));
    
    
 
    TL_fname=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TL_' fname '_'...
    num2str(freq/1000) 'kHz' ];
    load(TL_fname)
    GridX=A.meta{9,2};
    GridY=A.meta{10,2};

   figure
    av_name=['Z:\2013data\AVfiles\' fname '_' num2str(freq/1000) 'kHz_'...
        num2str(Min_SNR) 'dB' '.avi'];

    writerObj = VideoWriter(av_name);
    writerObj.FrameRate = 30;
    
    loops=length(NL_out);
%     open(writerObj);
    
    col_scale=round(linspace(min(NL_out), max(NL_out), length(jet)));
    temp=flipud(gray());
    for jj=1:length(NL_out)
        % Create the Basemap
        fh=imagesc(bath_grid);
        set(gca,'YDir','normal')
        hold on
        title([Dep_loc.Loc(fidx,:) '_' num2str(freq/1000) 'kHz' ])
        scatter(HydXYZ(1), HydXYZ(2), 'w', 'filled')
        ax=gca;
        xlim([HydXYZ(1)-300 HydXYZ(1)+300])
        ylim([HydXYZ(2)-300 HydXYZ(2)+300])
        
        xtick=round(linspace(HydXYZ(1)-300, HydXYZ(1)+300, 10));
        ytick=round(linspace(HydXYZ(2)-300, HydXYZ(2)+300, 10));
        
        xticklabels=xtick*A.meta{9,2};
        xticklabels=round((xticklabels-min(xticklabels))/1000);
        
        yticklabels=ytick*A.meta{10,2};
        yticklables=round((yticklabels-min(yticklabels))/1000);        
        
        set(gca, 'Xtick', xtick);
        set(gca, 'Ytick', ytick);
        set(gca, 'XtickLabel', xticklabels);
        set(gca, 'YtickLabel', yticklables);
        xlabel(datestr(Jdate_out(jj)))
        ylabel('Distance (km)')
        col_idx=temp(nearest2(col_scale,NL_out(jj)),:);
        ax.NextPlot = 'replaceChildren';
       
        % Add rose plot
        [area_monitored(jj), R]=SNRMonitored(TL_grids, Min_SNR, SL_click, NL_out(jj));
        Plot_Rose(R, map_angles, GridX, GridY, HydXYZ, col_idx);
        F=getframe(gcf);
        ax.NextPlot = 'replaceChildren';
%         writeVideo(writerObj,F);
        clf()
    end
  close(writerObj);
    
    toc
end





