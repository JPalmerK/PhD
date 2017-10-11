%% Listening Range as a function of SNR %%

% This code calculates the unconditional probability of detecting
% an echolcation click within the listening area of each of the ECOMASS
% deployments within the 40kHz 1/3rd octave band

% 2/12/2016 Change threshold to 120 again. The median detection radius was
% 5km for the entire site

% 2/12/2016 Nope, that drops the detection radius down to nothing. Move it
% back to 100dB


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
NL=[70:.1:150]; % Raw noise levels in the 40.1kHz 1/3rd octave levels are
% in the region of 80dB

Mean_AreavsNL=table(cellstr(Dep_loc.Loc), 'VariableNames', {'DepLoc'});
Mean_AreavsNL=[Mean_AreavsNL array2table(zeros(30,length(NL)))];

% Click SL
% SL_click=213; 

SL=220; %pp from Principles of Mar Bio Acs pp 503


% Minimum detection threshold for CPODs
RL_thresh=122; %ptp from Dahne

tau=[1500 5000 8000]*10^-6; % Estimate of the integration period for ENR

SNR_thresh=[1 10 15]; % dB 29/11/2016. Potential SNR thresholds above which clicks will be detected

R_MAX=6000; % Maximum range over which to calculate the probability of detection 

Mean_AreavsNL=table(cellstr(Dep_loc.Loc), 'VariableNames', {'DepLoc'});
Mean_AreavsNL=[Mean_AreavsNL array2table(zeros(30,length(NL)))];

Max_Range=Mean_AreavsNL;
Pdet_41kHzBand_5km=Mean_AreavsNL;

area_monitored=zeros(30, length(NL));
max_Range=area_monitored;

%%

% For each SNR threshold
for jj=1:length(SNR_thresh)
    
    % and each deployment location
    for ii=1:30

                % Load the bathymetry grid
                fname=Dep_loc.Loc(ii,:);

                % Use this if grids are unprocessed or single
                [bath_grid, A]=getbathgrid(bath_floc, [fname '.asc']);
                
                % index of the hydrophone in MATLAB space (Place the
                % Hydrophone on the map- mostly idiot checking)
                HydXYZ=[nearest2(A.Lat, Dep_loc.UTMX(ii,:)) nearest2(A.Lon, Dep_loc.UTMY(ii,:))];
                HydXYZ(3)=bath_grid(HydXYZ(2), HydXYZ(1));

%                  imagesc((bath_grid))
%                  hold on
%                  scatter(HydXYZ(2), HydXYZ(1),'filled', 'w' )

                GridX=A.meta{9,2};
                GridY=A.meta{10,2};
                
                % Load the  transmission loss file (average TL as a function
                % of range and theta fromt the CPOD)
                MeanTL_fname=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TL_' Dep_loc.Loc(ii,:) '_'...
                 '40_1kHz_Oct20rad']
                load(MeanTL_fname)

                % Area monitored as a function of deployment location (TL),
                % detector performance (tau and RL_thresh), source level of
                % the dolphin (SL) and ambient noise level (NL-assume
                % animals are not responding to ambient noise by changing
                % their clicking behaviour) 
                
                [area_monitored(ii,:), max_Range(ii,:)]=ENRMonitored(MeanTL, SNR_thresh(jj),... 
                RL_thresh, SL, NL, tau(jj));

                % Use this code if to investigate the relationship between
                % RL thresh, tao and area monitored- make pretty pictures
                %[Area_monitored, max_Range]=ENRMonitored_demoMode(MeanTL, RL_thresh, SL, 90)
                
%                 
%                 figure
%                 hold on
%                 plot(NL, area_monitored(ii,:))
%                 
%                 title(Dep_loc.Loc(ii,:))
%                 ylim([0 50])
                

    end
    
    Mean_AreavsNL(:, 2:end)=array2table(area_monitored);
    Max_Range(:, 2:end)=array2table(max_Range);
    Pdet_41kHzBand_5km(:, 2:end)=array2table(area_monitored/(pi*R_MAX^2));
    
    for ll=2:width(Mean_AreavsNL)
        Mean_AreavsNL.Properties.VariableNames(ll)= cellstr(['dB_'  strrep(num2str(NL(ll-1),4), '.', '_')]);
        Max_Range.Properties.VariableNames(ll)=cellstr(['dB_'  strrep(num2str(NL(ll-1),4), '.', '_')]);
        Pdet_41kHzBand_5km.Properties.VariableNames(ll)= cellstr(['dB_'  strrep(num2str(NL(ll-1),4), '.', '_')]);    
    end

    
     csv_loc='W:\KJP PHD\3-Detection Function\Propagation Model\Area_vs_NL\';
     csv_nam1=['AreaVsNL_41KhzBand_', num2str(SNR_thresh(jj)), 'SNR_',...
     num2str(tau(jj)), 'Tao.csv'];
     writetable(Mean_AreavsNL, [csv_loc csv_nam1])
        
     csv_nam2=['MaxRangeVsNL_41KhzBand_', num2str(SNR_thresh(jj)), 'SNR_',...
            num2str(tau(jj)), 'Tao.csv'];
     writetable(Max_Range, [csv_loc csv_nam2])
     
     
     csv_nam3=['Pdet_41kHzBand_5km', num2str(SNR_thresh(jj)), 'SNR_',...
     num2str(tau(jj)), 'Tao.csv'];
     writetable(Pdet_41kHzBand_5km, [csv_loc csv_nam3])
    
    
   sprintf(num2str(jj))
end
%%
% figure
% plot(NL,  max_Range./1000)
% title(['Max Range vs NL at tau = ', num2str((tau)),'sec'])
% figure
% plot(NL, area_monitored./1000000)
% title(['Area monitored (km^2) vs NL at tau = ', num2str((tau)),'sec'])

% %%
% Mean_AreavsNL(:, 2:end)=array2table(area_monitored);
% Max_Range(:, 2:end)=array2table(max_Range);
% Pdet_41kHzBand_5km(:, 2:end)=array2table(area_monitored/(pi*5000^2));
% %    figure
% %    hold on
% %    plot(area_monitored')
% %    plot(kk, '*k')
% %    l1=legend([cellstr(num2str(freq'/1000)); {'Weighted mean'}])
% %    mm=get(l1, 'Title')
% %    set(mm, 'String', 'kHz', 'EdgeColor', [0,0,0])
% %    xlabel('NL dB')
% %    ylabel('km^2')
% %    title(cellstr(Dep_loc.Loc(ii,:)))
%    
%  
%     for ii=2:width(Mean_AreavsNL)
%     Mean_AreavsNL.Properties.VariableNames(ii)=cellstr(['dB_'  strrep(num2str(NL(ii-1),3), '.', '_')]);
%     Max_Range.Properties.VariableNames(ii)=cellstr(['dB_'  strrep(num2str(NL(ii-1),3), '.', '_')]);
%     Pdet_41kHzBand_5km.Properties.VariableNames(ii)=cellstr(['dB_'  strrep(num2str(NL(ii-1),3), '.', '_')]);
% end
% 
% %%
% 
% 
%  csv_loc='W:\KJP PHD\3-Detection Function\Propagation Model\Area_vs_NL\';
%  csv_nam=['AreaVsNL_41KhzBand'];
%  writetable(Mean_AreavsNL, [csv_loc csv_nam])
%  writetable(Max_Range, [csv_loc 'MaximumRange41kHzBand'])
%  writetable(Pdet_41kHzBand_5km, [csv_loc 'Pdet_41kHzBand_5km'])
% %%
% 
%  
% idx=find(Dep_loc.SMNumber>0)
% 
% 
% plot(sqrt(area_monitored(idx,:)/pi)')
% xlim([120 150])
% ylabel('Maximum Detection Radius (km)')
% subplot(2,1,2)
% plot((area_monitored(idx,:))')
% xlim([70 200])
% ylabel('Maximum Area Monitored (km)')
% 
% 
% 
% idx=find(Dep_loc.SMNumber>0)
% plot(NL, Mean_AreavsNL{idx, 2:end}')
% legend(cellstr(Dep_loc.Loc(idx,:)))
% title('Area vs NL')
% xlabel('NL')
% ylabel('Area Monitored (km^2)')
% xlim([70 100])
% 
% 
% 
% 
% figure
% subplot(2,1,1)
% plot(Mean_AreavsNL{idx, 2:end}')
% legend(cellstr(Dep_loc.Loc(idx,:)))
% title('Area vs NL')
% xlabel('NL')
% xlim([120 150])
% ylabel('Area Monitored (km^2)')
% 
% 
% subplot(2,1,2)
% plot((area_monitored(idx,:))')
% xlim([120 150])
% ylabel('Detection Radius Monitored (km)')
% 
