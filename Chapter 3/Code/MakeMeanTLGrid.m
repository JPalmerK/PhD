%% Listening Range as a function of SNR %%

close all; clear all; clc;
cd('W:\KJP PHD\3-Detection Function\Propagation Model\Code')
floc='W:\KJP PHD\3-Detection Function\Propagation Model\BathData\';
bath_floc='Z:\PropModelBathData\';
set(0, 'DefaulttextInterpreter', 'none')
% Load F1 and F2 bands
FBands=readtable('W:\KJP PHD\SM2M Processing\PAMGuardCenterFrequencyNLBands.csv');

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
% Third octave band kHz
Third_oct40khz=[35000 40000 45000];

%% Read in the TL grids and make mean TL for 3rd octave bands

for ii=1:30
    
 
    
    for ff=1:length(Third_oct40khz)
            % Load the bathymetry grid;
            fname=Dep_loc.Loc(ii,:);
            TL_fname=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TL_' fname '_'...
            num2str(Third_oct40khz(ff)/1000) 'kHz_20rad'];
            load(TL_fname)
            
        if ff==1 % If it's the first mean create a place to store the values
          MeanTL=TL_grids;
            
        else
            % If it's not do the sum of the TL grids
            for rr=1:length(TL_grids)
                
                % Some of the little shits are different size
                if size(MeanTL(rr).pa,1)== size(TL_grids(rr).pa,1)
                    MeanTL(rr).pa= MeanTL(rr).pa+TL_grids(rr).pa;
                    
                elseif size(MeanTL(rr).pa,1) < size(TL_grids(rr).pa,1)
                    MeanTL(rr).pa= MeanTL(rr).pa+TL_grids(rr).pa(1:size(MeanTL(rr).pa,1),:);
                    
                elseif size(MeanTL(rr).pa,1) > size(TL_grids(rr).pa,1)
                    MeanTL(rr).pa(1:size(TL_grids(rr).pa,1),:)= MeanTL(rr).pa(1:size(TL_grids(rr).pa,1),:)+TL_grids(rr).pa;
                    MeanTL(rr).pa(size(TL_grids(rr).pa,1):end,:)=MeanTL(rr).pa(size(TL_grids(rr).pa,1):end,:);
                end
            end
        end
    end
    
    % Now divide the MeanTL by the number of bands available in the freq
        for rr=1:length(TL_grids)
                MeanTL(rr).pa= MeanTL(rr).pa./length(Third_oct40khz)
                      
        end  
   
    MeanTL_fname=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TL_' Dep_loc.Loc(ii,:) '_'...
      '40_1kHz_Oct20rad' ];
        save(MeanTL_fname, 'MeanTL');
   
    clear MeanTL
   
 
    
end
%%

for ii=2:width(Mean_AreavsNL)
    Mean_AreavsNL.Properties.VariableNames(ii)=cellstr(['dB_' num2str(NL(ii-1))]);
end
%  csv_loc='W:\KJP PHD\3-Detection Function\Propagation Model\Area_vs_NL\';
%  csv_nam=['AreaVsNL'];
%  writetable(Mean_AreavsNL, [csv_loc csv_nam])

    
idx=find(Dep_loc.SMNumber>0)
plot(Mean_AreavsNL{idx, 2:end}')
legend(cellstr(Dep_loc.Loc(idx,:)))
title('Area vs NL')
xlabel('NL')
ylabel('Area Monitored (km^2)')


