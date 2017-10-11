%% Listening Range as a function of SNR %%

% This code specifically looks at the projected SNR of the source locations
% within the Cromarty firth. The prupose of this code is to link the
% propagation modelling to the detectability as a function of SNR. 

% Currently this code uses the sonar equation and hypothesized detector
% characteristics to relate the received noise level at the sensor to
% transmisison loss modelling to, ultimately, the area monitored by the
% sensor

% Since the ultimate goal of this project is to look ONLY at animals at or
% near the surface (such as were viewed from the South Suitor in Cromarty)
% the depth of the receivers is limited to 0-1.5m waterdepth. 

% Moreover, unlike the ESC prop model code- here I am using a much finer
% scale model to capture the complexities of the suitor channel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 1- Initilizations
% Clear the works space and set the file directories of the bathymetric
% data (downloaded from somewhere), 

close all; clear all; clc;
cd('W:\KJP PHD\3-Detection Function\Propagation Model\Code')
bath_floc='Z:\PropModelBathData';
set(0, 'DefaulttextInterpreter', 'none')

% Number of radial lines to draw
nmapx=72;

% Frequency to run the bellhop model at
freq=7000;

% maximum distance to investigate the transmission loss (km)
rmax=3;

% Load the Cromarty Bathymetry Grid
[bath_grid, A]=getbathgrid(bath_floc, ['Cro_05.asc']);



%% Deployment lat/long
% deployment location of suitors d
Dep_loc=table(['Cro_05'; 'Cro_05';'Cro_05';...
    'Cro_05';	'Cro_05';	'Cro_05';...
    'Cro_05'],...
[57.683; 56.45947; 56.55413; 57.68911; 57.70653; 57.67517;...
    57.38005],...
   [ -3.988104;	-2.29861;	-2.48347;	-3.88155;	-3.81024;	-3.98799;...
    -1.73735],...
    zeros(7,1), zeros(7,1), ones(7,1)*15, 'VariableNames',{'Loc', 'Lat', 'Lon', 'UTMX', 'UTMY', 'SMNumber'});


  [Dep_loc.UTMX, Dep_loc.UTMY] = deg2utm(Dep_loc.Lat,Dep_loc.Lon);


%% Load the Bathymetry File or create it
  
  for ii=1:1 %height(Dep_loc) 
      
     % Set the hydrophone in XYZ space
     HydXYZ=[nearest2(A.Lat, Dep_loc.UTMX(ii,:)) nearest2(A.Lon, Dep_loc.UTMY(ii,:))];
     HydXYZ(3)=bath_grid(HydXYZ(2), HydXYZ(1));
     
     

     % Get the grid spacing
     GridX=A.meta{9,2};
     GridY=A.meta{10,2};

    % Set or determine the name of the transmission loss grid
    TL_fname=['W:\KJP PHD\3-Detection Function\Propagation Model\TL_Grids\TL_'...
    'Sut' '_0' num2str(ii) '_' num2str(freq/1000) 'kHz' ];
      
    % Look for the transmission loss grid
    % If the file exists, load it if not create it
    if exist([TL_fname '.mat'])>0
        load(TL_fname)
    
    else 
        sprintf('Creating new Bathymetry Grid')
        % Call out the deployment lat and long
        lat=Dep_loc.Lat(ii);
        lon=Dep_loc.Lon(ii);
        
        % Run the transmission loss model
        % This produces mean values for attenuation loss at each distance 
        % from the sensor to a depth of 1.5 meters, approximately the depth
        % of animals observed at the surface
        [TL_grids]=AttenuationAtSurface_Bellhop(nmapx, bath_grid, GridX, GridY,...
        HydXYZ, lat, lon, freq, rmax)
     
        % Save the grid for next time
        save(TL_fname, 'TL_grids');
    end
 
    
    % Assume a Sourece Level (dB)
    SL=120;
    % Noise Level (dB)- Assume for now, measure later
    NL=60;
    

    
    % Create an empty bath grid
    SNR_grid=zeros(size(bath_grid));   
    
    % Make a 3d contour plot of the SNR

    
    Z=[];Y=[];X=[]; IDX=[];TL1=[];
    map_angles=[0:floor(360/nmapx):360];
    
    for jj=1:length(TL_grids)
        
        
        x=cosd((360/nmapx)*jj)*TL_grids(jj).range;
        y=sind((360/nmapx)*jj)*TL_grids(jj).range;
        z=SL-TL_grids(jj).meanp(1:length(y))-NL;
        
        % Help with map scaling, anything with an SNR <0 won't be heard and
        % anything with an SNR>10 will definitely be heard
        TL1=[TL1; z'];
        z(z<0)=0;
        z(z>=10)=10;
        Z=[Z;z']; Y=[Y;y']; X=[X; x'];
        
        x_grid=round(x/GridX)+HydXYZ(1);
        y_grid=round(y/GridY)+HydXYZ(2);
        

        idx=sub2ind(size(bath_grid), y_grid, x_grid);
        IDX=[IDX idx];
        % Set SNR grid values
        SNR_grid(idx)=z;
    
    end
    
    SNR_grid(SNR_grid==0)=nan;
    SNR_grid=flipud(SNR_grid);
    
 
    
    % Create a 3d surface of what the detection function might look like
    figure (1)
    surf_from_scatter(X,Y,TL1)
    title('Transmission Loss Water Surface');
    zlabel('dB Loss') 
    hold on
    scatter3(0,0, max(TL1)-20, 'filled', 'w' )
    
    hold off

    figure (2)
    surf_from_scatter(X,Y,Z)
    colorbar
    title(sprintf('Expected SNR (dB) at Water Surface\r SL=120 (dB) NL= 60 (dB)'));
    hold on
    scatter3(0,0, 10.2, 'filled', 'w' )
    hold off
    
    
    figure (3)
    contourf(bath_grid);
    hold on
    contourf(flipud(SNR_grid));
    hold off
   
    
   
    
    
    
  end
  
%% Link the SNR to the Pdet
figure(4)
subplot(1,2,1)
   xx=0:.1:10;
  theta=.7;
   b=4;
   yy=1-exp(-(xx/theta).^-b)
   plot(xx,yy)
   xlabel('SNR')
   ylabel('Pdet')
   title(sprintf('Hypothesized Relationship Between SNR and Pdet\r for the PAMGuard Whistle Moan Detector'))


subplot(1,2,2)
ZZ=exp(-(Z/theta).^-b);
Pdet_grid=zeros(size(bath_grid));
Pdet_grid(IDX)=ZZ;
Pdet_grid(Pdet_grid==0)=nan;
Pdet_grid=(Pdet_grid);
scatter(Z, ZZ)
   xlabel('SNR')
   ylabel('Pdet')
   title(sprintf('Pdet as a Function of SNR'))
   
figure (5)
surf_from_scatter(X,Y,ZZ)
zlabel('Pdet')
title(sprintf('Probability of Detection in Cromarty Channel\rSL=120 NL=60 (dB)'))

% Create a 3d surface of what the detection function might look like
figure
subplot(2,1,1)

xlabel('Distance from sensor');
ylabel('Distane from sensor');
zlabel('Probability of Detection in Cromarty Channel');


subplot(2,1,2)
temp=bath_grid;
temp(temp<0)=0;
temp(temp>1)=.5;
contourf(temp);
colormap gray
cmap=colormap;
cmap=flipud(cmap);
colormap(cmap);


hold on
contourf(Pdet_grid);
xlim([1700 2000])
ylim([560 740])




  
  
  



