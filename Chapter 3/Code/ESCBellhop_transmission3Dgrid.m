
function [TL]= ESCBellhop_transmission(nmapx, bath_grid, GridX, GridY,...
     HydXYZ, lat, lon, freq, rmax)
 
% Plot Parameters
ncol=20;
dBs=10;
cmap=gray(ncol);
ang=-89:89; 
% map_angles=[0:floor(360/8):360]; % this has been updated to hopefully
% produce a more accurate estimation for the validation study
map_angles=[0:360];

r=rmax*1000; % Radius from hydrophone from the source in meters
rng=linspace(0, rmax, 501)';
mon=8; %Month
 
%plot parameters
ncol=20;
dBs=10;
cmap=gray(ncol);
 

    
% Get SoundSpeed Profile
     [sv, T, S]=getSoundSpeed2013(mon, round(lat), round(lon));
     sv(isnan(sv(:,2)),:)=[];
 
% For each Map Angle Create a Transmission Loss Grid
map_angles=[0:floor(360/nmapx):360];

  % get soundspeed profile
    [sv, T, S]=getSoundSpeed2013(mon, round(lat), round(lon));
    sv(isnan(sv(:,2)),:)=[];
    
% Create a blank transmission loss structure
TL=struct();

  for ii=1:length(map_angles)-1
          copyfile('C:\Program Files\MATLAB\PAM_BOOK_PROGRAMS\SHDFIL',...
        'W:\KJP PHD\3-Detection Function\Propagation Model\Code')
        
%      % look for the bathymetry file
%      if exist(fname)~=2
    [bty, x, y]=getBottomContour(bath_grid, map_angles(ii),r,GridX, GridY,...
                   HydXYZ(1), HydXYZ(2));
               
    rng=linspace(0, rmax, 501)';

    % truncate the bathymethry if the something pokes above 2 m
    trunk_dist=find(bty(:,2)<=2);
    if length(trunk_dist)>0
        bty=bty(1:trunk_dist(1),:);
        rng=rng(1:trunk_dist(1));
    end
     
    % Set the source depth
    sd=bty(1,2)-1.5;
     
     
    % receiver depths
    rd=0:.5:round(max(bty(:,2))+10);  
    
    %BOT.depth=sd+3; %  BOT.depth needs to be shallower than rd (?)
    BOT.depth=max(bty(:,2));
    BOT.dens=2.4;
    BOT.pSpeed= 1575; 
    BOT.pAtt=0.1; % dB/lambda 
    
  
    % extend sv if necessary
    if rd(end)>BOT.depth && rd(end)>max(sv(:,1))
%         sv(end+1,1)=BOT.depth;
        sv(end+1,1)=rd(end);
        sv(end,2)=sv(end-1,2)+(sv(end,1)-sv(end-1,1))/64.1;
        sv_adj=1;
    else
        sv_adj=0;
        
    end
    
    % Run the Bellhop Model
   [pressure, Pos]=doBellhop(freq, sv, BOT, sd, rd, bty, ang);
    
    % Create the Figures
    pressure(pressure==0)=1e-38;
    RL.pa=abs(squeeze(pressure));
    Db=-20*log10(abs(squeeze(pressure)));
    RL.bath=bty;
    RL.range=linspace(0,rmax*1000,size(Db,2));
    RL.freq=freq;
    RL.depth=rd;
    RL.meanp=getmeanp(RL.pa, bty, rd);
    
    
%     save(sprintf('TL%.0f',freq),'pressure','Pos');
% 
%     figure
%     clf
%     set(gcf,'PaperOrientation','landscape');
%     set(gcf,'PaperPositionMode','auto');
%     set(gcf,'position',[100 300 860 330]);
%     %sv profile -----------------------------
%     ax1=axes('units','pixel','position',[60,40,100,250]);
%     plot(sv(:,2),sv(:,1),'k.-'), axis ij,grid on
%     xlabel('Sound speed [m/s]')
%     ylabel('Depth [m]')
%     ylim(rd([1 end]))
%     %RL image -------------------------------
%     ax2=axes('units','pixel','position',[200,40,610,250]);
%     hold off
%     hi=imagesc(Pos.r.range/1000,Pos.r.depth,RL.db);
%     caxis( [0 ncol*dBs] )
%     set(gca,'yticklabel',[])
%     title('Transmission Loss','interpreter','none')
%     colormap(cmap)
%     xlabel('Range [km]')
%     %
%     hc=colorbar;
%     set(hc,'units','pixel','position',[820 40 10 250])
%     set(get(hc,'title'),'string',' dB')
%     set(hc,'ydir','reverse')
%     %
%     hold on
%     plot(bty(:,1),bty(:,2),'k','linewidth',2)
%     hold off
% %     
%     
 

    if sv_adj==1
        sv(end,:)=[];
    end
    
    if ii==1
        TL=RL;
    else
    TL=[TL, RL];
    end
  
 end
 
 
 
 
 
 
 
 
 end