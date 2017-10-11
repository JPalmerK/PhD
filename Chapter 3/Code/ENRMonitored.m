
function [Area_monitored, max_Range]=ENRMonitored(TL_grids, SNR_thresh, RL_thresh, SL, NL, tau, R_max)
% This function calculates the area monitored as a function of ambient
% noise level

% Use ENRMonitored_demoMode.m to make plots of the relationship between 
% tao, SNR thresh and detection range/area


%%
% Input %
% TL_grids- mean transmission loss grids created using PropModel.m and or
%           MakeMeanTLGrid
% SNR_thresh- threshold SNR value above which clicks are detected
% RL_thresh- threshold RL value above which clicks are detected (from
%   Dhane) THIS IS THE PEAK-TO-PEAK VALUE
% NL- ambient noise level (dB re 1upa) octave levels are fine
% R_max- maximum range over which to calculate the area monitored
% Tau- integration period over which the C-POD is looking for drops in zero
% crossing. Not released by the manufacturer.
% SL- peak-to-peak source level

% R max is the maximum range over which to calculate the SNR monitord
% if empty than ignore

if ~exist('R_max', 'var')
    R_max=0;
end


Area_monitored=zeros(size(NL));
R=zeros(1,length(TL_grids));
max_Range=Area_monitored;
    
SLe=SL+10*log10(tau);
    
    
    
% Calculate the area in each ray 
for hh=1:length(NL)

     for ii=1:length(TL_grids)
        % Transmission loss
        TL=TL_grids(ii).meanp(2:end);
             
         % Received level (peak-to-peak)
         RL=SL-TL;

         % SNR values   
         SNR=SL+10*log10(15e-6/tau)-NL(hh)-TL;
         idx1=find(SNR(2:end)<SNR_thresh,1);

         if isempty(idx1)
              idx1=length(SNR);
         end

         % Index of the first range where the detection SNR drops below the
         % RL threshold
         idx2=find(RL(2:end)<RL_thresh,1);

         % Grab the smaller of the two idx values
         idx=min([idx1, idx2]);
              
         %R(ii)=TL_grids(ii).range(idx);
         R(ii)=TL_grids(ii).bath(idx,1);
      end


    Area_monitored(hh)=sum(pi*R.^2/length(TL_grids)); % in meters squared
    max_Range(hh)=max(R);
    sprintf(num2str(hh));

end


    

end