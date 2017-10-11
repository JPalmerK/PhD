
function [Area_monitored, max_Range]=ENRMonitored_demoMode(TL_grids, RL_thresh, SL, NL)
% This function calculates the area monitored as a function of ambient
% noise level and the integration period (t_e)

% Takes in TL_grids, RL_thres, SLpp and creates a graph of the relationship
% between the maximum detection range, tau, and SNR_thresh



C={'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colors.
tau=[50:20:10000 ]*10^-6;
SNR_thresh=[0:.1:20];
Tao_snr_ARAEA=zeros(length(SNR_thresh),length(tau));
Tao_snr_RANGE=Tao_snr_ARAEA;
NL_temp=NL; 

   


for kk=1:length(SNR_thresh)
    for jj= 1:length(tau)
        R=[];
        for ii=1:length(TL_grids)
            TL=TL_grids(ii).meanp(2:end);
            RL=SL-TL;
            SNR=202+10*log10(40e-6./tau(jj))-NL_temp-TL;
           % SNR=SL+10*log10(15e-6./tau(jj))-NL_temp-TL;
             idx1=find(SNR(2:end)<SNR_thresh(kk),1);
             idx2=find(RL(2:end)<RL_thresh,1);
             idx=min([idx1, idx2]);
             R(ii)=TL_grids(ii).range(idx);
        end

        Tao_snr_ARAEA(kk,jj)=sum(pi*R.^2/length(TL_grids))/1000000; % in kilometers squared
        Tao_snr_RANGE(kk,jj)=max(R)/1000;
        
    end
      sprintf(num2str(kk))
end
%% Maximum range plot
n_tics=6;
figure
surf(Tao_snr_RANGE)
colormap jet
shading interp
title('Maximum Detection Range')
xlabel('Integration Period (\tau_e) in \musec','Interpreter', 'TeX')
set(gca, 'Xtick', [1:length(tau)/n_tics:length(tau)]-1)
set(gca, 'XtickLabel',  cellstr(num2str(round(quantile(tau/1e-6, [0: 1/n_tics:1])'),4)))

ylabel('SNR Threshold (dB) above Noise')
set(gca, 'Ytick', [1:length(SNR_thresh)/n_tics:length(SNR_thresh)]-1)
set(gca, 'YtickLabel', cellstr(num2str(quantile(SNR_thresh, [0,.2,.4,.6,.8, 1])')))

zlabel('Maximum Detection Range (km)')
%% Maximum Area Plot
figure
surf(Tao_snr_ARAEA)
colormap jet
shading interp
title('Area Detected')

xlabel('Integration Period (\tau_e) in \musec','Interpreter', 'TeX')
set(gca, 'Xtick', [1:length(tau)/n_tics:length(tau)]-1)
set(gca, 'XtickLabel',  cellstr(num2str(round(quantile(tau/1e-6, [0: 1/n_tics:1])'),4)))

ylabel('SNR Threshold (dB) above Noise')
set(gca, 'Ytick', [1:length(SNR_thresh)/n_tics:length(SNR_thresh)]-1)
set(gca, 'YtickLabel', cellstr(num2str(quantile(SNR_thresh, [0,.2,.4,.6,.8, 1])')))

zlabel('Area Monitored (km^2)', 'Interpreter', 'TeX')

%% Maximum range in 2d for picking out  useful variables
figure
imagesc(Tao_snr_RANGE)

title('Maximum Detection Range')
xlabel('Integration Period (\tau_e) in \musec','Interpreter', 'TeX')
set(gca, 'Xtick', [1:length(tau)/n_tics:length(tau)]-1)
set(gca, 'XtickLabel',  cellstr(num2str(round(quantile(tau/1e-6, [0: 1/n_tics:1])'),4)))

ylabel('SNR Threshold (dB) above Noise')
set(gca, 'Ytick', [1:length(SNR_thresh)/n_tics:length(SNR_thresh)]-1)
set(gca, 'YtickLabel', cellstr(num2str(quantile(SNR_thresh, [0,.2,.4,.6,.8, 1])')))

hcb=colorbar
title(hcb,'M')


end