
function [area_monitored, R]=SNRMonitored(TL_grids, Min_SNR, SL, NL, R_max)
% This function take in the the transmission loss grids and static valuse
% for the mininimum detectable SNR, estimated source level and ambient noise 
% level. The output is the area monitored for a given noise level

% Range is the Range in meters

% R max is the maximum range over which to calculate the SNR monitord
% if empty than ignore

if ~exist('R_max', 'var')
    R_max=0;
end

% SNR=SL-TL-NL;
area_monitored=[];

if ~isfield(TL_grids, 'meanp')
    for ii=1:length(TL_grids)

        rd=0:.5:round(max(TL_grids(ii).bath(:,2))+10);  
        TL_grids(ii).meanp=getmeanp(TL_grids(ii).pa, TL_grids(ii).bath, rd);

    end
end


% Calculate the area monitored as a function of NL


for jj=1:length(NL)
    for ii=1:length(TL_grids)
        
        SNR=SL-TL_grids(ii).meanp(2:end)-NL(jj);
        R_ind= find(SNR>=Min_SNR);
        R_max=nearest2(TL_grids(ii).bath(:,1), R_max/1000);
        
        if ~isempty(R_ind)     
            r_indfin=min(R_max, max(R_ind));
            R(ii)=TL_grids(ii).bath((r_indfin),1); %now in km
        
        else
            R(ii)=0;
        end
    end
    

area_monitored(jj)=sum(pi*R.^2/length(TL_grids)); % in kilometers squared
sprintf(num2str(jj));

end


end