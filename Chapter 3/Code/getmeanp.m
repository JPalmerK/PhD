function [meanRL]=getmeanp(RL_Pa, bty, rd)
% This fucntion calculates the transmission loss at all distances from the
% reciver (SM/CPOD).

warning('off')
meanRL=zeros(1,length(bty));

if size(RL_Pa,1)>length(rd)
    RL_Pa=RL_Pa(1:length(rd),:);
end

for ii=1:length(bty)
    bottom_level=bty(ii,2);
    bottom_idx=max(find(rd<=bottom_level));
    meanRL(ii)= -20*log10(mean(RL_Pa(1:bottom_idx,ii)));
end




end