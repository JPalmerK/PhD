
function [area_monitored, R_70]=DbMonitored(TL_grids, db_down)

for ii=1:length(TL_grids)
    R_70ind= find(TL_grids(ii).meanp(2:end)<=db_down);
    R_70(ii)=TL_grids(ii).range(max(R_70ind));
end

area_monitored=sum(pi*R_70.^2/length(TL_grids));
end