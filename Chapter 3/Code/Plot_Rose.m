function Plot_Rose(R, map_angles, GridX, GridY, HydXYZ, col_idx)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
rose_width=round(length(map_angles)*.7);

if max(R)<10
    R=R*1000;
end
for ii=1:length(R)
 %%   
    % grid spacing of the bathymetry map
    x1=HydXYZ(1); 
    y1=HydXYZ(2);
    
    x2=x1+cosd(map_angles(ii)).*(R(ii)/GridX);
    y2=y1+sind(map_angles(ii)).*(R(ii)/GridY);
    
    x3=x1+cosd([map_angles(ii)+rose_width])*(R(ii)/GridX);
    y3=y1+sind([map_angles(ii)+rose_width])*(R(ii)/GridY);

%     hold on
    temp=flipud(gray);
    p=patch([x1 x2 x3], [y1 y2 y3], col_idx);
%
end


end

