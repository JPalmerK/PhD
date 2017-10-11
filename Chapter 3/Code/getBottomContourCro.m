function [bty, x, y]=getBottomContourCro(bath_grid, map_angles,r,GridX, GridY, x1, y1,n_spaces)


    x2=x1+cosd(map_angles)*(r/GridX);
    y2=y1+sind(map_angles)*(r/GridY);

    x=round(linspace(x2,x1,n_spaces));
    y=round(linspace(y2,y1,n_spaces));

    idx=sub2ind(size(bath_grid), y, x)';
    
    bath_contour=-bath_grid(idx);
    bty(:,2)=flipud(bath_contour); %bathymetry depth
    bty(:,1)=linspace(0, r(end), n_spaces)/1000';
    
%     figure (1)
%     hold on
%     plot(x, y, 'w')
%     hold off

end