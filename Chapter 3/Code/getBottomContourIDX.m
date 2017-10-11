function [idx]=getBottomContourIDX(bath_grid, map_angle ,r,GridX, GridY, x1, y1,n_spaces)

% this fucntion gets the bottom contours for the ranges

    x2=x1+cosd(map_angle)*(r/GridX);
    y2=y1+sind(map_angle)*(r/GridY);

    x=round(linspace(x2,x1,n_spaces));
    y=round(linspace(y2,y1,n_spaces));

    idx=sub2ind(size(bath_grid), y, x)';
    
%     bath_contour=-bath_grid(idx);
%     bty(:,2)=flipud(bath_contour); %bathymetry depth
%     bty(:,1)=linspace(0, r(end), n_spaces)/1000';
    
%     figure (1)
%     hold on
%     plot(x, y, 'w')
%     hold off

end