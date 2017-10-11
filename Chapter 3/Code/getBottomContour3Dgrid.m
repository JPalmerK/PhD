function [bty, x, y]=getBottomContour3Dgrid(bath_grid, grid_idx,r,GridX, GridY,x1,y1)
    
    %%%% Input %%%%
    % bath_grid - bathymetry grid
    % grid_idx - the index of the bathyemetry grid from where to calculate
    %   the transmission loss values
    % GridX GridY - real spacing (in meters) between X and Y spacing
    % x1, y1 location of hydrophone on map space (index)
    % r - maximum distance (in meters) over which to pull the bathymetry
    
    %%%% Output %%%%%
    % bty - batheymetry (range depth) 
    % x - x indexes of the map
    % y- y indexes of the map
    

    x2=x1+cosd(map_angles)*(r/GridX);
    y2=y1+sind(map_angles)*(r/GridY);

    x=round(linspace(x2,x1,501));
    y=round(linspace(y2,y1,501));

    idx=sub2ind(size(bath_grid), y, x)';
    
    bath_contour=-bath_grid(idx);
    bty(:,2)=flipud(bath_contour); %bathymetry depth
    bty(:,1)=linspace(0, r(end), 501)/1000';
    
    % Add a gentle smoothing of 14 m
    bty(:,2)=smooth(bty(:,2),10, 'moving');
    
    
    
%     figure (1)
%     hold on
%     plot(x, y, 'w')
%     hold off

end