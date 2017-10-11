function make_3dcontourSNR(bath_grid, TL_grids, HydXYZ, GridX, GridY, NL, SL)

    % This function makes a 3d contour map of the Transmission loss grids
    
    % Create an empty bath grid
    SNR_grid=zeros(size(bath_grid)); 
    
    nmapx=length(TL_grids)
    
    Z=[];Y=[];X=[]; IDX=[];
    map_angles=[0:floor(360/nmapx):360];
    
        for jj=1:length(TL_grids)
        
        
        x=cosd((360/nmapx)*jj)*TL_grids(jj).range;
        y=sind((360/nmapx)*jj)*TL_grids(jj).range;
        z=SL-TL_grids(jj).meanp(1:length(y))-NL;
        
%         % Help with map scaling, anything with an SNR <0 won't be heard and
%         % anything with an SNR>10 will definitely be heard
%         z(z<0)=0;
%         z(z>=10)=10;
        Z=[Z;z']; Y=[Y;y']; X=[X; x'];
        
        x_grid=round(x/GridX)+HydXYZ(1);
        y_grid=round(y/GridY)+HydXYZ(2);
        

        idx=sub2ind(size(bath_grid), y_grid, x_grid);
        IDX=[IDX idx];

    
        end
    % Set SNR grid values
    SNR_grid(IDX)=Z;
    SNR_grid(SNR_grid<0)=0;
%     SNR_grid(SNR_grid==0)=nan;
    SNR_grid=flipud(SNR_grid);
   
    
    
    contourf(bath_grid);
    hold on
    contourf(flipud(SNR_grid));
    
        % Create a 3d surface of what the detection function might look like
    surf_from_scatter(X,Y,Z)
    xlabel('Distance from sensor');
    ylabel('Distane from sensor');
    zlabel('Projected SNR (dB) at Water Surface');
    title('Projected SNR From Animal Producing a 120db Call in 80db Ambient Noise')
    

end