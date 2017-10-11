function [depth]=getbathvalues(A, unitXY, R_length, n)

R_lenght=5000;
thetas=[1:n]*(360/n);
r_dists=zeros(size(A.data));

XYUtm=[A.meta{7,2} A.meta{8,2}];


for ii=1:n
    ang=thetas(ii);
    
    m=cosd(ang)/sind(ang);
    
    y=m*A.UTMInd(:,1)+A.meta{8,2};
    
    
    A.meta{10,2};


end




depth=[];
end