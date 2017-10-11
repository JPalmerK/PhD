% This code is used for combining ASC bathymethry grids from the eastern
% scottish coast

[bath_grid, B]=getbathgrid(bath_floc, ['Arb_05a' '.asc']);
[bath_grid, A]=getbathgrid(bath_floc, ['Arb_05c' '.asc']);




C=A;

% if the x's are the same
if str2num(A.meta{3,2})== str2num(B.meta{3,2})
    
    %Determine which is bigger for stacking
    if str2num(A.meta{4,2})< str2num(B.meta{4,2})
        C.data=[B.data; A.data];
        C.UtmInd=[B.UTMInd; A.UTMInd];

        
    else
        C.data=[A.data; B.data]
        C.Lat=sort([A.Lat B.Lat]);
        C.UtmInd=[A.UTMInd; B.UTMInd];
    end
        C.Lat=A.Lat;
        C.Lon=sort([A.Lon B.Lon]);
        C.meta(1,2)={num2str(size(C.data,2))};
        C.meta(2,2)={num2str(size(C.data,1))};
        C.meta(4,2)={num2str(min([str2num(C.meta{4,2}),...
            str2num(B.meta{4,2})]))};
        C.meta(5,2)={num2str(mean([str2num(C.meta{5,2}),...
            str2num(B.meta{5,2})]))};
        [a b]=deg2utm(str2num(C.meta{4,2}),str2num(C.meta{3,2}));
        C.meta{7,2}=num2str(a);
        C.meta{8,2}=num2str(b);
        C.meta{9,2}=mean([C.meta{9,2} B.meta{9,2}]);
        C.meta{10,2}=mean([C.meta{10,2} B.meta{10,2}]);
        
else
     if str2num(A.meta{3,2})< str2num(B.meta{3,2})
         C.data=[A.data B.data];         
         C.UtmInd=[A.UTMInd; B.UTMInd];
         C.meta{7,2}=A.meta{7,2};
     else 
         C.data=[B.data A.data];         
         C.UtmInd=[B.UTMInd; A.UTMInd];
         C.meta{7,2}=B.meta{7,2};
     end
        C.Lon=A.Lon;
        C.Lat=sort([A.Lat B.Lat]);
        C.meta(1,2)={num2str(size(C.data,2))};
        C.meta(2,2)={num2str(size(C.data,1))};
        C.meta(4,2)={num2str(min([str2num(C.meta{4,2}),...
            str2num(B.meta{4,2})]))};
        C.meta(5,2)={num2str(mean([str2num(C.meta{5,2}),...
            str2num(B.meta{5,2})]))};
        [a b]=deg2utm(str2num(C.meta{4,2}),str2num(C.meta{3,2}));
        
        C.meta{8,2}=num2str(b);
        C.meta{9,2}=mean([C.meta{9,2} B.meta{9,2}]);
        C.meta{10,2}=mean([C.meta{10,2} B.meta{10,2}]);
end



 C.data=flipud(C.data)
clear A
BB=C;
save([bath_floc 'Arb_05'], 'A')

% if the y's are the same






