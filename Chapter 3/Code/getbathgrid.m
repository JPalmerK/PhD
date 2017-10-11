function [bath_grid A]=getbathgrid(fin, fname)

if exist(['Z:\PropModelBathData\Functional MAT files\' fname(1:end-4) '.mat'])>0
    load(['Z:\PropModelBathData\Functional MAT files\' fname(1:end-4) '.mat']);

    bath_grid=A.data;
else

    delimiterIn = ' ';
    A = importdata([fin fname],delimiterIn,6);

    for ii=1:length(A.textdata)
    A.meta(ii,:)=strsplit(A.textdata{ii}, ' ');
    end

    bath_grid=A.data;
    A.meta(end+1,1)=cellstr('XLLUTM');
    A.meta(end+1,1)=cellstr('YLLUTM');

    % bottome left corner
    [a b]=deg2utm(str2double(A.meta{4,2}),str2double(A.meta{3,2}));

    % next square up
    [a1 b1]=deg2utm(str2double(A.meta{4,2})+str2double(A.meta{5,2}),...
        str2double(A.meta{3,2})+str2double(A.meta{5,2}));

    % Delta XY in utm
    deltaLat=abs(a-a1);
    deltaLon=abs(b-b1);

    A.meta(7,2)=num2cell(a);
    A.meta(8,2)=num2cell(b);

    A.LatInd(:,1)= str2double(A.meta{5,2})*[1:1800]+str2double(A.meta{3,2});
    A.LatInd(:,2)= str2double(A.meta(5,2))*[1:1800]+str2double(A.meta{4,2});

    A.UTMInd(:,1)=[1:1800]*[a1-a]+a1;
    A.UTMInd(:,2)=[1:1800]*[b1-b]+b1;

    % A.Lat(:,1)=str2double(A.meta{5,2})*[1:1800]+str2double(A.meta{4,2});
    A.meta(9,1)=cellstr('Xspacing'); A.meta(9,2)=num2cell(deltaLon);
    A.meta(10,1)=cellstr('Yspacing'); A.meta(10,2)=num2cell(deltaLat);

    A.Lat=linspace(a, a+(A.meta{10,2}*str2double(A.meta{2,2})), str2double(A.meta{2,2}));
    % A.Lat(:,1)=ones(size(A.Lat(:,1)));


    A.Lon=linspace(b, b+(A.meta{9,2}*str2double(A.meta{1,2})), str2double(A.meta{1,2}));
    % A.Lon(:,1)=A.Lat(:,1);


    A=rmfield(A,'textdata');


% 
%     save([fname(1:end-4)], 'A');
end
end