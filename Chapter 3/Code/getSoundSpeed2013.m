function [sv,T,S]=getSoundSpeed(mon,lat0,lon0)
warning off
%mon=8;
%lat0=41;
%lng0=5;
%    
root0='../climatology/';
%
Leroy = @(T,S,Z,L) ...
        1402.5 + 5*T -(T/4.288).^2 +(T/16.8).^3 ...
        +(1.33 -(T/81.3) +(T/107.2).^2).*S ...
        +(Z/64.1) +(Z/1980.3).^2 -(Z/5155).^3 ...
        +(Z/18519).*(L/45-1)-T.*(Z/10172).^3 ...
        +((T/57.74).^2 +(S/69.93)).*(Z/1000);
    
% 2013 World Ocean Atlas Data. Salinity and Temperature Values downloaded
% for each month. Each file has it's own folder. http://www.nodc.noaa.gov/OC5/woa13/woa13data.html
lfn=sprintf('woa13_decav_t%02dmn04',mon);
if exist([lfn '.mat'],'file')~=2
    % extract data from climatology
    %
    temp=[];
    root1=sprintf('woa13_decav_t%02dmn04/', mon);
    dirs=dir([root0 root1 '*.csv']);
    %
    for ii=1:length(dirs)
        fname=[root0 root1 dirs(ii).name];
        fprintf('loading %s\n',fname)
        temp(ii).dat=xlsread(fname);
    end

    sal=[];
    root1=sprintf('woa13_decav_s%02dmn01/', mon);
    dirs=dir([root0 root1 '*.csv']);
    %
    for ii=1:length(dirs)
        fname=[root0 root1 dirs(ii).name];
        fprintf('loading %s\n',fname)
        sal(ii).dat=xlsread(fname);
    end
    save(lfn,'temp','sal')
end
%
%load local mat file to access data
load(lfn)
% temp=temp.dat;
% sal=sal.dat;
kk=1;
S=nan;
T=nan;
lat01=lat0;
while length(S(~isnan(S)))<3 || length(T(~isnan(T)))<3

%extract temperature and salinity vector
    gd=find(abs(temp.dat(:,1)-lat0)<=1 & abs(temp.dat(:,2)-lon0)<=1);

         lon01=lon0+(.05*kk);
    while isempty(gd)
         lon01=lon01+(.05*kk);
        gd=find(abs(temp.dat(:,1)-lat01)<=1 & abs(temp.dat(:,2)-lon01)<=1);
    end
    
    % temperature
    T=nanmean(temp.dat(gd,3:end),1);
    % depth
    D=[0 temp.dat(1,4:end)];
    % salinity
    gd_sal=find(abs(sal.dat(:,1)-lat01)<=1 & abs(sal.dat(:,2)-lon01)<=1);
    S=nanmean(sal.dat(gd_sal,3:end),1);
    lon01=lon0+(.05)*kk;
    kk=kk*-1;
end
%estimate sound speed



sv=[D',Leroy(T,S,D,lat0)'];
