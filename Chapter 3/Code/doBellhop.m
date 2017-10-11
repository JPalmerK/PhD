function [pressure,Pos]=doBellhop(freq,sv,BOT,sd,rd,bty,ang)
% ready for TL estimation
% clean up temp files
warning off
delete ENVFIL
delete BTYFIL
delete LOGFIL
delete SHDFIL
warning on

% prepare input files for Bellhop
fid=fopen('ENVFIL','w');
% TITLE
    fprintf(fid,'''PAM_WMXZ''\n');        
% FREQ (Hz)
    fprintf(fid,'%f,\n',freq);            
% NMEDIA
    fprintf(fid,'1,\n');                  
% C-linear, Vacuum, db/lambda, Thorpe
    fprintf(fid,'''CVWT'',\n');           
% ignored, ignored, DEPTH of bottom (m)
    fprintf(fid,'0 0.0 %f \n',sv(end,1)); 
    for jj=1:size(sv,1)
        fprintf(fid,'%f %f /\n',sv(jj,1),sv(jj,2));
    end
%0.0 is bottom roughness in m
    fprintf(fid,'''A*'' 0.0\n');
    fprintf(fid,'%f %f 0.0 %f %f 0 \n',...
        BOT.depth,BOT.pSpeed,BOT.dens,BOT.pAtt);
% the following lines are specific to Bellhop
% NSD 
    fprintf(fid,'%d\n',length(sd));
% SD(1:NSD) (m)
    fprintf(fid,'%f %f /\n',sd([1 end]));
% NRD
    fprintf(fid,'%d\n',length(rd));
% RD(1:NRD) (m)
    fprintf(fid,'%f %f /\n',rd([1 end]));
% NR,
    fprintf(fid,'%d\n',size(bty,1));
% R(1:NR) (km)    
    fprintf(fid,'%f  %f /\n',bty([1 end],1));
% ''R/C/I/S''
    fprintf(fid,'''I''\n'); %incoherent
% NBEAMS
    fprintf(fid,'%d\n',length(ang)); 
% ALPHA1,2 (degrees)
    fprintf(fid,'%f %f /\n',ang([1 end]));
% STEP (m), ZBOX (m), RBOX (km)
    fprintf(fid,'0.0 %f %f,\n',...
        1.01*sv(end,1),1.01*bty(end,1));
fclose(fid);

%%%%%%%%%%%%% Bathy
fid = fopen( 'BTYFIL', 'w' );
fprintf(fid, '''%c'' \n', 'C');
  fprintf(fid,'%d\n',size(bty,1));
  for ii=1:size(bty,1)
      fprintf(fid,'%f %f\n',bty(ii,1),bty(ii,2));
  end
fclose(fid);

%execute Bellhop
tic
!bellhop <ENVFIL >LOGFIL
toc

%read TL
[titleText, plottype, freq, atten, Pos, pressure ]=...
    read_shd_bin( 'SHDFIL');
