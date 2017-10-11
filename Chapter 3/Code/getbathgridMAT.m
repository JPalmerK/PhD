function [bath_grid A]=getbathgrid(fin, name)

load([fin name '.mat']);

bath_grid=A.data;




end