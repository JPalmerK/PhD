function surf_from_scatter(x,y,z)

%% Little triangles
% The solution is to use Delaunay triangulation. Let's look at some
% info about the "tri" variable.

x=round(x); y=round(y);
tri = delaunay(x,y);
plot(x,y,'.')

%%
% How many triangles are there?

[r,c] = size(tri);
disp(r)

%% Plot it with TRISURF

h = trisurf(tri, x, y, z);
axis vis3d

%% Clean it up
set(gca, 'zgrid', 'Off')
set(gca, 'xgrid', 'Off')
set(gca, 'ygrid', 'Off')
set(gca, 'color', 'none')
lighting phong
shading interp
colorbar('off')


end