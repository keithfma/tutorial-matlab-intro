% Convert Quabbin Resevoir bathymetric contours to a gridded DEM.

% definitions
xres = 20; % m
yres = 20; % m

% shapefile to (x,y,z) points
[S, A] = shaperead('qrbath/QRBATH_ARC.shp');
%+ allocate
nv = -length(S);
for i = 1:length(S) 
    nv = nv+sum(~isnan(S(i).X));
end
x = nan(nv,1); y = nan(nv,1); z = nan(nv,1);
%+ populate
k = 1;
for i = 1:length(S)
    n = length(S(i).X)-1;
    x(k:k+n-1) = S(i).X(1:n);
    y(k:k+n-1) = S(i).Y(1:n);
    z(k:k+n-1) = A(i).DEPTH_M;
    k = k+n;
end

% grid
[xx,yy] = meshgrid( (min(x)-xres):xres:(max(x)+xres), (min(y)-yres):yres:(max(y)+yres) ); 
zz = griddata(x,y,z,xx,yy);
zz(isnan(zz)) = 0; % outside convex hull is outside lake, so depth=0

% save 
csvwrite('../data/quabbin.csv', zz);


