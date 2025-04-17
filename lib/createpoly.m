function [xv, yv] = createpoly(xv, yv)
% Create a convex hull polygon to enclose all 2D points
%
% Author: Hong-Hsi Lee, HLEE84@mgh.harvard.edu

idx = boundary(xv, yv, 0);
xv = xv(idx);
yv = yv(idx);

end