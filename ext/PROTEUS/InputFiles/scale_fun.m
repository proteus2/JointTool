function [diff] = scale_fun(airfoil1,airfoil2,scale)

airfoil2 = scale*airfoil2;

lmcam1 = sum(sqrt((airfoil1(1:2:end-2)-airfoil1(3:2:end)).^2+(airfoil1(2:2:end-2)-airfoil1(4:2:end)).^2));
lmcam2 = sum(sqrt((airfoil2(1:2:end-2)-airfoil2(3:2:end)).^2+(airfoil2(2:2:end-2)-airfoil2(4:2:end)).^2));

diff = lmcam2-lmcam1;
