function [lampar] = lampar_constr(lamparvec,sens)

V1A = lamparvec(1);
V2A = lamparvec(2);
V3A = lamparvec(3);
V4A = lamparvec(4);
V1D = lamparvec(5);
V2D = lamparvec(6);
V3D = lamparvec(7);
V4D = lamparvec(8);

lampar.g1 = 1-2*V1A^2*(1-V3A)-2*V2A^2*(1+V3A)-V3A^2-V4A^2+4*V1A*V2A*V4A;
lampar.g2 = 1- V1A^2 - V2A^2;
lampar.g3 = 1-2*V1D^2*(1-V3D)-2*V2D^2*(1+V3D)-V3D^2-V4D^2+4*V1D*V2D*V4D;
lampar.g4 = 1- V1D^2 - V2D^2;

lampar.g5 = -16 + 32*V3A + 40*V2D^2 + 40*V1D^2 + 16*V4A^2 - 80*V2A*V2D + ...
 72*V2A^2 - 80*V1A*V1D + 72*V1A^2 - 80*V4A*V1D*V2D + 80*V2A*V4A*V1D - ...
 40*V3A*V2D^2 - 120*V3A*V1D^2 - 32*V3A*V4A^2 + 80*V3A*V2A*V2D - ...
 72*V3A*V2A^2 - 32*V3A^3 + 80*V1A*V4A*V2D - 144*V1A*V2A*V4A + ...
 240*V1A*V3A*V1D - 216*V1A^2*V3A - 25*V1D^2*V2D^2 + ...
 50*V2A*V1D^2*V2D - 105*V2A^2*V1D^2 + 160*V3A*V4A*V1D*V2D - ...
 160*V3A*V2A*V4A*V1D - 40*V3A^2*V2D^2 + 120*V3A^2*V1D^2 + ...
 16*V3A^2*V4A^2 + 80*V3A^2*V2A*V2D - 72*V3A^2*V2A^2 + 16*V3A^4 + ...
 50*V1A*V1D*V2D^2 + 20*V1A*V2A*V1D*V2D + 90*V1A*V2A^2*V1D - ...
 160*V1A*V3A*V4A*V2D + 288*V1A*V3A*V2A*V4A - 240*V1A*V3A^2*V1D - ...
 105*V1A^2*V2D^2 + 90*V1A^2*V2A*V2D - 81*V1A^2*V2A^2 + ...
 216*V1A^2*V3A^2 + 50*V4A*V1D*V2D^3 - 150*V2A*V4A*V1D*V2D^2 + ...
 190*V2A^2*V4A*V1D*V2D - 90*V2A^3*V4A*V1D + 50*V3A*V1D^2*V2D^2 - ...
 100*V3A*V2A*V1D^2*V2D + 210*V3A*V2A^2*V1D^2 - 80*V3A^2*V4A*V1D*V2D + ...
 80*V3A^2*V2A*V4A*V1D + 40*V3A^3*V2D^2 - 40*V3A^3*V1D^2 - ...
 80*V3A^3*V2A*V2D + 72*V3A^3*V2A^2 - 50*V1A*V4A*V2D^3 + ...
 190*V1A*V2A*V4A*V2D^2 - 270*V1A*V2A^2*V4A*V2D + 162*V1A*V2A^3*V4A - ...
 100*V1A*V3A*V1D*V2D^2 - 40*V1A*V3A*V2A*V1D*V2D - ...
 180*V1A*V3A*V2A^2*V1D + 80*V1A*V3A^2*V4A*V2D - ...
 144*V1A*V3A^2*V2A*V4A + 80*V1A*V3A^3*V1D + 210*V1A^2*V3A*V2D^2 - ...
 180*V1A^2*V3A*V2A*V2D + 162*V1A^2*V3A*V2A^2 - 72*V1A^2*V3A^3 - ...
 25*V4A^2*V2D^4 + 100*V2A*V4A^2*V2D^3 - 190*V2A^2*V4A^2*V2D^2 + ...
 180*V2A^3*V4A^2*V2D - 81*V2A^4*V4A^2 - 50*V3A*V4A*V1D*V2D^3 + ...
 150*V3A*V2A*V4A*V1D*V2D^2 - 190*V3A*V2A^2*V4A*V1D*V2D + ...
 90*V3A*V2A^3*V4A*V1D - 25*V3A^2*V1D^2*V2D^2 + ...
  50*V3A^2*V2A*V1D^2*V2D - 105*V3A^2*V2A^2*V1D^2 + ...
  50*V1A*V3A*V4A*V2D^3 - 190*V1A*V3A*V2A*V4A*V2D^2 + ...
  270*V1A*V3A*V2A^2*V4A*V2D - 162*V1A*V3A*V2A^3*V4A + ...
  50*V1A*V3A^2*V1D*V2D^2 + 20*V1A*V3A^2*V2A*V1D*V2D + ...
  90*V1A*V3A^2*V2A^2*V1D - 105*V1A^2*V3A^2*V2D^2 + ...
  90*V1A^2*V3A^2*V2A*V2D - 81*V1A^2*V3A^2*V2A^2;

% To define the constraint in line with the other constraints, violation
% <=0 and normalise (32 was selected based on random laminate generation)
lampar.g5 = -lampar.g5/32;

lampar.g6 = 5*(V1A-V1D)^2-2*(1+V3A-2*(V1A)^2);
lampar.g6 = -lampar.g6/10;


if sens == 1
    lampar.dg1 = [-4*V1A*(1-V3A) + 4*V2A*V4A, -4*V2A*(1+V3A)+4*V1A*V4A, 2*V1A^2-2*V2A^2-2*V3A, -2*V4A+4*V1A*V2A, zeros(1,4), 0];
    lampar.dg2 = [-2*V1A,-2*V2A,0,0,zeros(1,4),0];
    lampar.dg3 = [zeros(1,4),-4*V1D*(1-V3D)+4*V2D*V4D,-4*V2D*(1+V3D)+4*V1D*V4D,2*V1D^2-2*V2D^2-2*V3D,-2*V4D+4*V1D*V2D,0];
    lampar.dg4 = [zeros(1,4),-2*V1D,-2*V2D,0,0,0];
    lampar.dg6 = -1/10*[10*(V1A-V1D)+8*V1A,0,-2,0,-10*(V1A-V1D),0,0,0,0];
    lampar.dg5 = -1/32*[-80*V1D+2*72*V1A+80*V4A*V2D-144*V2A*V4A+...
        240*V3A*V1D-2*216*V1A*V3A+...
        50*V1D*V2D^2+20*V2A*V1D*V2D+90*V2A^2*V1D-160*V3A*V4A*V2D+288*V3A*V2A*V4A-240*V3A^2*V1D-2*105*V1A*V2D^2+2*90*V1A*V2A*V2D-2*81*V1A*V2A^2+2*216*V1A*V3A^2-...
        50*V4A*V2D^3+...
        190*V2A*V4A*V2D^2-270*V2A^2*V4A*V2D+162*V2A^3*V4A-100*V3A*V1D*V2D^2-40*V3A*V2A*V1D*V2D-180*V3A*V2A^2*V1D+80*V3A^2*V4A*V2D-144*V3A^2*V2A*V4A+80*V3A^3*V1D+2*210*V1A*V3A*V2D^2-...
        2*180*V1A*V3A*V2A*V2D+2*162*V1A*V3A*V2A^2-2*72*V1A*V3A^3+...
        50*V3A*V4A*V2D^3-190*V3A*V2A*V4A*V2D^2+270*V3A*V2A^2*V4A*V2D-162*V3A*V2A^3*V4A+...
        50*V3A^2*V1D*V2D^2+20*V3A^2*V2A*V1D*V2D+90*V3A^2*V2A^2*V1D-2*105*V1A*V3A^2*V2D^2+2*90*V1A*V3A^2*V2A*V2D-2*81*V1A*V3A^2*V2A^2,...
        -80*V2D+2*72*V2A+80*V4A*V1D+80*V3A*V2D-2*72*V3A*V2A-144*V1A*V4A+...
    50*V1D^2*V2D-2*105*V2A*V1D^2-160*V3A*V4A*V1D+80*V3A^2*V2D-2*72*V3A^2*V2A+...
    20*V1A*V1D*V2D+2*90*V1A*V2A*V1D+288*V1A*V3A*V4A+90*V1A^2*V2D-2*81*V1A^2*V2A-150*V4A*V1D*V2D^2+...
    2*190*V2A*V4A*V1D*V2D-3*90*V2A^2*V4A*V1D-100*V3A*V1D^2*V2D+2*210*V3A*V2A*V1D^2+80*V3A^2*V4A*V1D-80*V3A^3*V2D+2*72*V3A^3*V2A+...
    190*V1A*V4A*V2D^2-2*270*V1A*V2A*V4A*V2D+3*162*V1A*V2A^2*V4A-40*V1A*V3A*V1D*V2D-2*180*V1A*V3A*V2A*V1D-144*V1A*V3A^2*V4A-...
    180*V1A^2*V3A*V2D+2*162*V1A^2*V3A*V2A+100*V4A^2*V2D^3-2*190*V2A*V4A^2*V2D^2+3*180*V2A^2*V4A^2*V2D-4*81*V2A^3*V4A^2+150*V3A*V4A*V1D*V2D^2-2*190*V3A*V2A*V4A*V1D*V2D+...
    3*90*V3A*V2A^2*V4A*V1D+50*V3A^2*V1D^2*V2D-2*105*V3A^2*V2A*V1D^2-190*V1A*V3A*V4A*V2D^2+2*270*V1A*V3A*V2A*V4A*V2D-3*162*V1A*V3A*V2A^2*V4A+...
    20*V1A*V3A^2*V1D*V2D+2*90*V1A*V3A^2*V2A*V1D+90*V1A^2*V3A^2*V2D-2*81*V1A^2*V3A^2*V2A,...
    32-40*V2D^2-120*V1D^2-32*V4A^2+80*V2A*V2D-72*V2A^2-3*32*V3A^2+...
    240*V1A*V1D-216*V1A^2+160*V4A*V1D*V2D-160*V2A*V4A*V1D-2*40*V3A*V2D^2+2*120*V3A*V1D^2+2*16*V3A*V4A^2+2*80*V3A*V2A*V2D-2*72*V3A*V2A^2+4*16*V3A^3-...
    160*V1A*V4A*V2D+288*V1A*V2A*V4A-2*240*V1A*V3A*V1D+2*216*V1A^2*V3A+...
    50*V1D^2*V2D^2-100*V2A*V1D^2*V2D+210*V2A^2*V1D^2-2*80*V3A*V4A*V1D*V2D+2*80*V3A*V2A*V4A*V1D+3*40*V3A^2*V2D^2-3*40*V3A^2*V1D^2-3*80*V3A^2*V2A*V2D+3*72*V3A^2*V2A^2-...
    100*V1A*V1D*V2D^2-40*V1A*V2A*V1D*V2D-180*V1A*V2A^2*V1D+2*80*V1A*V3A*V4A*V2D-2*144*V1A*V3A*V2A*V4A+3*80*V1A*V3A^2*V1D+210*V1A^2*V2D^2-...
    180*V1A^2*V2A*V2D+162*V1A^2*V2A^2-3*72*V1A^2*V3A^2-50*V4A*V1D*V2D^3+150*V2A*V4A*V1D*V2D^2-190*V2A^2*V4A*V1D*V2D+...
    90*V2A^3*V4A*V1D-2*25*V3A*V1D^2*V2D^2+2*50*V3A*V2A*V1D^2*V2D-2*105*V3A*V2A^2*V1D^2+50*V1A*V4A*V2D^3-190*V1A*V2A*V4A*V2D^2+270*V1A*V2A^2*V4A*V2D-162*V1A*V2A^3*V4A+...
    2*50*V1A*V3A*V1D*V2D^2+2*20*V1A*V3A*V2A*V1D*V2D+2*90*V1A*V3A*V2A^2*V1D-2*105*V1A^2*V3A*V2D^2+2*90*V1A^2*V3A*V2A*V2D-2*81*V1A^2*V3A*V2A^2,...
    2*16*V4A-80*V1D*V2D+80*V2A*V1D-2*32*V3A*V4A+80*V1A*V2D-144*V1A*V2A+...
    160*V3A*V1D*V2D-160*V3A*V2A*V1D+2*16*V3A^2*V4A-...
    160*V1A*V3A*V2D+288*V1A*V3A*V2A+50*V1D*V2D^3-150*V2A*V1D*V2D^2+...
    190*V2A^2*V1D*V2D-90*V2A^3*V1D-80*V3A^2*V1D*V2D+80*V3A^2*V2A*V1D-50*V1A*V2D^3+...
    190*V1A*V2A*V2D^2-270*V1A*V2A^2*V2D+162*V1A*V2A^3+80*V1A*V3A^2*V2D-144*V1A*V3A^2*V2A-...
    2*25*V4A*V2D^4+2*100*V2A*V4A*V2D^3-2*190*V2A^2*V4A*V2D^2+2*180*V2A^3*V4A*V2D-2*81*V2A^4*V4A-50*V3A*V1D*V2D^3+150*V3A*V2A*V1D*V2D^2-190*V3A*V2A^2*V1D*V2D+...
    90*V3A*V2A^3*V1D+50*V1A*V3A*V2D^3-190*V1A*V3A*V2A*V2D^2+270*V1A*V3A*V2A^2*V2D-162*V1A*V3A*V2A^3,...
    2*40*V1D-80*V1A-80*V4A*V2D+80*V2A*V4A-2*120*V3A*V1D+...
    240*V1A*V3A-2*25*V1D*V2D^2+2*50*V2A*V1D*V2D-2*105*V2A^2*V1D+160*V3A*V4A*V2D-160*V3A*V2A*V4A+2*120*V3A^2*V1D+...
    50*V1A*V2D^2+20*V1A*V2A*V2D+90*V1A*V2A^2-240*V1A*V3A^2+50*V4A*V2D^3-150*V2A*V4A*V2D^2+...
    190*V2A^2*V4A*V2D-90*V2A^3*V4A+2*50*V3A*V1D*V2D^2-2*100*V3A*V2A*V1D*V2D+2*210*V3A*V2A^2*V1D-80*V3A^2*V4A*V2D+80*V3A^2*V2A*V4A-2*40*V3A^3*V1D-...
    100*V1A*V3A*V2D^2-40*V1A*V3A*V2A*V2D-180*V1A*V3A*V2A^2+80*V1A*V3A^3-...
    50*V3A*V4A*V2D^3+150*V3A*V2A*V4A*V2D^2-190*V3A*V2A^2*V4A*V2D+...
    90*V3A*V2A^3*V4A-2*25*V3A^2*V1D*V2D^2+2*50*V3A^2*V2A*V1D*V2D-2*105*V3A^2*V2A^2*V1D+...
    50*V1A*V3A^2*V2D^2+20*V1A*V3A^2*V2A*V2D+90*V1A*V3A^2*V2A^2,...
    2*40*V2D-80*V2A-80*V4A*V1D-2*40*V3A*V2D+80*V3A*V2A+80*V1A*V4A-...
    2*25*V1D^2*V2D+50*V2A*V1D^2+160*V3A*V4A*V1D-2*40*V3A^2*V2D+80*V3A^2*V2A+...
    2*50*V1A*V1D*V2D+20*V1A*V2A*V1D-160*V1A*V3A*V4A-2*105*V1A^2*V2D+90*V1A^2*V2A+3*50*V4A*V1D*V2D^2-2*150*V2A*V4A*V1D*V2D+...
    190*V2A^2*V4A*V1D+2*50*V3A*V1D^2*V2D-100*V3A*V2A*V1D^2-80*V3A^2*V4A*V1D+2*40*V3A^3*V2D-80*V3A^3*V2A-3*50*V1A*V4A*V2D^2+...
    2*190*V1A*V2A*V4A*V2D-270*V1A*V2A^2*V4A-2*100*V1A*V3A*V1D*V2D-40*V1A*V3A*V2A*V1D+80*V1A*V3A^2*V4A+2*210*V1A^2*V3A*V2D-...
    180*V1A^2*V3A*V2A-4*25*V4A^2*V2D^3+3*100*V2A*V4A^2*V2D^2-2*190*V2A^2*V4A^2*V2D+180*V2A^3*V4A^2-3*50*V3A*V4A*V1D*V2D^2+2*150*V3A*V2A*V4A*V1D*V2D-190*V3A*V2A^2*V4A*V1D-...
    2*25*V3A^2*V1D^2*V2D+50*V3A^2*V2A*V1D^2+3*50*V1A*V3A*V4A*V2D^2-2*190*V1A*V3A*V2A*V4A*V2D+270*V1A*V3A*V2A^2*V4A+...
    2*50*V1A*V3A^2*V1D*V2D+20*V1A*V3A^2*V2A*V1D-2*105*V1A^2*V3A^2*V2D+90*V1A^2*V3A^2*V2A,...
    0,0,0];
end
%