%============================================================
% Auto generate putative inverse matrix based on shift invariance
% on pixel data based on the Sherwood Morrison Woodbury and
% T. Greville's method for partitioned matricies
% 1 NOV 04
% Latsavongsakda Sethaphong
% VUMC PISTON LAB
%============================================================
function hpinv = sps_hpinv(pxl_id,vxl_id,rtype);
% vxl_id selects the sub ro
load('hpinv_1en7.mat','x');
% x is hpinv/ U+ with size (186545,2)
% general format is such
% | U+.(At^n).(Dt^m) | for n = 0 to 432, and m = 0 to 432

switch rtype
    case 1 % single compressed row
    case 2 % single 
    case 3 % base inverse return
    otherwise
        msg = 'unknown return type'
end
return;