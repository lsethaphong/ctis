function Kterm = sps_Kcon(Upinv,Vm,Cpinv,Cm);
% constructs P matrix from Upinv
% 
prod = sps_matmul(fnameX,'x',fnameH,'h',fnameI,1,1,1,itol);

return