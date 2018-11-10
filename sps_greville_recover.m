function sps_greville_recover;
% from Jiri Rohn
% Institute of Computer Science
% Academy of Sciences of the Czech Republic
% Prehled nekterych dulezitych vet z teorie matic
% 5 Sep 2003
% Technical Report No. 893
% -------------------------
% modified for spares matrix 
% The inverse will be saved by rowId to expedite
% multiplication
% -------------------------
% Utilizing Greville with Sparse matrices and vectors
% implementation of Greville
% 12 OCT 04
% Latsavongsakda Sethaphong
% -------------------------
% the system matrix is row sparse, hence, it is indexed by column
% hence, the inverse is column sparse, and, it is indexed by row.
% -------------------------
fpa_x = 4096;
fpa_y = 2048;
src_x = 433;
src_y = 433;
nz_pts =127;
vxl_depth=41;
n = src_x*src_y*vxl_depth; % col id
% --------------------------
% --------------------------
% colId goes from 1 to src_x*src_y*vxl_depth
% rowId goes from 1 to fpa_x*fpa_y
% spsMat(colId,rowId,fpa_x,fpa_y,src_x,src_y,vxl_depth,nz_pts); % column format
% prod = sps_matmul(fnameA,fnameB,fnameC,func,b_ki,b_kn) % multiplication
% the 
%
% to_sgcol(c,fname) conversion to signed column row successive format
% to_sgrow(c,fname) conversion to signed row, column successive format
%[m,n]=size(A);
fnameA = 'c:/matlab_sv13/work/acol.mat';
fnameB = 'c:/matlab_sv13/work/bcol.mat';
fnameC = 'c:/matlab_sv13/work/ccol.mat';
fnameD = 'c:/matlab_sv13/work/dcol.mat';
fnameX = 'c:/matlab_sv13/work/xcol.mat';
fnameK = 'c:/matlab_sv13/work/kcol.mat';
try
tol=5.0e-7;
%d=A(:,1); % col vector
% X is row vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%d = spsMat(1,-3,4096,2048,433,433,41,127);
%if all(abs(d)<tol*ones(fpa_x*fpa_y,1))
%    x=zeros(1,fpa_x*fpa_y); 
%else
%    x=d'/(d'*d);
%end
%    to_sgcol(x,fnameX); % save variable X
for j=2:n %n% j is also pixel id
    % store X by rows and column
    % compute
    %d = X*A(:,j); % full vector expansion
    d = sps_HxV(fnameX,spsMat(j,-3,4096,2048,433,433,41,127)); % return full vector d
    % compute
%    return;
    % c= A(:,j)-A(:,1:(j-1))*d;
    c = spsMat(j,-3,4096,2048,433,433,41,127) - sps_MxV(1,(j-1),d); % return full vector c
    if all(abs(c)<tol*ones(fpa_x*fpa_y,1))
        kt=sps_VxH(d,fnameX)/(1+d'*d); % return full vector kt
        msg = 'not in tolerance'
    else 
        kt=c'/(c'*c); % return full vector kt
        clear c;
%        msg = 'in tolerance'
    end
    % use compressed X as well
    % check out nonzero d elements
    % take and save difference of Xprev and new column id in sparse format
    % signed column, row successive
    % convert new kt to signed column, row successive format and append to 
    % fnameX;   
    to_sgrow(d,fnameD); % turn d into signed row, col successive
    to_sgcol(kt,fnameK); % turn kt into signed col, row successive
    sps_VxVt(fnameD,fnameK,fnameC); % output signed column, row successive
    sps_sum(fnameX,fnameC,fnameX,-1); % return signed column, row-successive format
    %X=[X-d*kt; kt]; % index for the difference Write to a file and tag on the last row for good measure
    %tmpkt = kt(2:length(kt));
    X_append(kt,fnameX,-1);
    %save(fnameX, tmpkt, '-append -ASCII');
    % append signed column row successive format to new X file;
end
msg = '------------- finis ------------'
return;
%save('c:/matlab_sv13/work/fullspsH.dat','X','-ASCII');
catch
msg = '--------- memory error ---------'
j
return;
end