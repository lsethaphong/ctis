function sps_greville_rowexp(i_hcol,n_hcol,tol,ftol,func); %input tolerance level
% from Jiri Rohn
% Institute of Computer Science
% Academy of Sciences of the Czech Republic
% Prehled nekterych dulezitych vet z teorie matic
% 5 Sep 2003
% Technical Report No. 893
% -------------------------
% modified for sparse matrix 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first calculate hbasis inverse
% then compute the corresponding shifted matrix in sparse format
% then use the partition theorem to grow the inverse matrix by 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[m,n]=size(A);
dirname = 'c:/matlab_sv13/work/xdata/';
fnameA = [dirname 'ar' '.mat'];
fnameB = [dirname 'br' '.mat'];
fnameC = [dirname 'cr' '.mat'];
fnameD = [dirname 'dr' '.mat'];
fnameX = [dirname 'xinv' '.mat']; % inverse of basis projection at pixel 1
fnameK = [dirname 'kr' '.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% index files
fnameAi = [dirname 'ai' '.mat'];
fnameBi = [dirname 'bi' '.mat'];
fnameCi = [dirname 'ci' '.mat'];
fnameDi = [dirname 'di' '.mat'];
fnameXi = [dirname 'xinvi' '.mat']; % inverse of basis projection at pixel 1
fnameKi = [dirname 'ki' '.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% row expansion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% row shifting

%try
%tol=5.0e-6;
%d=A(:,1); % col vector
% X is row vector
if func == 1 % basis matrix calculation
   x=ones(1,2);
   d = spsMat(i_hcol,-1,4096,2048,433,433,41,127);
   x(1,1)=d(1,2); 
   x(1,2)=d(1,1);
if all(abs(d(2:length(d),1)))<tol*ones(fpa_x*fpa_y,1)
    % compressed form
    x(2,1)=0;
    x(2,2)=-1;
else
    tpv =d(2:length(d),1)'*d(2:length(d),1);
    for j= 2:length(d)
       x(j,1)=d(j,1)/tpv; % get the scaling
       x(j,2)=d(j,2);
    end
end
    save(fnameX,'x','-mat');
    
    tic; % start time
for j=(1+i_hcol):n_hcol %vxl_depth %n% j is also pixel id
    t=toc;
    msg = ['=== started column ' int2str(j) ' | ' int2str(t) ' secs ===']
    % store X by rows and column
    % compute
    %d = X*A(:,j); % full vector expansion
    % input current invx and column j of x
    msg = ['--- making d ---']
    d=sps_HxV(fnameX,spsMat(j,-1,4096,2048,433,433,41,127),fnameD,tol); % get a new d vector
    % c= A(:,j)-A(:,1:(j-1))*d;
    % signed row, col successive difference
    %return  as column vector compressed rows
    msg = ['--- making c ---'] 
    c = sps_A_MxV(i_hcol,(j-1),d,tol); % previously i_hcol would be the first column
    % 
    if all(abs(c(2:length(c),1))<tol*ones(length(c)-1,1))
        %     if all(abs(c)<tol*ones(m,1)), bt=d'*X/(1+d'*d); else
        %     bt=c'/(c'*c); end
        kt=sps_VxH(d,fnameX,(1+d(2:length(d),1)'*d(2:length(d),1))); % return compressed vector kt, row vector, compressed columns
        msg = 'not in tolerance'
        save(fnameK,'kt','-mat');
    else 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kt=zeros(length(c),2);
        kt(2:length(c),1)=c(2:length(c),1)/(c(2:length(c),1)'*c(2:length(c),1)); % return row vector, compressed rows
        kt(2:length(c),2)=c(2:length(c),2); % keep the row id's which have become col id's
        kt(1,1)=c(1,2);
        kt(1,2)=c(1,1);
        save(fnameK,'kt','-mat');
%        msg = 'in tolerance'
    end
    % use compressed X as well
    % check out nonzero d elements
    % take and save difference of Xprev and new column id in sparse format
    % signed column, row successive
    % convert new kt to signed column, row successive format and append to 
    % fnameX;
    
    sps_VxVt(fnameD,fnameK,fnameC,tol); % output signed column, row successive
    msg = '---- made subtractand ----'
    switch j % final tolerance setting
       case vxl_depth
          sps_sum(fnameX,fnameC,fnameX,-1,ftol); % return signed column, row-successive format
       otherwise
          sps_sum(fnameX,fnameC,fnameX,-1,tol); % return signed column, row-successive format
    end
    %X=[X-d*kt; kt]; % index for the difference Write to a file and tag on the last row for good measure
    X_append(fnameX,kt(2:length(kt),:)); % snippet off the header for kt
    t=toc; % measure elapse time
    msg = ['========= completed column ' int2str(j) ' | ' int2str(t) ' secs =========']
end
msg = '***** CALCULATED BASIS H PSEUDOINVERSE *****'
%%%%%%%%%%%%%%%%%%%%%%%%%
elseif func == 2 % full row expansion
% get initial subpartition inverse name    
fnameU = [dirname 'ur' '.mat'];
fnameR = [dirname 'rr' '.mat'];
fnameJ = [dirname 'jr' '.mat'];
fnameP = [dirname 'pr' '.mat'];
fnameV = [dirname 'vr' '.mat'];
fnameUp = [dirname 'uinv' '.mat'];
fnameWp = [dirname 'winv' '.mat']
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Row Basis inversion & storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fnameU,'h'); % starting partition
% (U:V) = (U+. - U+. V J ; J) = (U+. - EJ ; J)
% E = U+.V
% C = (I-UU+.)V = V - UE
% P = [U+. V(I - C+. C] = E - EC+.C
% K = inv{I + tran(P)*P}
% J = C+. + (I-C+.C)K tran(V)*tran(U+.)*U+.*(I-VC+.)
% dim(C+.) = dim(Vt.V);
for k=2:(src_x-1)/2 % right shifting
    % calculate K
    
end

elseif func == 3 % full down expansion
% get initial row inverse name    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full multi-row shift and inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=2:(src_y-1)/2 % down shifting
    % calculate C
    % calculate K
    % calculate J

end

else
    msg = '**** UNKNOWN FUNCTION ****'
end
% getfname(j,'c:/matlab_sv13/work/xinv');
return;
%catch
%msg = '--------- memory error ---------'
%j
%return;
%end