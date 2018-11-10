function sps_greville(n_hcol,tol,ftol,fpa_x,fpa_y,src_x,src_y,nz_pts,vxl_depth); %input tolerance level
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
% 19 OCT 04
% Latsavongsakda Sethaphong
% -------------------------
% the system matrix is row sparse, hence, it is indexed by column
% hence, the inverse is column sparse, and, it is indexed by row.
% -------------------------
%fpa_x = 4096;
%fpa_y = 2048;
%src_x = 433;
%src_y = 433;
%nz_pts =127;
%vxl_depth=41;
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
fnameA = [dirname 'av' '.mat'];
fnameB = [dirname 'bv' '.mat'];
fnameC = [dirname 'cv' '.mat'];
fnameD = [dirname 'dv' '.mat'];
fnameX = [dirname 'xinv' '.mat']; % inverse of basis projection at pixel 1
fnameK = [dirname 'kv' '.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% index files
fnameAi = [dirname 'ai' '.mat'];
fnameBi = [dirname 'bi' '.mat'];
fnameCi = [dirname 'ci' '.mat'];
fnameDi = [dirname 'di' '.mat'];
fnameXi = [dirname 'xinvi' '.mat']; % inverse of basis projection at pixel 1
fnameKi = [dirname 'ki' '.mat'];

%try
%tol=5.0e-6;
%d=A(:,1); % col vector
% X is row vector
x=ones(1,2);
d = spsMat(1,-1,4096,2048,433,433,41,127);
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
for j=2:n_hcol %vxl_depth %n% j is also pixel id
    t=toc;
    msg = ['=== started column ' int2str(j) ' | ' int2str(t) ' secs ===']
    % store X by rows and column
    % compute
    %d = X*A(:,j); % full vector expansion
    % input current invx and column j of x
    msg = ['--- making d ---']
    d=sps_HxV(fnameX,spsMat(j,-1,4096,2048,433,433,41,127),fnameD,tol);
    % c= A(:,j)-A(:,1:(j-1))*d;
    % signed row, col successive difference
    % return  as column vector compressed rows
    msg = ['--- making c ---'] 
    c = sps_A_MxV(1,(j-1),d,tol);
    if all(abs(c(2:length(c),1))<tol*ones(length(c)-1,1))
        %     if all(abs(c)<tol*ones(m,1)), bt=d'*X/(1+d'*d); else
        %     bt=c'/(c'*c); end
        % return compressed vector kt, row vector, compressed columns
        kt=sps_VxH(d,fnameX,(1+d(2:length(d),1)'*d(2:length(d),1))); 
        msg = 'not in tolerance'
        save(fnameK,'kt','-mat');
    else 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kt=zeros(length(c),2); % return row vector, compressed rows
        kt(2:length(c),1)=c(2:length(c),1)/(c(2:length(c),1)'*c(2:length(c),1)); 
        kt(2:length(c),2)=c(2:length(c),2); % keep the row id's which have become col id's
        kt(1,1)=c(1,2);
        kt(1,2)=c(1,1);
        save(fnameK,'kt','-mat');
        % msg = 'in tolerance'
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
        % return signed column, row-successive format
       case vxl_depth
          sps_sum(fnameX,fnameC,fnameX,-1,ftol); % final tolerance (may be a bad idea b/c kt)
       otherwise
          sps_sum(fnameX,fnameC,fnameX,-1,tol); % main tolerance
    end
    % X=[X-d*kt; kt];
    % index for the difference Write to a file and tag on the last row for good measure
    X_append(fnameX,kt(2:length(kt),:)); % snippet off the header for kt
    t=toc; % measure elapse time
    msg = ['========= completed column ' int2str(j) ' | ' int2str(t) ' secs =========']
end
msg = '***** CALCULATED BASIS H PSEUDOINVERSE *****'
%%%%%%%%%%%%%%%%%%%%%%%%%
% getfname(j,'c:/matlab_sv13/work/xinv');
return;
%catch
%msg = '--------- memory error ---------'
%j
%return;
%end