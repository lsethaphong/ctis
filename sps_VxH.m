function vtc = sps_VxH(comV,fnameX,scalef);% return compressed vector
% Latsavongsakda Sethaphong
% H is condensed to signed col, row successive format
%fid = fopen(fnameX); % for the maximum basis set with the modeled code
%[a] = fscanf(fid,'%e',[inf]);
%fclose(fid);
%basisX = a'; %
%clear a;
x = load(fnameX,'x');
vtc = zeros(1,2);
vtc(1,1) = 1;
vtc(1,2) = comV(1,1); % number of columns %% taken of the transpose
vt=zeros(4096*2048,1); % number of columns
% skip the first row
tol = 5.0e-7;
% fill temp vector with nonzero elements of comV
for i =2:length(comV)
    vt(abs(comV(i,2))) = comV(i,1);
end
%
summand = 0;
idx = 1;
colid = 0;
for j = 2:length(basisX) % n by m
  summand = summand + x(j,1)*vt(abs(x(j,2)));
  if basisX(j,2) < 0 % if the signed column is  negative
      % take the sum and increment the col id
      summand = summand/scalef;
      colid = colid + 1;
      if abs(summand) > tol
          idx = idx+1;
         vtc(idx,1)=summand;
         vtc(idx,2)=colid;
      end
      summand =0; % clear summand
  end
end
vtc(idx,2)=-vtc(idx,2);
return;