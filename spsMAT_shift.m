function smat = spsMAT_shift(imat,idx_i,rshift,dshift,fpa_x); 
% assume the input sparse matrix is in sparse row, column successive format
% transpose of shift is equivalent to inverse because the right
% and down operator are unitary
mdshift=dshift*fpa_x + rshift; % net displacement
for i=idx_i:length(imat)
    if imat(i,2) < 0
       imat(i,2)= -(abs(imat(i,2)) + mdshift);
    else
       imat(i,2)= imat(i,2) + mdshift;
    end
end
smat = imat;
return;