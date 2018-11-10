function valid_I = chk_Ident(src_size,tgt_sizeX,tgt_sizeY,vxl_depth); 
Ho=construct_basis(src_size,tgt_sizeX,tgt_sizeY,vxl_depth);
size(Ho)
%Hp=pinv(Ho);
Hp = greville(Ho);
Hi = Hp*Ho;
cnt1 = 0.0;
cnt2(src_size*src_size*vxl_depth) = 0.0;
for jj = 1:tgt_sizeX*tgt_sizeY % row
    for ii = 1:src_size*src_size*vxl_depth % column
            cnt1 = cnt1 + 1;
            cnt2(ii) = cnt2(ii) + Ho(jj,ii);
    end
end
%Ho(:,1)
%Hp(1:7,1:7)
clear Ho;
clear Hp;
cnt = 0.0;
Hi(1:7,1:7)
size(Hi)
cnt2
for ii = 1:src_size*src_size*vxl_depth
        cnt = cnt + Hi(ii,ii); % check the diagonals 
end
if cnt < double(src_size*src_size*vxl_depth)
    cnt
    message = 'failure'
else
    message = 'success'
end
clear Hi;
valid_I=cnt
