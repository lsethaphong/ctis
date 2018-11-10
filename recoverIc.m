function recoverIc(FPA_sizeX,FPA_sizeY,vxl_depth,obj_size);
  Ho = construct_basis(obj_size,FPA_sizeX,FPA_sizeY,vxl_depth);
%  Hp = pinv(Ho);
  Hp = greville(Ho);
  nnz(Ho)
  nnz(Hp)
  Hos = sparse(Ho);
  Hps = sparse(Hp);
  save('c:\matlab_sv13\work\Hos.mat','Hos','-ASCII');
  save('c:\matlab_sv13\work\Hps.mat','Hps','-ASCII');
  spy(Hps)
  return;