function m = recombineJR;
for i = 2:7
fid = fopen(['c:/matlab_sv13/work/JR/trc' int2str(i) '.dat']); % for the maximum basis set with the modeled code
[d] = fscanf(fid,'%e',[1 128]); % It has two rows now.
fclose(fid);
   if i == 2
      m= d;
   else
      m = [m ; d];
   end
end
fid = fopen('JR/ct.dat'); % for the maximum basis set with the modeled code
[a] = fscanf(fid,'%e',[7 150]); % It has two rows now.
fclose(fid);
A=a';
m = [A(1:128,1)';m]';
save('c:/matlab_sv13/work/JR/treat_cells.dat','m','-ASCII');
return;