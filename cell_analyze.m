%==================================
%  this program opens up raw data from the local directory
%  numbers them as an output
%  and outputs them to a new data set of separate files
%==================================
function cell_analyze;
% NON TREATED CELLS
fid = fopen('JR/cn.dat'); % for the maximum basis set with the modeled code
[a] = fscanf(fid,'%e',[6 150]); % It has two rows now.
fclose(fid);
A=a';
fid = fopen('JR/ncn.dat'); % for the maximum basis set with the modeled code
[b] = fscanf(fid,'%e',[6 150]); % It has two rows now.
fclose(fid);
B=b';
% FGF TREATED CELLS
fid = fopen('JR/ct.dat'); % for the maximum basis set with the modeled code
[c] = fscanf(fid,'%e',[7 150]); % It has two rows now.
fclose(fid);
C=c';
fid = fopen('JR/nct.dat'); % for the maximum basis set with the modeled code
[d] = fscanf(fid,'%e',[7 150]); % It has two rows now.
fclose(fid);
D=d';

%========== Denoising ===========
mdirf = 'c:/matlab_sv13/work/JR/';
for j=1:4
    switch j
        case 1
            tmpM = A; %unnormalised untreated
            srcf = 'raw_u_cell';
            outf = 'f_d6_raw_u_cell';
        case 2
            tmpM = B; %normalised untreated
            srcf = 'norm_u_cell';
            outf = 'f_d6_norm_u_cell';
        case 3
            tmpM = C; %unormalised FGF treated
            srcf = 'raw_t_cell';
            outf = 'f_d6_raw_t_cell';
        case 4 
            tmpM = D; %normalised FGF treated
            srcf = 'norm_t_cell';
            outf = 'f_d6_norm_t_cell';
        otherwise
            return;
    end
    
   if j < 3, stp = 6; else stp = 7; end; 
    
   for i=2:stp
   %rawsignal = D(20:147,i)'; % must be dyadic 2^number
   if j < 3
      rawsignal = tmpM(1:128,i)'; % must be dyadic 2^number
   else
      rawsignal = tmpM(20:147,i)'; % must be dyadic 2^number
   end
   QMF8 = MakeONFilter('Daubechies',6); %D8 filter coefficients low pass
   %QMF8 = MakeONFilter('Coiflet',2); %Coifman filter coefficients low pass
   %QMF8 = MakeONFilter('Symmlet',8); %symmetric filter filter coefficients low pass
   scaleraw = NormNoise(rawsignal,QMF8);
   y        = scaleraw*1.0; % make it a double
   [xh,wcoef]= WaveShrink(y,'Visu',4,QMF8);
   tsig = 1:length(rawsignal);
%
   clf;
   versaplot(211,tsig,y, [],' 1 (a) Islet Cell Flourscence',[],[])
   versaplot(212,tsig,xh,[],' 1 (b) Wavelet Shrinkage De-Noising',[],[])
   save([mdirf srcf int2str(i) '.dat'],'xh','-ASCII'); % saving the cleaned signal
   end

   for i = 2:stp
   fid = fopen([mdirf srcf int2str(i) '.dat']); % for the maximum basis set with the modeled code
   [e] = fscanf(fid,'%e',[1 128]); % It has two rows now.
   fclose(fid);
      if i == 2
         m= e;
      else
         m = [m ; e];
      end
   end
   if j < 3
      m = [tmpM(1:128,1)';m]';
   else 
      m = [tmpM(20:147,1)';m]';
   end
   save([mdirf outf '.dat'],'m','-ASCII');
end
% ================================
%PlotWaveCoeff(wc_i,32,0);
return;