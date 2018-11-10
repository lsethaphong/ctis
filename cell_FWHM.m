%========================================
% This function takes in filtered data from the cell_analyze function
% and produces FWHM information for each cell in the cohort
%========================================
function cell_FWHM(srcfname,idx_start,idx_stop,ncolumns,outfname);
% for John Rocheleau
% Preliminary code
% has all length as a dyad 2^n 
% in this case 128

  mdirf ='C:/matlab_sv13/work/JR/';
%  desf = 'fwhm';
  desf = outfname;
  dfwhm = zeros(1,20);
  for j = idx_start:idx_stop
      srcf= [srcfname int2str(j)];
      fid = fopen([mdirf srcf '.dat']); % for the maximum basis set with the modeled code
      [e] = fscanf(fid,'%e',[ncolumns inf]); % It has two rows now.
      fclose(fid);
      E = e';
      dimE =size(E);
   for i = 2:dimE(2)
      tmph = 0;
      tmph_idx = 0;
%      tmpl = 0;
      tmpl_idx = 0;
      dif = 0;
      ldif = 0;
      memh = [1 0 0]; % choose blank
      meml = [1 0 0]; % choose blank
      m_results = [0 0];
      avg_i = mean(E(:,i)); % use mean as a threshold
      hval = max(E(:,i)); % get the highest value
      lval = min(E(:,i)); % get the lowest value
      tval = (hval - lval)/2 + lval; % threshold
      tmpl = avg_i;
      for k = 2:128
          % E(k,i) % search through a column of data
             dif = E(k,i)-E(k-1,i); % derivative          
          % turning points
          if dif < 0 && ldif > 0 && E(k,i) > avg_i % looking for downtick
              % get highest point
              if tmpl_idx > 0
                 meml = [meml ; tmpl_idx E(tmpl_idx,1) tmpl]; % store index, time and value
                 tmpl = avg_i; % reset to average
                 tmpl_idx = 0;
              end
              % get highest point
              if E(k,i) > tmph
                  tmph = E(k,i); % record value
                  tmph_idx = k; % record index 
              end
          elseif dif > 0 && ldif < 0 && E(k,i) < avg_i % looking for uptick
              % get lowest point
              if tmph_idx > 0 
                  memh = [memh ; tmph_idx E(tmph_idx,1) tmph]; % store index, time and value
                  tmph = 0;
                  tmph_idx = 0;
              end
              % get lowest point
              if E(k,i) < tmpl
                  tmpl=E(k,i); % record value
                  tmpl_idx = k; % record index
              end
          end
             ldif=dif;
          % use as trigger to get lowest point
      end % end of for k
      % check for non matched points
      if tmph_idx > 0
         memh = [memh ; tmph_idx E(tmph_idx,1) tmph]; % store index, time and value
         tmph = 0;
         tmph_idx = 0;          
      end
      if tmpl_idx > 0
         meml = [meml ; tmpl_idx E(tmpl_idx,1) tmpl]; % store index, time and value
         tmpl = avg_i; % reset to average
         tmpl_idx = 0;                    
      end
      fwhm_i = [0];
      for t=2:length(memh)
          if length(memh)>length(meml) && t > length(meml)
             half_v = (memh(t,3) - meml(t-1,3))/2 + meml(t-1,3);
          else
             half_v = (memh(t,3) - meml(t,3))/2 + meml(t,3);
          end% half value at this point
          if length(memh) > length(meml)
              % typically differs by one
             if t <= length(meml)
                 ist = meml(t-1,1);
                 isp = meml(t,1);
             else
                 ist = meml(t-1,1);
                 isp = length(E(:,1));
             end
         elseif length(meml) > length(memh) % started with a low
             ist = meml(t,1);
             isp = meml(t+1,1);
         else
            ist = meml(t-1,1);
            isp = meml(t,1);             
         end
                str = 0;
                for u=ist+1:isp
                    % perform linear interpolation
                   if E(u,i) >= half_v && str == 0
                       dx_i = E(u,i)-E(u-1,i);
                       dt_i = E(u,1)-E(u-1,1);
                       mt_i = dt_i/dx_i;
                       ts_i = (E(u,i)-half_v)*(dt_i/dx_i) + E(u,1);
                       str = 1;
                   elseif E(u,i) <= half_v && str == 1
                       dx_i = E(u-1,i)-E(u,i);
                       dt_i = E(u-1,1)-E(u,1);
                       mt_i = dt_i/dx_i;
                       tr_i = (half_v-E(u,i))*(dt_i/dx_i) + E(u,1);             
                       str = 0;
                   end
                end   
                % store first fwhm
                fwhm_i = [fwhm_i;tr_i - ts_i];
      end
      % compile fwhm
      fwhm_i(1)=mean(fwhm_i(2:length(fwhm_i))); 
      dummy  =zeros(1,20-length(fwhm_i));
      fwhm_i = cat(1,fwhm_i, dummy');
      dfwhm = [dfwhm;fwhm_i'];
   end
      csvwrite([mdirf desf '.csv'],dfwhm);
  end
%  save();
return;