!======================================================================
! This program returns a compressed mapping matrix of dimensions
! spectral section by number of nonzero points: 41-57x127 for 
! third order detection which is not likely to be attainable.
!
!======================================================================
       subroutine SPSH(colId,rowId,fpa_x,fpa_y,src_x,src_y,vxl_depth,nz_pts)
!% 15 September 04
!% Latsavongsakda Sethaphong
!% Shifting of Basis Projection for expansion into full transfer matrix 
!% To make life easier, the source dimensions should be assumed to be odd
!% read in basisH file
        fid = fopen('matdata.dat');
!% for the maximum basis set with the modeled code
        [a] = fscanf(fid,'%e',[vxl_depth nz_pts*2]); !% It has two rows now.
        fclose(fid);
        basisH = transpose(a); % 254 rows, 41 columns = voxel depth
        clear a;
!% this function takes the shift invariant 
!% basis projection at the center of the FPA
!% and shifts the results to the correct 
!% representation for an fpa_x by fpa_y target
!% from a unit-source of a src_x by src_y dimensions
!% this is a column by column operation that shifts based on
!% groups and pixels
!% cg2voxel scans along the x axis

!% locate nonzero points per column of the basis projection
!% find nonzero components of basisM and record rowId per column
!%nsz=(343*2,vxl_depth); % this is what the basis M looks like
!% the basis matrix H is 343*2 by vxl_depth
!% two pointers are used to dereference the column data
!% pt & Id_pt, one for the value, the other for the column Id
!% calculated desired position and translate
!% calculates nonzero components row by row of the transfer
!% matrix H
!% stores it by row files with contents of column id
!% find basis column

!% assume that the source distance is odd
        ref_x = (src_x-1)/2 + 1;
        ref_y = (src_y-1)/2 + 1;
!% the colId gives us the voxel of interest
        vxl_basis_id = mod(colId,vxl_depth);

        if vxl_basis_id == 0
            vxl_basis_id = vxl_depth;
        end
!% determine source displacement
        pixel_x_id = mod((colId-mod(colId,vxl_depth)),vxl_depth*src_x)/vxl_depth + 1; % get the pixel 
        pixel_y_id = (colId-mod(colId,vxl_depth*src_x))/(vxl_depth*src_x) + 1; % take the whole number
!% id's of zero means it's at the bottom right corner or at the last column
!% displacement vectors
        shift_x = pixel_x_id-ref_x;
        shift_y = pixel_y_id-ref_y;
        shift_f = shift_x + shift_y*fpa_x; !% shifting along the pixels of the FPA which are row elements of the
!% transform matrix H
!% set shift conditions for the rowId's
        t_ptr = 0;
        for ptr = 1:nz_pts
           newrow = basisH(ptr+sz_H(1)/2,vxl_basis_id) + shift_f;
           if rowId == newrow
               t_ptr = ptr;
               ptr = nz_pts; % stop the for loop
           endif
        end
           if t_ptr > 0
               valM=basisH(t_ptr,vxl_basis_id);
           else
               valM = 0;
           endif
        return;
!% calculate nonzero row elements based on position
      end
      end
