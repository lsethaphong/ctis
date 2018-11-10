function H_inverse = gaussj(H_matrix);
% inversion 
%
%
% 1024x1024x3 matrix
dimension = size(H_matrix);
nn = dimension(1); % number of rows
mm = dimension(2); % number of columns
b_matrix = zeros
irow=1;
icol=1;
indxc=zeros(nn);
indxr=zeros(nn);
ipiv=zeros(nn);
%    ipiv(j)=0; %already initialized
for i=1:nn
    big=double(0.0);
    for j=1:nn
        if ipiv(j) ~= 1
            for k=1:nn
                if ipiv(k) == 0
                    if abs(H_matrix(j,k)) >= big
                        big=abs(H_matrix(j,k));
                        irow=j;
                        icol=k;
                    end
                end
            end
            ipiv[icol]=ipiv[icol]+1;
            if irow ~= icol
                for l =1:nn
                end
                for k =1:mm
                end
            end
        end
    end
end
