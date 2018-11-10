function C = borinfunc;
A=load('d:\2rr1col.txt');
reshape(A,1024,1);
for m=1:1024
    if A(m)<0
        A(m)=0
    end
end

for m=7:1018;
    B(m)=-A(m-6)+A(m-2)+A(m+2)-A(m+6);
    if B(m)<0
        B(m)=0;
    end
end
% no semicolon and we return the variable to something more permanent 
C = B
return;