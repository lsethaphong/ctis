function [X] = greville(A);
% from Jiri Rohn
% Institute of Computer Science
% Academy of Sciences of the Czech Republic
% Prehled nekterych dulezitych vet z teorie matic
% 5 Sep 2003
% Technical Report No. 893
[m,n]=size(A); tol=1.0e-10;
d=A(:,1);
if all(abs(d)<tol*ones(m,1)), X=zeros(1,m); else X=d'/(d'*d); end
for j=2:n
    d=X*A(:,j);
    c=A(:,j)-A(:,1:(j-1))*d;
    if all(abs(c)<tol*ones(m,1)), bt=d'*X/(1+d'*d); else bt=c'/(c'*c); end
    X=[X-d*bt; bt]; % Matlab automatically adds a row
end