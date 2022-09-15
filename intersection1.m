function T=intersection1(A,B);
[m,l]=size(A);
[n,l]=size(B);
c=zeros(m+n,1);
A0=[A',-B'];
a1=[ones(1,m),zeros(1,n)];
b1=[zeros(1,m),ones(1,n)];
Aeq=[A0;a1;b1];
beq=[zeros(l,1);1;1];
lb=zeros(m+n,1);
ub=ones(m+n,1);
s=linprog(c,[],[],Aeq,beq,lb,ub);
if length(s)==0
    T=1;
else 
    T=0;
end
end