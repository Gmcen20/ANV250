function distance=HausdorffDist(A,B)
[n1,m1]=size(A);
[n2,m2]=size(B);
min1=zeros(n1,1);
for i=1:n1
    min=norm(A(i,:)-B(1,:));
    for j=2:n2
        if norm(A(i,:)-B(j,:))<min
            min=norm(A(i,:)-B(j,:));
        end
    end
    min1(i)=min;
end 
a=max(min1);
min2=zeros(n2,1);
for i=1:n2
    min=norm(A(1,:)-B(i,:));
    for j=2:n1
        if norm(A(j,:)-B(i,:))<min
            min=norm(A(j,:)-B(i,:));
        end
    end
    min2(i)=min;
end 
b=max(min2);
if a<b
    distance=b;
else
    distance=a;
end
end