function flag_family=knn(x,AA,index,n1,k)
    [n,m]=size(AA);
    flag_family=zeros(75,1);%75
    P=zeros(n,2);
    for j=1:n1(1)
                P(j,1)=norm(x-AA(j,:));
                P(j,2)=1;
    end
    for i=2:75
            num=sum(n1(1:i-1));
            for j=1:n1(i)
                P(num+j,1)=norm(x-AA(num+j,:));
                P(num+j,2)=i;
            end
    end
    aa=myQuickSort1(P,1,size(P,1));
    for mm=1:k
        flag_family(aa(mm,2))=flag_family(aa(mm,2))+1;
    end
end