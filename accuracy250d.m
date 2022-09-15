name=importdata('name3.txt');
name0=replace(name,"\","");
n=zeros(75,1);
n_train=zeros(75,1);
n_test=zeros(75,1);
A=zeros(16010,250);
R=zeros(10727,250);
T=zeros(5283,250);
A1=importdata('Virgaviridae250d.csv');
[n(1),m1]=size(A1);
n_train(1)=round(0.67*n(1));
n_test(1)=n(1)-n_train(1);
A(1:n(1),:)=A1;
R(1:n_train(1),:)=A1(1:round(0.67*n(1)),:);
T(1:n_test(1),:)=A1(round(0.67*n(1))+1:n(1),:);
for i=2:74
    n0=sum(n);
    n00=sum(n_train);
    n000=sum(n_test);
    A2=importdata(name0{i+6});
    [n(i),m1]=size(A2);
    A(n0+1:sum(n),:)=A2;
    n_train(i)=round(0.67*n(i));
    n_test(i)=n(i)-n_train(i);
    R(n00+1:sum(n_train),:)=A2(1:round(0.67*n(i)),:);
    T(n000+1:sum(n_test),:)=A2(round(0.67*n(i))+1:n(i),:);
end
A75=importdata('Chrysoviridae250d.csv');
    n0=sum(n);
    n00=sum(n_train);
    n000=sum(n_test);
[n(75),m1]=size(A75);
n_train(75)=round(0.67*n(75));
n_test(75)=n(75)-n_train(75);
A(n0+1:sum(n),:)=A75;
R(n00+1:sum(n_train),:)=A75(1:round(0.67*n(75)),:);
T(n000+1:sum(n_test),:)=A75(round(0.67*n(75))+1:n(75),:);
flag=zeros(75,75);
index=1;
for i=1:n_test(index)
    dis=zeros(75,1);
    dis(1)=distance_convex(T(i,:)',R(1:n_train(1),:));
    for j=2:75
        N1=sum(n_train(1:j-1));
        dis(j)=distance_convex(T(i,:)',R(N1+1:N1+n_train(j),:));
    end
    if dis(index)==min(dis)
        flag(index,i)=1;
    end
end
accuracy=zeros(75,1);

for index=2:75
N0=sum(n_test(1:index-1));
for i=1:n_test(index)
    dis=zeros(75,1);
    dis(1)=distance_convex(T(N0+i,:)',R(1:n_train(1),:));
    for j=2:75
        N1=sum(n_train(1:j-1));
        dis(j)=distance_convex(T(N0+i,:)',R(N1+1:N1+n_train(j),:));
    end
    if dis(index)==min(dis)
        flag(index,i)=1;
    end
end
accuracy=zeros(75,1);
csvwrite('flag250daccu.csv',flag)
index
end
for i=1:75
    accuracy(i)=sum(flag(i,:))/n_test(i);
end
mean(accuracy)
