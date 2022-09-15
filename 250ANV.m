name=importdata('name0.txt');
name0=replace(name,"\","");
name0=replace(name0,"Virgaviridae.fasta}","Virgaviridae.fasta");
name0=replace(name0,"f0fs24 cf0 adenoviridae.fasta","adenoviridae.fasta");
name1=replace(name,".fasta\","");
for ww=61:61
s1=fastaread(name0{ww+6});
n1=length(s1);
feature1=zeros(n1,250);
for m=1:n1
    flagg=1;
   N=length(s1(m).Sequence);
   miu=zeros(20,N);
   t=zeros(20,N);
   t0=zeros(N,1);
   theta=zeros(20,1);
   sigma=zeros(20,1);
   n=zeros(20,1);
   D=zeros(20,1);
   kesai=zeros(20,1);
   cov=zeros(20,20);
   for i=1:N
   if s1(m).Sequence(i)=='A'
        n(1)=n(1)+1;
        t(1,n(1))=i;
        miu(1,i)=miu(1,i)+1;
        elseif s1(m).Sequence(i)=='R'
        n(2)=n(2)+1;
        t(2,n(2))=i;
        miu(2,i)=miu(2,i)+1;
        elseif s1(m).Sequence(i)=='N' 
        n(3)=n(3)+1;
        t(3,n(3))=i;
        miu(3,i)=miu(3,i)+1;
        elseif s1(m).Sequence(i)=='D'
        n(4)=n(4)+1;
        t(4,n(4))=i;
        miu(4,i)=miu(4,i)+1;
        elseif s1(m).Sequence(i)=='C'
        n(5)=n(5)+1;
        t(5,n(5))=i;
        miu(5,i)=miu(5,i)+1;
        elseif s1(m).Sequence(i)=='Q'
        n(6)=n(6)+1;
        t(6,n(6))=i;
        miu(6,i)=miu(6,i)+1;
        elseif s1(m).Sequence(i)=='E'
        n(7)=n(7)+1;
        t(7,n(7))=i;
        miu(7,i)=miu(7,i)+1;
        elseif s1(m).Sequence(i)=='G'
        n(8)=n(8)+1;
        t(8,n(8))=i;
        miu(8,i)=miu(8,i)+1;
         elseif s1(m).Sequence(i)=='H'
        n(9)=n(9)+1;
        t(9,n(9))=i;
        miu(9,i)=miu(9,i)+1;
         elseif s1(m).Sequence(i)=='I'
        n(10)=n(10)+1;
         t(10,n(10))=i;
        miu(10,i)=miu(10,i)+1;
         elseif s1(m).Sequence(i)=='L'
        n(11)=n(11)+1;
         t(11,n(11))=i;
        miu(11,i)=miu(11,i)+1;
         elseif s1(m).Sequence(i)=='K'
        n(12)=n(12)+1;
         t(12,n(12))=i;
        miu(12,i)=miu(12,i)+1;
         elseif s1(m).Sequence(i)=='M'
        n(13)=n(13)+1;
         t(13,n(13))=i;
        miu(13,i)=miu(13,i)+1;
         elseif s1(m).Sequence(i)=='F'
        n(14)=n(14)+1;
         t(14,n(14))=i;
        miu(14,i)=miu(14,i)+1;
        elseif s1(m).Sequence(i)=='P'
        n(15)=n(15)+1;
         t(15,n(15))=i;
        miu(15,i)=miu(15,i)+1;
        elseif s1(m).Sequence(i)=='S'
        n(16)=n(16)+1;
         t(16,n(16))=i;
        miu(16,i)=miu(16,i)+1;
        elseif s1(m).Sequence(i)=='T'
        n(17)=n(17)+1;
         t(17,n(17))=i;
        miu(17,i)=miu(17,i)+1;
        elseif s1(m).Sequence(i)=='W'
        n(18)=n(18)+1;
         t(18,n(18))=i;
        miu(18,i)=miu(18,i)+1;
        elseif s1(m).Sequence(i)=='Y'
        n(19)=n(19)+1;
         t(19,n(19))=i;
        miu(19,i)=miu(19,i)+1;
        elseif s1(m).Sequence(i)=='V'
        n(20)=n(20)+1;
         t(20,n(20))=i;
        miu(20,i)=miu(20,i)+1;
       else
      flagg=0;
       continue;
   end
   end
   if flagg==0 
       m
       continue;
   end
   for i=1:20
       if n(i)>0&&t(i,1)-1>0
       for j=1:t(i,1)-1
           miu(i,j)=0;
       end
       for j=2:n(i)
           for k=t(i,j-1):(t(i,j)-1)
               miu(i,k)=j-1;
           end
       end
       for k=t(i,n(i)):N
           miu(i,k)=n(i);
       end
       elseif n(i)>0&&t(i,1)-1==0
           for j=2:n(i)
           for k=t(i,j-1):(t(i,j)-1)
               miu(i,k)=j-1;
           end
           end
       for k=t(i,n(i)):N
           miu(i,k)=n(i);
       end
       end
   end
   for i=1:20
       for j=1:N
       theta(i)=theta(i)+miu(i,j)/N;
       end
   end
   for i=1:20
       for j=1:N
           sigma(i)=sigma(i)+miu(i,j);
       end
   end
   for i=1:20
       for j=1:N
           if n(i)>0
           D(i)=D(i)+(miu(i,j)-theta(i))*(miu(i,j)-theta(i))/(n(i)*n(i));
           else
               D(i)=0;
           end
       end
   end
   for i=1:20
       if n(i)>0
       kesai(i)=sigma(i)/n(i);
       else 
           kesai(i)=0;
       end
   end
   for i=1:20
       for j=1:20
           if i>j
               for k=1:N
                   if n(i)>0&&n(j)>0
               cov(i,j)=cov(i,j)+(miu(i,k)-theta(i))*(miu(j,k)-theta(j))/(n(i)*n(j));
                   else
                      cov(i,j)=0; 
                   end
               end
           else
               cov(i,j)=0;
           end
       end
   end
   for i=1:20
       feature1(m,i)=n(i);
   end
   for i=21:40
       feature1(m,i)=kesai(i-20);
   end
   for i=41:60
       feature1(m,i)=D(i-40);
   end
   mm=0;
   for i=1:20
       for j=1:20
           if i>j
               mm=mm+1;
           feature1(m,mm+60)=cov(i,j);
           end
       end
   end
end
feature1(all(feature1==0,2),:)=[];
c=[name1{ww+6},'_cor250d.csv'];
end
