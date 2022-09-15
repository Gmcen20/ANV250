%name=importdata('name2.txt');
%name0=replace(name,"\","");
%name1=replace(name,".fasta\","");
cd('/Users/guanmengcen/Desktop/virus/giant virus')
s1=fastaread('Yersinia phage phiR1-37.fasta');
n1=length(s1)
%feature39=zeros(n1,12);
for m=1:n1
   N=length(s1(m).Sequence);
   miu=zeros(4,N);
   t=zeros(4,N);
   t0=zeros(N,1);
   theta=zeros(4,1);
   sigma=zeros(4,1);
   n=zeros(4,1);
   D=zeros(4,1);
   kesai=zeros(4,1);
   for i=1:N
   if s1(m).Sequence(i)=='A'
        n(1)=n(1)+1;
        t(1,n(1))=i;
        miu(1,i)=miu(1,i)+1;
        elseif s1(m).Sequence(i)=='T'
        n(2)=n(2)+1;
        t(2,n(2))=i;
        miu(2,i)=miu(2,i)+1;
        elseif s1(m).Sequence(i)=='C' 
        n(3)=n(3)+1;
        t(3,n(3))=i;
        miu(3,i)=miu(3,i)+1;
        elseif s1(m).Sequence(i)=='G'
        n(4)=n(4)+1;
        t(4,n(4))=i;
        miu(4,i)=miu(4,i)+1;
   end
   end
   for i=1:4
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
   for i=1:4
       for j=1:N
       theta(i)=theta(i)+miu(i,j)/N;
       end
   end
   for i=1:4
       for j=1:N
           sigma(i)=sigma(i)+miu(i,j);
       end
   end
   for i=1:4
       for j=1:N
           if n(i)>0
           D(i)=D(i)+(miu(i,j)-theta(i))*(miu(i,j)-theta(i))/(n(i)*n(i));
           else
               D(i)=0;
           end
       end
   end
   %feature39(m,:)=[n',theta',D']
end
Feature=[feature1;feature2;feature3;feature4;feature5;feature6;feature7;feature8;feature9;feature10;feature11;feature12;feature13;feature14;feature15;feature16;feature17;feature18;feature19;feature20;feature21;feature22;feature23;feature24;feature25;feature26;feature27;feature28;feature29;feature30;feature31;feature32;feature33;feature34;feature35;feature36;feature37;feature38;feature39]
c1=['Myoviridae_GV','12d.csv'];
csvwrite(c1,Feature)


