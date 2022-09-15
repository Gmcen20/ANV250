name=importdata('name7.txt');
name0=replace(name,"\","");
name0=replace(name0,"f0fs20 cf0 adenoviridae270d.csv","adenoviridae270d.csv");
name0=replace(name0,"cf0 Virgaviridae270d.csv","Virgaviridae270d.csv");
name1=replace(name,".csv\","");
flag_lowdim250pcalda=zeros(75,75);
for i=1:74
        mutset0=importdata(name0{i+7});
      for j=i+1:74
        mutset00=importdata(name0{j+7});
try
[coh0x,coh00x,coh1y,coh2y,coeff0_00]=MMC_test(mutset0,mutset00);
A=[coh0x,coh1y];
B=[coh00x,coh2y];
flag_lowdim250pcalda(i,j)=intersection1(A,B);
catch
    continue;
end
     end
end

for i=1:74
        mutset0=importdata(name0{i+7});
      for j=76
        mutset00=importdata(name0{j+7});
try
[coh0x,coh00x,coh1y,coh2y,coeff0_00]=MMC_test(mutset0,mutset00);
A=[coh0x,coh1y];
B=[coh00x,coh2y];
flag_lowdim250pcalda(i,j-1)=intersection1(A,B);
catch
    continue;
end
     end
end
csvwrite('flag_lowdim_of_270.csv',flag_lowdim250pcalda)