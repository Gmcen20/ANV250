% name=importdata('name6.txt');
% for i=1:74
%     name0{i}=name{75-i};
% end
% name=importdata('name_100dcsv.txt');
% name1=replace(name0,".csv","");
flag_lowdim250correct=zeros(74,74);
for i=1:73
    A=importdata(name0{i});
    mutset0=A;
    for j=i+1:74
    A=importdata(name0{j});
    mutset00=A;
try
    [coh0x,coh00x,coh1y,coh2y,coeff0_00]=MMC_test(mutset0,mutset00);
A=[coh0x,coh1y];
B=[coh00x,coh2y];
flag_lowdim250correct(i,j)=intersection1(A,B);
catch
    continue;
end
     end
end
csvwrite('flag_highdim_of_250d.csv',flag_lowdim250correct)