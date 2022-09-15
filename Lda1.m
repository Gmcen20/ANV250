function coeff = Lda1(training,group)
    class = sort(unique(group));%group��ѵ������ÿһ�������ı�ǩ��class�ǲ�ͬ����
    C = arrayfun(@(x)training(group==x,:),class,'UniformOutput',false);
    %����arrayfun�������Ա�����ν��ѭ�����Ӷ������ߴ���ļ���ԡ�
    mu  = mean(training);%��������ƽ��ֵ
    col = size(training,2);
    Sb = zeros(col,col);%���ɢ�Ⱦ���
    Sw = zeros(col,col);%����ɢ�Ⱦ���
    for i = class'
        c = C{i};
        x = mean(c) - mu;
        Sb = Sb + x'*x*length(c);
        Sw = Sw - cov(c)*length(training);
    end
    [V,D] = eigs(Sw\Sb);%��Sw\Sb���������ֽ�,V������ֵ���Ӵ�С���У���D������������ɵľ���
    coeff=real(V(:,1:2));%����������ȡʵ��������������Ͷ��2ά�У�����ȡǰ������������ֵ����
end