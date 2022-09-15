function coeff = Lda1(training,group)
    class = sort(unique(group));%group是训练集中每一个样本的标签，class是不同的类
    C = arrayfun(@(x)training(group==x,:),class,'UniformOutput',false);
    %利用arrayfun函数可以避免无谓的循环，从而大大提高代码的简介性。
    mu  = mean(training);%计算样本平均值
    col = size(training,2);
    Sb = zeros(col,col);%类间散度矩阵
    Sw = zeros(col,col);%类内散度矩阵
    for i = class'
        c = C{i};
        x = mean(c) - mu;
        Sb = Sb + x'*x*length(c);
        Sw = Sw - cov(c)*length(training);
    end
    [V,D] = eigs(Sw\Sb);%对Sw\Sb进行特征分解,V是特征值（从大到小排列），D是特征向量组成的矩阵
    coeff=real(V(:,1:2));%对特征向量取实部，而且由于想投到2维中，所有取前两个最大的特征值就行
end