function [coh1x,coh2x,coh1y,coh2y,coeff2]=LDA_test(mutsetA,mutsetB)
        n1 = size(mutsetA,1);
        n2 = size(mutsetB,1);
        training = [mutsetA; mutsetB];
        one = zeros(n1,1)+1;
        two = zeros(n2,1) + 2;
        group = [one;two];
        coeff2 = Lda1(training, group);
        coh1x = mutsetA * coeff2(:,2);
	    coh2x = mutsetB * coeff2(:,2);
        coh1y = mutsetA * coeff2(:,1);
        coh2y = mutsetB * coeff2(:,1);       
end
        