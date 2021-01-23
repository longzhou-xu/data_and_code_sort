function [kappa,betak,center,mava,FNA_OUT,F_OUT] = kappa_lz(ava_size_serise,min_size,max_size,m)
%calculate the kappa
%input is the avalanches size serise
    max_size_alter=max(ava_size_serise);
    if max_size_alter<max_size
       max_size=max_size_alter;
    end
    betak=logspace(log10(min_size), log10(max_size), m);
    ava_size_serise(ava_size_serise>max_size)=[];
    ava_size_serise(ava_size_serise<min_size)=[];
    [size_hist_logspace,center]=hist(ava_size_serise,betak);
    mava = size_hist_logspace./sum(size_hist_logspace);

    L=max_size;
    size_hist=hist(ava_size_serise,min_size:max_size);
    size_hist_p=size_hist./sum(size_hist);
    syms x;
    equ = (sqrt(max_size)/(sqrt(max_size)-sqrt(x)))*(1-sqrt(x)) == size_hist_p(1);
    m_FNA=solve(equ,x);
    l=eval(m_FNA);
    C = sqrt(L)./(sqrt(L)-sqrt(l));
    
    kappa=0;

    for j = 1:length(betak)
        a=sqrt(betak(j));
        b=sqrt(l);
        FNA = C*(1-b./a);
        F = sum(mava(1:j));
        kappa = kappa + FNA - F;
        FNA_OUT(j)=FNA;
        F_OUT(j)=F;
    end
    kappa = 1 + kappa / m;
    
end