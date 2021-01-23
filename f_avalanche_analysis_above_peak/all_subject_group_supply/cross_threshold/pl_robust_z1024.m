%% Z1024 power law fit robustness for allsubject
clc; clear; close all
addpath('D:\supplyment_dataANDcode\z_function\CritAnalysisSoftwarePackage2016-04-25')
load('z_1024\Avalanches_3.2.mat')
for Size_min = 1 : 10
    for Size_max = 30 : 60
        I = Size_min
        J = Size_max-29
        [Alpha(I, J), Smin2fit(I, J), Smax2fit(I, J)] = plmle(AvaS_AllSubject, 'xmin', Size_min, 'xmax', Size_max);
        [P_alpha(I, J)] = pvcalc(AvaS_AllSubject, Alpha(I, J), 'xmin',Smin2fit(I, J),'xmax',Smax2fit(I, J));
    end
end
[Lifetime, AverageSize]=lifetimeAverageSize(AvaD_AllSubject, AvaS_AllSubject);
for Dura_min = 1 : 5
    for Dura_max = 9 : 20
        X = Dura_min
        Y = Dura_max - 8
        [Tau(X, Y), Tmin2fit(X, Y), Tmax2fit(X, Y)] = ...
            plmle(AvaD_AllSubject,'xmin', Dura_min, 'xmax', Dura_max);
        
        [P_tau(X, Y)] = ...
            pvcalc(AvaD_AllSubject, Tau(X, Y), 'xmin',Tmin2fit(X, Y),'xmax',Tmax2fit(X, Y));
        
        [Gamma_nlsfit] = ...
            createFit_pl_nonlinearLS(Lifetime(Tmin2fit(X, Y):Tmax2fit(X, Y)), AverageSize(Tmin2fit(X, Y):Tmax2fit(X, Y)));
        
        Gamma_fit(X, Y) = Gamma_nlsfit.b;
        Gamma_fit_a(X, Y) = Gamma_nlsfit.a;
    end
end
save('z_1024\Z1024_powerlaw_robust_3.2_allsubject.mat');