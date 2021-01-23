%% power law analysis robustness for HY96 LMS
clc; clear; close all
addpath('CritAnalysisSoftwarePackage2016-04-25')
load('avalanches_LMS_1.4.mat', 'ava_size_LMS')
load('avalanches_LMS_1.4.mat', 'ava_duration_LMS')
AvaS_LMS = ava_size_LMS;
AvaD_LMS = ava_duration_LMS;
for Size_min = 1 : 10
    for Size_max = 30 : 60
        I = Size_min
        J = Size_max-29
        [Alpha(I, J), Smin2fit(I, J), Smax2fit(I, J)] = plmle(AvaS_LMS, 'xmin', Size_min, 'xmax', Size_max);
        [P_alpha(I, J)] = pvcalc(AvaS_LMS, Alpha(I, J), 'xmin',Smin2fit(I, J),'xmax',Smax2fit(I, J));
    end
end
[Lifetime, AverageSize]=lifetimeAverageSize(AvaD_LMS, AvaS_LMS);
for Dura_min = 1 : 5
    for Dura_max = 9 : 20
        X = Dura_min
        Y = Dura_max-8
        [Tau(X, Y), Tmin2fit(X, Y), Tmax2fit(X, Y)] = ...
            plmle(AvaD_LMS,'xmin', Dura_min, 'xmax', Dura_max);
        
        [P_tau(X, Y)] = ...
            pvcalc(AvaD_LMS, Tau(X, Y), 'xmin',Tmin2fit(X, Y),'xmax',Tmax2fit(X, Y));
        
        [Gamma_nlsfit] = ...
            createFit_pl_nonlinearLS(Lifetime(Tmin2fit(X, Y):Tmax2fit(X, Y)), AverageSize(Tmin2fit(X, Y):Tmax2fit(X, Y)));
        
        Gamma_fit(X, Y) = Gamma_nlsfit.b;
        Gamma_fit_a(X, Y) = Gamma_nlsfit.a;
    end
end
save('HY96_powerlaw_robust_1.4_LMS.mat');
%% power law analysis robustness for HY96 MMS
clc; clear; close all
addpath('CritAnalysisSoftwarePackage2016-04-25')
load('avalanches_MMS_1.4.mat', 'ava_size_MMS')
load('avalanches_MMS_1.4.mat', 'ava_duration_MMS')
AvaS_MMS = ava_size_MMS;
AvaD_MMS = ava_duration_MMS;
for Size_min = 1 : 10
    for Size_max = 30 : 60
        I = Size_min
        J = Size_max-29
        [Alpha(I, J), Smin2fit(I, J), Smax2fit(I, J)] = plmle(AvaS_MMS, 'xmin', Size_min, 'xmax', Size_max);
        [P_alpha(I, J)] = pvcalc(AvaS_MMS, Alpha(I, J), 'xmin',Smin2fit(I, J),'xmax',Smax2fit(I, J));
    end
end
[Lifetime, AverageSize]=lifetimeAverageSize(AvaD_MMS, AvaS_MMS);
for Dura_min = 1 : 5
    for Dura_max = 9 : 20
        X = Dura_min
        Y = Dura_max-8
        [Tau(X, Y), Tmin2fit(X, Y), Tmax2fit(X, Y)] = ...
            plmle(AvaD_MMS,'xmin', Dura_min, 'xmax', Dura_max);
        
        [P_tau(X, Y)] = ...
            pvcalc(AvaD_MMS, Tau(X, Y), 'xmin',Tmin2fit(X, Y),'xmax',Tmax2fit(X, Y));
        
        [Gamma_nlsfit] = ...
            createFit_pl_nonlinearLS(Lifetime(Tmin2fit(X, Y):Tmax2fit(X, Y)), AverageSize(Tmin2fit(X, Y):Tmax2fit(X, Y)));
        
        Gamma_fit(X, Y) = Gamma_nlsfit.b;
        Gamma_fit_a(X, Y) = Gamma_nlsfit.a;
    end
end
save('HY96_powerlaw_robust_1.4_MMS.mat');
%% power law analysis robustness for HY96 HMS
clc; clear; close all
addpath('CritAnalysisSoftwarePackage2016-04-25')
load('avalanches_HMS_1.4.mat', 'ava_size_HMS')
load('avalanches_HMS_1.4.mat', 'ava_duration_HMS')
AvaS_HMS = ava_size_HMS;
AvaD_HMS = ava_duration_HMS;
for Size_min = 1 : 10
    for Size_max = 30 : 60
        I = Size_min
        J = Size_max-29
        [Alpha(I, J), Smin2fit(I, J), Smax2fit(I, J)] = plmle(AvaS_HMS, 'xmin', Size_min, 'xmax', Size_max);
        [P_alpha(I, J)] = pvcalc(AvaS_HMS, Alpha(I, J), 'xmin',Smin2fit(I, J),'xmax',Smax2fit(I, J));
    end
end
[Lifetime, AverageSize]=lifetimeAverageSize(AvaD_HMS, AvaS_HMS);
for Dura_min = 1 : 5
    for Dura_max = 9 : 20
        X = Dura_min
        Y = Dura_max-8
        [Tau(X, Y), Tmin2fit(X, Y), Tmax2fit(X, Y)] = ...
            plmle(AvaD_HMS,'xmin', Dura_min, 'xmax', Dura_max);
        
        [P_tau(X, Y)] = ...
            pvcalc(AvaD_HMS, Tau(X, Y), 'xmin',Tmin2fit(X, Y),'xmax',Tmax2fit(X, Y));
        
        [Gamma_nlsfit] = ...
            createFit_pl_nonlinearLS(Lifetime(Tmin2fit(X, Y):Tmax2fit(X, Y)), AverageSize(Tmin2fit(X, Y):Tmax2fit(X, Y)));
        
        Gamma_fit(X, Y) = Gamma_nlsfit.b;
        Gamma_fit_a(X, Y) = Gamma_nlsfit.a;
    end
end
save('HY96_powerlaw_robust_1.4_HMS.mat')