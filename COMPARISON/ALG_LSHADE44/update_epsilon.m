function [epsilon]=update_epsilon(epsilon,fea_percent,alpha,cp,Tc,ratio_FES, max_cv)
if epsilon==inf && fea_percent<1
    epsilon=max_cv;
end
if ratio_FES<Tc
    if fea_percent<alpha
        epsilon = epsilon*(1-ratio_FES/Tc)^cp;
    else
        epsilon=1.1*max_cv;
    end
else
    epsilon=0;
end
end