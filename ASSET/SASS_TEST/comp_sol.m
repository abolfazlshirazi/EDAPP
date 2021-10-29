function I = comp_sol(f_new,f,g_new,g,h_new,h,Eg,Eh)
%% old solution
G = sum(max(0,g-Eg*ones(size(g))),1);
H = sum(max(0,abs(h)-(Eh*ones(size(h))+0.0001)),1);
conv = (G+H);
%% new solution
G_new = sum(max(0,g_new-Eg*ones(size(g))),1);
H_new = sum(max(0,abs(h_new)-(Eh*ones(size(h))+0.0001)),1);
conv_new = (G_new+H_new);
%% comparision
for i = 1:length(f)
    I(i) = eps_rank(f_new(i),conv_new(i),f(i),conv(i),0);
end
end