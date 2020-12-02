function [Par] = Introd_Par(I_fno)
%%
Par = Cal_par(I_fno);  

%% loading
if Par.n <= 10
        max_nfes = 1e5;
elseif Par.n > 10 && Par.n <= 30
        max_nfes = 2e5;
elseif Par.n > 30 && Par.n <= 50
        max_nfes = 4e5;
elseif Par.n > 50 && Par.n <=150
        max_nfes = 8e5;
else
        max_nfes = 1e6;
end
Par.Max_FES = max_nfes;

   Par.PopSize = 60; 

%% printing the detailed results- this will increase the computational time
Par.Printing=1; %% 1 to print; 0 otherwise
Par.disp = 1;
end