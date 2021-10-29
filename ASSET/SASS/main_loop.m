%%----------------------------------------------------------------------------------------------------%%
% Authors: Abhishek Kumar, Swagatam Das, and Ivan Zelinka
% Title: A Self-Adaptive Spherical Search Algorithm for Real-World Constrained Optimization Problems
%%----------------------------------------------------------------------------------------------------%%

clc;
format short;
global  initial_flag
global initial_flag2
%%  introductory Definitions
num_prbs = 57;                      %% number of test problems
max_runs = 25;                      %% number of runs
outcome=zeros(max_runs,1);          %% to save the solutions of each run
com_time=zeros(max_runs,1);         %% Computational time
SR=zeros(max_runs,1);               %% How many times the optimal solution is obtained
Avg_FES=zeros(max_runs,1);          %% average fitness evaluations to reach f*+1e08
Final_results=zeros(num_prbs,8);    %% to save the final results
warning('off');
%% ========================= main loop ====================================
for I_fno = 1:57
    initial_flag = 0;
    initial_flag2 = 0;
    Par= Introd_Par(I_fno);
    
    for run=1:max_runs
        
       [tab,dta] = SASS_BIN(run,I_fno);%tab(abs(tab) < 1e-8) = 0;
        OB(:,run) = dta(:,1);
        CONV(:,run) = dta(:,2);
        FitT(run,:)=[tab(1,:) tab(2,:) tab(3,:)];
        %% to print the convergence of ech run % set 0 if not
        if Par.disp==1
           disp([I_fno, run, OB(10,run), CONV(10,run)]);
        end
    end
    if Par.Printing == 1
       filename1 = strcat(['Result_Record\SASS_RC' num2str(I_fno) '_F' '.txt']);
       filename2 = strcat(['Result_Record\SASS_RC' num2str(I_fno) '_CV' '.txt']);
       save(filename1,'OB','-ascii','-double','-tabs');
       save(filename2,'CONV','-ascii','-double','-tabs');
       Tab=build_stats(FitT,max_runs);
       Stats = xlsread('Result_Record\SASS_temp.xlsx',1); Stats = Stats';
       Stats(I_fno,:)=[I_fno Tab];
       disp([I_fno Tab]);
       xlswrite('Result_Record\SASS_temp.xlsx', Stats', 1 );
    end
end


