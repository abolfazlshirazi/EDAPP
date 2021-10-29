function [tab,dta] = SASS_BIN(run,I_fno)
Par    = Introd_Par(I_fno);
iter   = 0;          
stream = RandStream('mt19937ar','Seed',sum(100*clock)*run);
RandStream.setGlobalStream(stream);
const_num = Par.gn(I_fno) + Par.hn(I_fno); %% constraint number
%% define variables
current_eval=0;             %% current fitness evaluations
PS1 = Par.PopSize;            %% define PS1
%% ====================== Initalize x ==================================
x   = repmat(Par.xmin,Par.PopSize,1)+repmat((Par.xmax-Par.xmin),Par.PopSize,1).*rand(Par.PopSize,Par.n);
%% calc. fit. and update FES
[fx,gx,hx]   = cec20_func(x,I_fno); fx = fx';
current_eval = current_eval + Par.PopSize;
convx = (sum([sum(gx.*(gx>0),1); 
sum(abs(hx).*(abs(hx)>1e-4),1)],1)./const_num);
%% ================== Initalization of population for SS ===================================
EA_1= x(1:PS1,:);    EA_obj1= fx(1:PS1);   
EA_gx1 = gx(:,1:PS1); EA_hx1 = hx(:,1:PS1); EA_cx1 = convx(1:PS1);
%%%
Pop.x = EA_1';
Pop.f = EA_obj1;
Pop.g = EA_gx1;
Pop.h = EA_hx1;
Pop.conv = EA_cx1;
%% ================== Initalization of \epsilon-constrained ================================
[TC,Eg,Eh,CPg,CPh] = ETAintialization(Pop);
 Eg0 = Eg;
 Eh0 = Eh;  
 

%% ====================== store the best ==================
ranking = eta_sort(Pop,Eg,Eh);
best_sol = Pop.x(:,ranking(1));
best_of  = Pop.f(ranking(1));
best_conv = Pop.conv(ranking(1));
%% ===================== archive data ====================================
arch_rate = 2.6;
archive.NP = arch_rate * PS1; % the maximum size of the archive
archive.pop = zeros(0, Par.n); % the solutions stored in te archive
archive.funvalues = zeros(0, 1); % the function value of the archived solutions
%% ==================== to adapt ci and ri =================================
hist_pos    = 1;
memory_size = 6;
archive_c  = ones(1,memory_size).*0.5;
archive_rp  = ones(1,memory_size).*0.5;
%%
stop_con=0; 
flag = zeros(1,10);
flag1 = 0;
flag2 = 0;
%% main loop
while stop_con==0
    %% update Eg and Eh
        if(iter>1 && iter<TC)
          Eg=Eg0.*((1-iter./TC).^CPg);
          Eh=Eh0.*(1-iter./TC).^CPh;
        elseif(iter+1>=TC)
          Eg = 0;
          Eh = 0;
        end 
    iter=iter+1;
    %% ====================== SS main algorithm ============================
    if (current_eval<Par.Max_FES)
            %% apply MODE
            [EA_1, EA_obj1,EA_gx1, EA_hx1, EA_cx1,best_sol,best_of,best_conv,archive,hist_pos,memory_size, archive_c,archive_rp, current_eval] = ...
                SASS( EA_1, EA_obj1,EA_gx1, EA_hx1, EA_cx1,best_sol, best_of, best_conv,archive,hist_pos,memory_size, archive_c,archive_rp,....
                 Par.xmin, Par.xmax,  Par.n,  PS1,  current_eval, I_fno, iter, const_num, Eg, Eh, Par);
    end
    %% ====================== Repair best solution ========================
    problem.constr_fun_name = @cec20_func;
    problem.gn              = Par.gn;
    problem.hn              = Par.hn;        
    problem.I_fno           = I_fno;
    problem.xmin            = Par.xmin;
    problem.xmax            = Par.xmax;
    problem.n               = Par.n;
    problem.maxiter         = 5000;
    if  best_conv ~= 0 && mod(iter,Par.n) == 0 && Par.n >= 20 && problem.hn(I_fno) > 0 
          [~,gx,hx]           = cec20_func(best_sol(:)',I_fno); 
          [new_mutant,fes]     = QNrepair(problem, best_sol', gx, hx);
          best_new_sol        = han_boun(new_mutant', Par.xmax, Par.xmin, new_mutant',1,3);
          current_eval        = current_eval+fes;
          [fval, gv, hv]      = cec20_func(best_new_sol,I_fno);
          best_new_of         = fval;                                                               
          best_new_conv       = (sum([sum(gv.*(gv>0),1); sum(abs(hv).*(abs(hv)>1e-4),1)])./const_num);
          if ((best_of-best_new_of)>0 && best_conv == best_new_conv)|| best_new_conv < best_conv
              best_sol  = best_new_sol;
              best_of   = best_new_of;
              best_conv = best_new_conv;
          end
    end 
    %% ====================== stopping criterion check ====================
    if (current_eval>=Par.Max_FES)
        stop_con=1;
    end
    %% calculation of Data
    if flag(1) == 0 && current_eval >= 0.1*Par.Max_FES
       dta(1,:) = [best_of best_conv];
       flag(1) = 1;
    elseif flag(2) == 0 && current_eval >= 0.2*Par.Max_FES
       dta(2,:) = [best_of best_conv];
       flag(2) = 1; 
    elseif flag(3) == 0 && current_eval >= 0.3*Par.Max_FES
       dta(3,:) = [best_of best_conv];
       flag(3) = 1;
    elseif flag(4) == 0 && current_eval >= 0.4*Par.Max_FES
       dta(4,:) = [best_of best_conv];
       flag(4) = 1; 
    elseif flag(5) == 0 && current_eval >= 0.5*Par.Max_FES
       dta(5,:) = [best_of best_conv];
       flag(5) = 1;
    elseif flag(6) == 0 && current_eval >= 0.6*Par.Max_FES
       dta(6,:) = [best_of best_conv];
       flag(6) = 1;
    elseif flag(7) == 0 && current_eval >= 0.7*Par.Max_FES
       dta(7,:) = [best_of best_conv];
       flag(7) = 1; 
    elseif flag(8) == 0 && current_eval >= 0.8*Par.Max_FES
       dta(8,:) = [best_of best_conv];
       flag(8) = 1;
    elseif flag(9) == 0 && current_eval >= 0.9*Par.Max_FES
       dta(9,:) = [best_of best_conv];
       flag(9) = 1; 
    elseif flag(10) == 0 && current_eval >= 1*Par.Max_FES
       dta(10,:) = [best_of best_conv];
       flag(10) = 1;
    end
    %% calculation of Tab
    % log global best after having used 10%, and 50% of the evaluation budget
        if current_eval>= Par.Max_FES*10/100 && flag1==0
            fit10=best_of;
            con10=best_conv;
            [ff,gg,hh]=feval(@cec20_func,best_sol,I_fno);
            c10_1    = sum(gg>1)                  + sum(abs(hh)>1);
            c10_2    = sum((gg>0.01) & (gg<1))    + sum(abs(hh)>0.01 & abs(hh)<1);
            c10_3    = sum((gg>0.0001)&(gg<0.01)) + sum(abs(hh)>0.0001 & abs(hh)<0.01);  
            flag1=1;
        elseif current_eval>=Par.Max_FES*50/100 && flag2==0
            fit50=best_of;
            con50=best_conv;
            [ff,gg,hh]=feval(@cec20_func,best_sol,I_fno);
            c50_1    = sum(gg>1)                  + sum(abs(hh)>1);
            c50_2    = sum((gg>0.01)&(gg<1))      + sum(abs(hh)>0.01 & abs(hh)<1);
            c50_3    = sum((gg>0.0001)&(gg<0.01)) + sum(abs(hh)>0.0001 & abs(hh)<0.01)  ;
            flag2=1;
        end
%% diversity
if (std(EA_cx1) < 1e-8 && min(EA_cx1) > 0) || (std(EA_obj1) < 1e-8)
%% ====================== Initalize x ==================================
x   = repmat(Par.xmin,Par.PopSize,1)+repmat((Par.xmax-Par.xmin),Par.PopSize,1).*rand(Par.PopSize,Par.n);
%% calc. fit. and update FES
[fx,gx,hx] = cec20_func(x,I_fno); fx = fx';
current_eval = current_eval+Par.PopSize;
convx = (sum([sum(gx.*(gx>0),1); sum(abs(hx).*(abs(hx)>1e-4),1)],1)./const_num);
%% ================== fill in for each  individual ===================================
%% SS
EA_1= x(1:PS1,:);    EA_obj1= fx(1:PS1); 
EA_gx1 = gx(:,1:PS1); EA_hx1 = hx(:,1:PS1); EA_cx1 = convx(1:PS1);           
end
end
%% dta size
while size(dta,1) < 10
    dta(end+1,:) = dta(end,:);
end
    fit100=best_of;
    con100=best_conv;
    [ff,gg,hh]=feval(@cec20_func,best_sol,I_fno);
    c100_1    = sum(gg>1)                   + sum(abs(hh)>1);
    c100_2    = sum((gg>0.01)&(gg<1))       + sum(abs(hh)>0.01 & abs(hh)<1);
    c100_3    = sum((gg>0.0001)&(gg<0.01))  + sum(abs(hh)>0.0001 &abs(hh)<0.01);


    tab = [fit10 con10 c10_1 c10_2 c10_3;
             fit50 con50 c50_1 c50_2 c50_3;
             fit100 con100 c100_1 c100_2 c100_3];
end
