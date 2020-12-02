function [Sol_arr,FE_arr] = LSHADE44_run(I_fno)

global initial_flag;
initial_flag=0;

%%%%%%%%%%%%%%%%%%%%%%%%%% A. Shirazi: FE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%    
global FCounter;
FCounter = 0;
Sol_arr = [];
FE_arr = [];
   
setup = Introd_Par(I_fno);
    nvar  = setup.n;
%     lb    = setup.xmin;
%     ub    = setup.xmax;
    FEMax = setup.Max_FES;
    const_num = setup.gn(I_fno) + setup.hn(I_fno);
    
% max_run=25;
% Dimension_set=[10,30,50,100];
% I_fno=1;        %question number
ini_CR=0.5;             
ini_F=0.5;
%%%%%%%%%%%%%%% The directory to save results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
% file_path = strcat(pwd,'/Results/');
% if ~isdir(file_path)
%     mkdir(file_path);
% end
%%%%%%%%%%%%%%% Strategy selection parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
n0=2;
K=4;
delta=1/K/5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for Dim_num = 1:4


%     D=Dimension_set(Dim_num);  %variable dimension
    D=nvar;  %variable dimension
    
    
%     for run=1:max_run
        
        strategy_prob=ones(K,1)/K;
        strategy_success=zeros(1,K);
%         output_rate=0;
        %%%%% initialize MCR & MF  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        H=10;            %Archive length
        H_label=zeros(1,K);        %Archive location
        H_length=zeros(1,K);
        MCR_archive=zeros(H,K);
        MF_archive=zeros(H,K);
        flag_success=zeros(1,K);   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ini_popsize=D*5;  %maximal population size
        min_popsize=5;    %minimal population size
        popsize=ini_popsize;
        Max_FES=FEMax;
        
%         lu = decision_range(I_fno,D)';        
        lu = [[setup.xmin]' [setup.xmax]']';
        
        %%%%%%%%%%    Parameters for epsilon  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Tc=0.8;             
        cp=2;               
        alpha=0.8;        
        theta= 0.2*popsize; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = ones(popsize, 1) * lu(1, :) + rand(popsize, D) .* (ones(popsize, 1) * (lu(2, :) - lu(1, :))); 
        [f,g,h]=cec20_func(x,I_fno);
        g=g';
        h=h';
        
        cv=overall_cv(g,h);
        population=[x,g,h,f,cv];
        population=epsilon_sort(population,0);
        result=population(1,:);
        FES=popsize;
        %%%%%%%%%%%  Initialize epsilon  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        max_cv = -1e30;
        infea_label=population(:,end)>0;
        cv=population(infea_label,end);
        fea_percent=1-sum(infea_label)/popsize; %feasibility ratio
        temp_infea_cv=sort(population(:,end),'descend');
        if temp_infea_cv(1) > max_cv
            max_cv = temp_infea_cv(1);
        end
        if temp_infea_cv(theta)>0
            epsilon=temp_infea_cv(theta);
        else
            epsilon=inf;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while FES<Max_FES
            
            population=epsilon_sort(population,0);
            population=population(1:popsize,:);
            temp_population=population;
            MCR=zeros(K,1); 
            MF=zeros(K,2);  
            label_success=zeros(1,K);      
            w_sum=zeros(1,K);
            if Max_FES-FES<popsize
                count=Max_FES-FES;
            else
                count=popsize;
            end
            for s=1:count

                indiv=population(s,:);
        %%%%%%%%%%%%%%% select strategy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
                num_strategy=select_strategy(strategy_prob);
        %%%%%%%%%%%%%%% Product F and CR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if flag_success(num_strategy)==0
                    CR=CR_product(ini_CR);
                    F=F_product(ini_F);
                else
                    ri=randi(H_length(num_strategy));
                    CR=CR_product(MCR_archive(ri,num_strategy));
                    F=F_product(MF_archive(ri,num_strategy));
                end
        %%%%%%%%%%%%%%%%  Product an offspring solution y %%%%%%%%%%%%%%%%%%%%%%%%%
                y=select_reproduct(indiv,s,lu,D,CR,F,temp_population,popsize,num_strategy);
                [f_y,g_y,h_y]=cec20_func(y,I_fno);
                
                
                g_y = g_y';
                h_y = h_y';
                cv=overall_cv(g_y,h_y);
                Y=[y,g_y,h_y,f_y,cv];
                FES=FES+1;
                
                if max_cv < cv
                    max_cv = cv;
                end
        %%%%%%%%%%%%%%%%%% Update the population %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if max(0,temp_population(s,end)-epsilon)>max(0,Y(end)-epsilon)
                    population(s,:)=Y;
                    [strategy_success,strategy_prob]=strategy_prob_update(strategy_success,strategy_prob,num_strategy,K,n0,delta);
                        w_k=temp_population(s,end)-Y(end);
                        label_success(num_strategy)=1;
                        w_sum(num_strategy)=w_k+w_sum(num_strategy);
                        MCR(num_strategy)=MCR(num_strategy)+CR*w_k;
                        MF(num_strategy,1)=MF(num_strategy,1)+F^2*w_k;
                        MF(num_strategy,2)=MF(num_strategy,2)+F*w_k;
                    if max(0,result(end))>max(0,Y(end))
                        result=Y;
                    elseif max(0,result(end))==max(0,Y(end))&&result(end-1)>Y(end-1)
                        result=Y;
                    end
                end
                if (max(0,temp_population(s,end)-epsilon)==max(0,Y(end)-epsilon))&& (temp_population(s,end-1)>Y(end-1))
                    population(s,:)=Y;
                    [strategy_success,strategy_prob]=strategy_prob_update(strategy_success,strategy_prob,num_strategy,K,n0,delta);
                        w_k=temp_population(s,end-1)-Y(end-1);
                        label_success(num_strategy)=1;
                        w_sum(num_strategy)=w_k+w_sum(num_strategy);
                        MCR(num_strategy)=MCR(num_strategy)+CR*w_k;
                        MF(num_strategy,1)=MF(num_strategy,1)+F^2*w_k;
                        MF(num_strategy,2)=MF(num_strategy,2)+F*w_k;
                    if max(0,result(end))>max(0,Y(end))
                        result=Y;
                    elseif max(0,result(end))==max(0,Y(end))&&result(end-1)>Y(end-1)
                        result=Y;
                    end
                end
%                 if FES==2000*D
%                     result_stage1=result;
%                 end
%                 if FES==10000*D
%                     result_stage2=result;
%                 end
            end
            
         %%%%%%%%%%%%%%%%%% Update the archives of CR and F %%%%%%%%%%%%%%%%%%%%%%%
           for hh=1:K
            if label_success(hh)>0
                flag_success(hh)=1;
                if H_length(hh)<H
                    H_length(hh)=H_length(hh)+1;
                end
                if H_label(hh)<H
                    H_label(hh)=H_label(hh)+1;
                    MCR_archive(H_label(hh),hh)= MCR(hh)/w_sum(hh);
                    MF_archive(H_label(hh),hh) = MF(hh,1)/MF(hh,2);

                else
                    H_label(hh)=1;
                    MCR_archive(H_label(hh),hh)= MCR(hh)/w_sum(hh);
                    MF_archive(H_label(hh),hh) = MF(hh,1)/MF(hh,2);     
                end
            end
           end
        %%%%%%%%%%%%%%%%%  Update epsilon  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fea_percent=sum(population(:,end)==0)/popsize; 
            epsilon=update_epsilon(epsilon,fea_percent,alpha,cp,Tc,FES/Max_FES, max_cv);%Update epsilon

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%% Display the progress of the algorithm %%%%%%%%%%%%%%%%%%%
%             if FES/Max_FES>=output_rate
%                output_rate=output_rate+0.1;
%                disp(strcat('The progress of the algorithm for Function ',num2str(I_fno),'_',num2str(D),'D on NO. ',num2str(run),' run', '£º  ',num2str(FES/Max_FES*100),'%'));
%             end
        %%%%%%%%%%%%%%% Reset the size of population %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            popsize=round(ini_popsize-(FES/Max_FES)*(ini_popsize-min_popsize));
            
            % Fix FES
            FES = FCounter;
            
%             % testing
%             [qf_y,qg_y,qh_y]=cec20_func(result(1:nvar),I_fno);
%             out = (sum([sum(qg_y.*(qg_y>0),1); sum(abs(qh_y).*(abs(qh_y)>1e-4),1)])./const_num);            
%             disp([num2str(FES) '   '   num2str(out) '     ' num2str(qf_y)]);

        %%%%%%%%%%%%%%%%%%%%%%%%%% A. Shirazi: FE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%
            Sol_arr(end+1,:) = result(1:nvar);
            FE_arr(end+1,1) = FCounter;    
            
        end
        
%         if run==1
%             datasize=size(result,2);
%             data=zeros(25*3,datasize);
%         end
%         data((run-1)*3+1,:)=result_stage1;
%         data((run-1)*3+2,:)=result_stage2;
%         data(run*3,:)=result;
%     end

end



