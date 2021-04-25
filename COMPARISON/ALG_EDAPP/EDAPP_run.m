function [Sol_arr,FE_arr] = EDAPP_run(prbNum)

    global initial_flag
    initial_flag = 0;    
    
    global fcnNum const_num;
    fcnNum = prbNum;
    
%%%%%%%%%%%%%%%%%%%%%%%%%% A. Shirazi: FE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%    
    global FCounter;
    FCounter = 0;
    Sol_arr = [];
    FE_arr = [];
    
    setup = Introd_Par(fcnNum);
        nvar  = setup.n;
        lb    = setup.xmin;
        ub    = setup.xmax;
        FEMax = setup.Max_FES;
        const_num = setup.gn(fcnNum) + setup.hn(fcnNum);

    npop = min(200,max(10*nvar,50));
    
    mappingtype = 'lin';
    mappingmethod = 'det';    

    objFcn = @objFcnCEC2020;
    conFcn = @conFcnCEC2020;
    
    disp('---------------------')
    disp('------ SEEDING ------')
    disp('---------------------')
    
    [x,c,FE] = seeding(conFcn,lb,ub,npop,FEMax); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% A. Shirazi: FE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%    
        Sol_arr(end+1,:) = x(1,:);
        FE_arr(end+1,1) = FCounter;  
    
    if nnz(c)~=0
        disp('-----------------------------------------------------------------------');
        disp('UNABLE TO FIND SUFFICIENT NUMBER OF FEASIBLE SOLUTIONS FOR OPTIMIZATION');
        disp('-----------------------------------------------------------------------');
        return;
    end

    f = evalFcn(objFcn,x);

    disp('--------------------------')
    disp('------ OPTIMIZATION ------')
    disp('--------------------------')
    
    iter = 0;
    while FE < FEMax
        iter=iter+1;
        
            [x_sel,f_sel]=selection(x,f);

            [parentclusters,smartclusters,curFE] = clustering(x_sel,f_sel,conFcn);    
            FE = FE + curFE;

            parentclusters = refinement(parentclusters,x_sel,f_sel);

            [parentclusters,smartclusters] = learning(parentclusters,smartclusters);
        
            [parentclusters,smartclusters] = sampling(parentclusters,smartclusters,npop);
        
            [parentclusters,smartclusters] = repairing(parentclusters,smartclusters,lb,ub);
        
            [parentclusters,smartclusters,~,overall_iter2correct,curFE] = mapping(parentclusters,smartclusters,conFcn,mappingtype,mappingmethod);
            FE = FE + overall_iter2correct + curFE;
        
            [parentclusters,smartclusters] = evaluation(parentclusters,smartclusters,objFcn);

            [x_new,f_new,x_best,f_best] = replacement(parentclusters,smartclusters,x,f);

            % Refine for cec2020 (FE limit)
            FE = FCounter;            
            % fprintf('Iter: %3.0f    Best Obj: %10.5g   FE/FEMax: %% %3.2f \n',[iter f_best min(100,100*FE/FEMax)])    
            x = x_new;
            f = f_new;
                    
            %%%%%%%%%%%%%%%%%%%%%%%%%% A. Shirazi: FE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%
                Sol_arr(end+1,:) = x_best;
                FE_arr(end+1,1) = FCounter;            
    end
   
end









