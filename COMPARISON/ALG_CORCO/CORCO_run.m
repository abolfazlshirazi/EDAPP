function [Sol_arr,FE_arr] = CORCO_run(prbNum)

%***********************************************************************************
% Author: Y Wang, JP Li, XH Xue, and BC Wang
% Last Edited: 112/03/2019
% Eamil: bingcwang3-c@my.cityu.edu.hk; ywang@csu.edu.cn; ljpcsu@csu.edu.cn;
% Reference: Utilizing the Correlation between Constraints and Objective Function 
% for Constrained Evolutionary Optimization, IEEE Transactions on
% Evolutionary Computation, in press.
%***********************************************************************************

% clc; 
% clear all; 
% format long g;
% tic;

global gen maxGen CorIndex learningGenNum DivIndex ExFlag diversityDerta feasibleRatioInitial

problemSet = prbNum;

    global initial_flag
    initial_flag = 0;    
    
    global const_num;
    
    global FCounter;
    FCounter = 0;
    Sol_arr = [];
    FE_arr = [];    
    
    setup = Introd_Par(problemSet);
    totalFES = setup.Max_FES;
    const_num = setup.gn(problemSet) + setup.hn(problemSet);
    Par = Cal_par(problemSet);        
    minVar = Par.xmin;
    maxVar = Par.xmax;
    aaa=0;
    n = Par.n;  
    popsize = 100; % default for CORCO
    



% if problemSetNum==2010
%     % problem set of 18 benchmark functions to be tested
%     problemSet=[1:18];
%     
%     % the number of dimensions which can be set as 10 or 30
%     n=30;
%     
%     % the size of evolving population
%     if n==10
%         popsize=60;
%     elseif n==30
%         popsize=100;
%     end
%     
%     % the maximum number of function evaluations
%     if n==10
%         totalFES=200000;
%     elseif n==30
%         totalFES=600000; 
%     end
% else
%     problemSet=1;
%     popsize=100;
%     totalFES=500000; 
% end

maxGen=ceil(totalFES/popsize);

% choose the function to be tested
for problem = problemSet
    % fprintf('problem:%d\n',problem);
    % the value used to record the number of runs
    time=1;
    
    % the value which represents the number of total runs
    totalTime=1;
    
    % the array used to record the best solution of each run
    minValue=[];
    
    while time <= totalTime
%         if problemSetNum==2010
%             [minVar,maxVar] = problemSelection2010(problem,n);aaa=0;
%         else
%             [minVar,maxVar,n,aaa] = problemSelection2006(problem);
%         end

%         rand('seed',sum(100*clock));
        
        % initialization
        p=repmat(minVar,popsize,1)+rand(popsize,n).*repmat((maxVar-minVar),popsize,1);
        
%         if problemSetNum==2010
%             [objF, conV]=fitness2010(p, problem); 
%         else
%             [objF, conV]=fitness2006(p, problem,aaa);
%         end

        % A. Shirazi        
        [objF, conV]=fitness2020(p, problem); 
        
        FES=popsize;
       
        nondominateIndex=1:popsize;
        archive=p(nondominateIndex,:); archiveobjF=objF(nondominateIndex); archiveconV=conV(nondominateIndex);
       

        gen=0; X=0; ExFlag=0; CorIndex=0; diversityDerta = 0;
        betterRecord=[];learningGenNum=maxGen/20;
        betterRecord1 = [];betterRecord2 = [];
        diversityInitial=sum(std(p)./(maxVar-minVar))/n;
        feasibleNum = size(find(conV==0),1);
        feasibleRatioInitial=feasibleNum/popsize;
        
        % learning stage
        while gen<=learningGenNum         
            
            % weights generation and assignment
            weights = WeightGenerator(X,popsize,conV,objF);
            
            % Evolutionary operation
            [p,objF,conV,archive,archiveobjF,archiveconV,FES]=EvolutionaryOperation(p,objF,conV,archive,archiveobjF,archiveconV,weights,FES,minVar,maxVar,problem,popsize,aaa);
            
            % record the better number
            [con_obj_betterNum,obj_con_betterNum]=InterCompare(archiveobjF,archiveconV,objF,conV);
            
            % Calculate the DivIndex
            DivIndex=sum(std(p)./(maxVar-minVar))/n;
            
            % Record
            betterRecord1 = [betterRecord1;con_obj_betterNum];
            betterRecord2 = [betterRecord2;obj_con_betterNum];
            
            gen=gen+1;
        end
        
       
        % 计算0占比
      %  indicator1 = find(betterRecord ~= 0);
      %  betterLength = length(indicator1);
      %  CorIndex=betterLength/learningGenNum ;
        
        indicator1 = find(betterRecord1 ~= 0);
        indicator2 = find(betterRecord2 ~= 0);

        betterLength1 = length(indicator1);        
        betterLength2 = length(indicator2);
        betterLength = min(betterLength1,betterLength2);
        CorIndex=betterLength/learningGenNum ;
        
        diversityDerta=diversityInitial - DivIndex;
        Diver = DivIndex;
            
        while FES<totalFES
        
            % 产生并分配权重
            weights = WeightGenerator(X,popsize,conV,objF);
            
            % Evolutionary operation
            [p,objF,conV,archive,archiveobjF,archiveconV,FES]=EvolutionaryOperation(p,objF,conV,archive,archiveobjF,archiveconV,weights,FES,minVar,maxVar,problem,popsize,aaa);

            X=X+1/maxGen; gen=gen+1;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % A. Shirazi: Save best X + FCounter
            feasiIndex=find(archiveconV==0);
            if ~isempty(feasiIndex)
                % disp('feasible solution is found')            
                myF = archiveobjF(feasiIndex);
                myX = archive(feasiIndex,:);

                [~,myInd] = sort(myF);
                bestX = myX(myInd(1),:);
            else
                % disp('No feasible solutions');
                [~,myInd] = sort(archiveconV);
                bestX = archive(myInd(1),:);
            end
            Sol_arr(end+1,:) = bestX;
            FE_arr(end+1,1) = FCounter;             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            
            
        end
        
        
%         % calculate the index of feasible solutions
%         feasiIndex=find(archiveconV==0);
%         
%         if ~isempty(feasiIndex)
%             
%             % sort the feasible solutions according to their objective
%             % function values
%             [sortedFeasiVec,~]=sort(archiveobjF(feasiIndex));
%             
%             % the best solution is the feasible solution with minimum
%             % objective function value
%             if ExFlag==1
%                 Ex_bestSolution=sortedFeasiVec(1)
%                 % record the best solution of each run
%                 minValue=[minValue Ex_bestSolution];
%             else
%                 NonEx_bestSolution=sortedFeasiVec(1)
%                 % record the best solution of each run
%                 minValue=[minValue NonEx_bestSolution];
%             end
%             
%         else
%             
%             % if the population is infeasible, the best solution is
%             % assigned the value of NaN
%             [~,minconvIndex]=min(archiveconV);
%             bestSolution=archiveobjF(minconvIndex);
%             fprintf('*bestSolution=%d ExFlag=%d\n',bestSolution,ExFlag);            
%         end
        
        % modify the value of "time"
        time=time+1;
    end

    % calculate the best value, mean value, worst value and standard error
    % of the array of minValue
%     statisticValue=[min(minValue) mean(minValue) max(minValue) std(minValue)]
end
% run_time = toc

end