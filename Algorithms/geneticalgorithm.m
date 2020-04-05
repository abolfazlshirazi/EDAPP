function [x,fval,x_,fval_]=geneticalgorithm(CostFunction,ConstraintFunction,nVar,LBound,UBound,x_init,penaltyfactor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% GENETIC ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    By: Abolfazl Shirazi (ashirazi@bcamath.org)   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Date: 2020-04-02   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: 
%       This code is a modified version for evaluation of its
%       performance against MGD and it is not meant for other
%       purposes. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------
% ------ Begin ------
% ------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check Number of Inputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% narginchk(2,9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Construct default options %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defaultoptions.LBound=-inf.*ones(1,nVar);
% defaultoptions.UBound=+inf.*ones(1,nVar);
defaultoptions.nPop=max(20*nVar,100);
defaultoptions.MaxIt=30*nVar;
defaultoptions.StopCr=-inf;
defaultoptions.Display='iter';
defaultoptions.SelfAdjustment=2;
defaultoptions.SocialAdjustment=2;
defaultoptions.AlgParams=[0.7 0.4 0.3 0.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Options Refinement %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
    AlgParams=defaultoptions.AlgParams;
    Display=defaultoptions.Display;
    StopCr=defaultoptions.StopCr;
    MaxIt=defaultoptions.MaxIt;
    nPop=defaultoptions.nPop;
    
%%%%%%%%%%%%%%%%%%%%%
%%% GA Parameters %%%
%%%%%%%%%%%%%%%%%%%%%
pc      =   AlgParams(1);         % Crossover Percentage
gamma   =   AlgParams(2);         % Extra Range Factor for Crossover
pm      =   AlgParams(3);         % Mutation Percentage
mu      =   AlgParams(4);         % Mutation Rate

nm=round(pm*nPop);      % Number of Mutants
nc=2*round(pc*nPop/2);  % Number of Offsprings (also Parnets)

% pause(0.01); % Due to a bug in older versions of MATLAB

% Initialization
FunEval=0;
empty_individual.Position=[];
empty_individual.Cost=[];

pop=repmat(empty_individual,nPop,1);

for i=1:nPop    
    pop(i).Position=x_init(i,:);
    
    % Evaluation
    pop(i).Cost=penaltyFcn(CostFunction,ConstraintFunction,pop(i).Position,penaltyfactor);
    FunEval=FunEval+1;
end

% Sort Population
Costs=[pop.Cost];
[~, SortOrder]=sort(Costs);
pop=pop(SortOrder);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

% Store Cost
WorstCost=pop(end).Cost;

% Main Loop
AllStat=zeros(MaxIt,nVar+1);
for it=1:MaxIt
    
    % Crossover
    popc=repmat(empty_individual,nc/2,2);
    for k=1:nc/2
        
        % Select Parents Indices
            i1=randi([1 nPop]);
            i2=randi([1 nPop]);

        % Select Parents
        p1=pop(i1);
        p2=pop(i2);
        
        % Apply Crossover               
        [popc(k,1).Position, popc(k,2).Position]=ga_crossover(p1.Position,p2.Position,gamma,LBound,UBound);
        
        % Evaluate Offsprings      
        popc(k,1).Cost=penaltyFcn(CostFunction,ConstraintFunction,popc(k,1).Position,penaltyfactor);
        popc(k,2).Cost=penaltyFcn(CostFunction,ConstraintFunction,popc(k,2).Position,penaltyfactor);
        
        FunEval=FunEval+2;
        
    end
    popc=popc(:);
    
    % Mutation
    popm=repmat(empty_individual,nm,1);
    for k=1:nm
        
        % Select Parent
        i=randi([1 nPop]);
        p=pop(i);
        
        % Apply Mutation
        popm(k).Position=ga_mutate(p.Position,mu,LBound,UBound);       
        
        % Evaluate Mutant
        popm(k).Cost=penaltyFcn(CostFunction,ConstraintFunction,popm(k).Position,penaltyfactor);
        
        FunEval=FunEval+1;
        
    end
    % Create Merged Population
    pop=[pop
         popc
         popm]; %#ok
     
    % Sort Population
    Costs=[pop.Cost];
    [~, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
    
    % Update Worst Cost
    WorstCost=max(WorstCost,pop(end).Cost);
    
    % Truncation
    pop=pop(1:nPop);
    
    % Store Best Solution Ever Found
    BestSol=pop(1);
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Save Statistics %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    AllStat(it,:)=[BestCost(it) BestSol.Position];    

    %%%%%%%%%%%%%%%
    %%% Display %%%
    %%%%%%%%%%%%%%%
    if strcmp(Display,'iter')
        fprintf('%5.0f         %7.0f    %12.4g            %5.0f\n', ...
            it, FunEval, BestSol.Cost, 0);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Stopping Criteria %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if BestCost(it) <= StopCr
        AllStat=AllStat(1:it,:);
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Construct the outputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x =         AllStat(end,2:end);
fval =      AllStat(end,1);
x_ =        AllStat(:,2:end);
fval_ =     AllStat(:,1);

end

function [y1, y2]=ga_crossover(x1,x2,gamma,LBound,UBound)
    alpha=unifrnd(-gamma,1+gamma,size(x1));
    
    y1=alpha.*x1+(1-alpha).*x2;
    y2=alpha.*x2+(1-alpha).*x1;
    
    y1=max(y1,LBound);
    y1=min(y1,UBound);
    
    y2=max(y2,LBound);
    y2=min(y2,UBound);
end

function y=ga_mutate(x,mu,LBound,UBound)
    nVar=numel(x);
    nmu=ceil(mu*nVar);
    
    j=randsample(nVar,nmu)';
    sigma=0.1*(UBound-LBound);
    
    y=x;
    
    y(j)=x(j)+sigma(j).*randn(size(j));
    
    y=max(y,LBound);
    y=min(y,UBound);
end

function out = penaltyFcn(objFcn,conFcn,x,k)
    obj = objFcn(x);
    con = conFcn(x);
    
    violation = sum(con.*(con>0));
    out = obj + k*violation;
end

