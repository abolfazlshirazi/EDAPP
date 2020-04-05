function [x,fval,x_,fval_,flag]=cmaes(CostFunction,ConstraintFunction,nVar,LBound,UBound,x_init)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% COVARIANCE MATRIX ADAPTATION EVOLUTIONARY STRATEGY %%%%%%%
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

% Population + MaxIt
nPop=20*nVar;
MaxIt=30*nVar;

% Parameters Settings
VarSize=[1 nVar]; 
lambda=nPop;
mu=round(lambda/2);
w=log(mu+0.5)-log(1:mu);
w=w/sum(w);
mu_eff=1/sum(w.^2);

% Control Step
sigma0=0.001*mean((UBound-LBound));
cs=(mu_eff+2)/(nVar+mu_eff+5);
ds=1+cs+2*max(sqrt((mu_eff-1)/(nVar+1))-1,0);
ENN=sqrt(nVar)*(1-1/(4*nVar)+1/(21*nVar^2));

% Parameters Update 
cc=(4+mu_eff/nVar)/(4+nVar+2*mu_eff/nVar);
c1=2/((nVar+1.3)^2+mu_eff);
alpha_mu=2;
cmu=min(1-c1,alpha_mu*(mu_eff-2+1/mu_eff)/((nVar+2)^2+alpha_mu*mu_eff/2));
hth=(1.4+2/(nVar+1))*ENN;

% Initialization
ps=cell(MaxIt,1);
pc=cell(MaxIt,1);
C=cell(MaxIt,1);
sigma=cell(MaxIt,1);

ps{1}=zeros(VarSize);
pc{1}=zeros(VarSize);
C{1}=eye(nVar);
sigma{1}=sigma0;

empty_individual.Position=[];
empty_individual.Step=[];
empty_individual.Cost=[];

M=repmat(empty_individual,MaxIt,1);

    M(1).Position=x_init(1,:); % only one
    M(1).Step=zeros(VarSize);
    M(1).Cost=CostFunction(M(1).Position);
    BestSol=M(1);
    
BestCost=zeros(MaxIt,1);

% Main Loop
AllStat=zeros(MaxIt,nVar+1);

for ii=1:MaxIt
    pop=repmat(empty_individual,lambda,1);      
    
    k = 10; % iteration for resampling
    thesteps = mvnrnd(zeros(VarSize),C{ii},k*lambda);
    thepos = repmat(M(ii).Position,k*lambda,1) + thesteps;
        % Repairing
        thepos = max(thepos,repmat(LBound,k*lambda,1));
        thepos = min(thepos,repmat(UBound,k*lambda,1));
        
        
    num_of_feasible = 0;
    for i=1:size(thepos,1)
        curstep = thesteps(i,:);
        curpos = thepos(i,:);
        thisconst = ConstraintFunction(curpos);
        isfeasible = isempty(find(thisconst > 0,1));
        
        if isfeasible
            num_of_feasible = num_of_feasible+1;
            pop(num_of_feasible).Step=curstep;
            pop(num_of_feasible).Position=curpos;
            pop(num_of_feasible).Cost=CostFunction(pop(num_of_feasible).Position);
            % Update Best Solution Ever Found
            if pop(num_of_feasible).Cost<BestSol.Cost
                BestSol=pop(num_of_feasible);
            end
        end
        
        if num_of_feasible>= lambda
            % disp('Sufficient feasible solutions are sampled.')            
            break;
        end
    end
    
    if num_of_feasible<lambda
        % disp('Insufficient fesible solutions')
        % Use infisible for the remaining
        
        c = num_of_feasible;
        
        thesteps = mvnrnd(zeros(VarSize),C{ii},lambda-num_of_feasible);
        thepos = repmat(M(ii).Position,lambda-num_of_feasible,1) + thesteps;
            % Repairing
            thepos = max(thepos,repmat(LBound,size(thepos,1),1));
            thepos = min(thepos,repmat(UBound,size(thepos,1),1));
        
        for jj=1:size(thepos,1)
            c = c+1;

            curstep = thesteps(jj,:);
            curpos = thepos(jj,:);
            
            pop(c).Step=curstep;
            pop(c).Position=curpos;
            pop(c).Cost=CostFunction(pop(c).Position);
            % we dont compare the objective with the best solution found so far because it is infeasible
        end
    end
    
    % Sort Population
    Costs=[pop.Cost];
    [~, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
    
    % Save Results
    BestCost(ii)=BestSol.Cost;
    
    % A quick check the center point
    thisconst = ConstraintFunction(M(ii).Position);
    ismeanfeasible = isempty(find(thisconst > 0,1));    
    
    % Display Results
    disp(['Iteration ' num2str(ii) ': Best Cost = ' num2str(BestCost(ii),20) '; P:' num2str(num_of_feasible/lambda,10), ' Is mean feasible:' num2str(ismeanfeasible)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Save Statistics %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    AllStat(ii,:)=[BestCost(ii) BestSol.Position];
    
    % Update Mean
    M(ii+1).Step=0;
    for j=1:mu
        M(ii+1).Step=M(ii+1).Step+w(j)*pop(j).Step;
    end    
    
    M(ii+1).Position=M(ii).Position+sigma{ii}*M(ii+1).Step;
    
    % Repair mean (added), still does not satisfy the constraint
        M(ii+1).Position = max(M(ii+1).Position,LBound);
        M(ii+1).Position = min(M(ii+1).Position,UBound);
    
    M(ii+1).Cost=CostFunction(M(ii+1).Position);
    
    % Update Step Size
    ps{ii+1}=(1-cs)*ps{ii}+sqrt(cs*(2-cs)*mu_eff)*M(ii+1).Step/chol(C{ii})';
    sigma{ii+1}=sigma{ii}*exp(cs/ds*(norm(ps{ii+1})/ENN-1))^0.3;
    
    % Update Covariance Matrix
    if norm(ps{ii+1})/sqrt(1-(1-cs)^(2*(ii+1)))<hth
        hs=1;
    else
        hs=0;
    end
    delta=(1-hs)*cc*(2-cc);
    pc{ii+1}=(1-cc)*pc{ii}+hs*sqrt(cc*(2-cc)*mu_eff)*M(ii+1).Step;
    C{ii+1}=(1-c1-cmu)*C{ii}+c1*(pc{ii+1}'*pc{ii+1}+delta*C{ii});
    for j=1:mu
        C{ii+1}=C{ii+1}+cmu*w(j)*pop(j).Step'*pop(j).Step;
    end
    
    % If Covariance Matrix is not Positive Defenite or Near Singular
    try
        chol(C{ii+1});
    catch flag
        C{ii+1} = C{ii};
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
