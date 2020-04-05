function [x,fval,x_,fval_]=pso(CostFunction,ConstraintFunction,nVar,LBound,UBound,x_init,penaltyfactor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PARTICLE SWARM OPTIMIZATION %%%%%%%%%%%%%%%%%%%
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
% narginchk(2,11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct default options
%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultoptions.nPop=max(20*nVar,100);
defaultoptions.MaxIt=30*nVar;
defaultoptions.StopCr=-inf;
defaultoptions.Display='iter';
defaultoptions.AlgParams=[1 0.99 2 2];
defaultoptions.InitialPop=x_init;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Options Refinement %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
    InitialPop=defaultoptions.InitialPop;
    AlgParams=defaultoptions.AlgParams;
    Display=defaultoptions.Display;
    StopCr=defaultoptions.StopCr;
    MaxIt=defaultoptions.MaxIt;
    nPop=defaultoptions.nPop;


%%%%%%%%%%%%%%%%%%%%%%
%%% PSO Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%
w       =   AlgParams(1);     % Inertia Weight
wdamp   =   AlgParams(2);     % Inertia Weight Damping Ratio
c1      =   AlgParams(3);     % Personal Learning Coefficient
c2      =   AlgParams(4);     % Global Learning Coefficient

% Convert Position Limits to Matrices to avoid using repmat too much
LBound=repmat(LBound,[nPop 1]);
UBound=repmat(UBound,[nPop 1]);

% Velocity Limits
VelMax=0.1*(UBound-LBound);
VelMin=-VelMax;

%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization %%%
%%%%%%%%%%%%%%%%%%%%%%
FunEval=0;

% Initialize Position
particlePosition=myunifrnd(LBound,UBound);
if ~isempty(InitialPop)
    particlePosition(1:size(InitialPop,1),:)=InitialPop;
end

% Initialize Velocity
particleVelocity=zeros(nPop,nVar);
% Evaluation
particleCost=evalfunc(CostFunction,ConstraintFunction,particlePosition,penaltyfactor);

FunEval=FunEval+nPop;

% Update Personal Best
particleBestPosition=particlePosition;
particleBestCost=particleCost;
% Update Global Best
[GlobalBestCost,ind]=min(particleBestCost);
GlobalBestPosition=particleBestPosition(ind,:);

% Setup display header 
if strcmp(Display,'iter')
    fprintf('\n                                 Best             Stall\n');
    fprintf(  'Iteration     f-count            f(x)           Iterations\n');
    fprintf('%5.0f         %7.0f    %12.4g            %5.0f\n', ...
        0, FunEval, GlobalBestCost, 0);
end

%%%%%%%%%%%%%%%%%%%%%%
%%% PSO Main Loop %%%%
%%%%%%%%%%%%%%%%%%%%%%
AllStat=zeros(MaxIt,nVar+1);

for it=1:MaxIt
    % Update Velocity
    particleVelocity = w*particleVelocity ...
        +c1*rand([nPop nVar]).*(particleBestPosition-particlePosition) ...
        +c2*rand([nPop nVar]).*(repmat(GlobalBestPosition,[nPop 1])-particlePosition);

    % Apply Velocity Limits
    particleVelocity = max(particleVelocity,VelMin);
    particleVelocity = min(particleVelocity,VelMax);
    
    % Update Position
    particlePosition = particlePosition + particleVelocity;    
    
    % Velocity Mirror Effect
    IsOutside=(particlePosition<LBound) | (particlePosition>UBound);
    particleVelocity(IsOutside)=-particleVelocity(IsOutside);
    
    % Apply Position Limits
    particlePosition = max(particlePosition,LBound);
    particlePosition = min(particlePosition,UBound);
    
    % Evaluation    
    particleCost=evalfunc(CostFunction,ConstraintFunction,particlePosition,penaltyfactor);
    
    FunEval=FunEval+nPop;
    
    % Update Personal Best Cost
    ind=particleCost<particleBestCost;
    particleBestCost = (ind).*particleCost + (1-ind).*particleBestCost;
    
    % Update Personal Best Position
    repind=repmat(ind,[1 nVar]);
    particleBestPosition = repind.*particlePosition + (1-repind).*particleBestPosition;
    
    % Update Global Best
    [GlobalBestCost,ind]=min(particleBestCost);
    GlobalBestPosition=particleBestPosition(ind,:);    
    
    w=w*wdamp;
    %%%%%%%%%%%%%%%%%%%%%%
    %%%% Save Statistics %%%
    %%%%%%%%%%%%%%%%%%%%%%
    AllStat(it,:)=[GlobalBestCost GlobalBestPosition];    

    %%%%%%%%%%%%%%
    %%%% Display %%%
    %%%%%%%%%%%%%%
    if strcmp(Display,'iter')
        fprintf('%5.0f         %7.0f    %12.4g            %5.0f \n', ...
            it, FunEval, GlobalBestCost, 0);
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Stopping Criteria %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if GlobalBestCost <= StopCr
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

function r = myunifrnd(a,b)
    % Replacement for MATLAB unifrnd fnction
    a2 = a/2;
    b2 = b/2;
    mu = a2+b2;
    sig = b2-a2;

    r = mu + sig .* (2*rand(size(a))-1);

    % Fill in elements corresponding to illegal parameter values
    if ~isscalar(a) || ~isscalar(b)
        r(a > b) = NaN;
    elseif a > b
        r(:) = NaN;
    end
end

function varout=evalfunc(func,conFcn,varin,penaltyfactor)
    n=size(varin,1); % rows are for each input, columns for variables    
    varout=zeros(n,1);    
        for ii=1:n
            varout(ii)=penaltyFcn(func,conFcn,varin(ii,:),penaltyfactor);
        end
end

function out = penaltyFcn(objFcn,conFcn,x,k)
    obj = objFcn(x);
    con = conFcn(x);
    
    violation = sum(con.*(con>0));
    out = obj + k*violation;
end

