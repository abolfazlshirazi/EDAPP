function  [p,objF,conV,recordIndex,betterIndex1,betterIndex2]=Debselect(p,objF,conV,trial,objFtrial,conVtrial)

%%
global gen learningGenNum 
% temporary value to record the objective function value and degree of
% constraint violation
   termconV=conV; 
   termobjF=objF;
   [popsize,n]=size(p);
   
%  selection operator according to feasibility rule
  
  % feasible solutions are better than infeasible solutions and infeasible
  % solution with less degree of constraint violation is better than the one 
  % with larger degree of constraint violation 
   betterIndex1=find(conVtrial < termconV);
   p(betterIndex1,:)=trial(betterIndex1,:);
   objF(betterIndex1)=objFtrial(betterIndex1);
   conV(betterIndex1)=conVtrial(betterIndex1);

   if gen>learningGenNum 
    % feasible solution with less objective function value is better than the
    % one with larger objective function value
    betterIndex2=find(conVtrial==termconV & objFtrial < termobjF);      
    p(betterIndex2,:)=trial(betterIndex2,:);
    objF(betterIndex2)=objFtrial(betterIndex2);
    conV(betterIndex2)=conVtrial(betterIndex2);
   end
