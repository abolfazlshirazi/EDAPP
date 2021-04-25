function [p,objF,conV,archive,archiveobjF,archiveconV,FES]=EvolutionaryOperation(p,objF,conV,archive,archiveobjF,archiveconV,weights,FES,minVar,maxVar,problem,popsize,aaa)
% global problemSetNum
    
    [p,objF,conV]=Exchange(p,objF,conV,archive,archiveobjF,archiveconV);
                    
    % DE operation
    trial=DEgenerator(p,objF,conV,archive,weights,minVar,maxVar);
%     if problemSetNum==2010
%     	[objFtrial,conVtrial]=fitness2010(trial,problem);
%     else
%     	[objFtrial,conVtrial]=fitness2006(trial,problem,aaa);
%     end
  	[objFtrial,conVtrial]=fitness2020(trial,problem);

    FES=FES+popsize;
            
    % record the frontier population information for comparison
    objF1=objF;
    conV1=conV;
            
    % update evolution population
    [p,objF,conV]=environmentSelect(p,objF1,conV1,trial,objFtrial,conVtrial,weights);
            
    % update archive population
    [archive,archiveobjF,archiveconV]=Debselect(archive,archiveobjF,archiveconV,trial,objFtrial,conVtrial);

