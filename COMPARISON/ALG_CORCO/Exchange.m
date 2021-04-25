function [p,objF,conV]=Exchange(p,objF,conV,archive,archiveobjF,archiveconV)
global gen learningGenNum DivIndex  feasibleRatioInitial CorIndex
if (gen>learningGenNum && rand<(1-CorIndex))

    if (gen>learningGenNum)
    	[~,index]=max(conV);
    	if archiveconV(index,:)<conV(index,:) || (archiveconV(index,:)==conV(index,:)&&archiveobjF(index,:)<objF(index,:))
        	p(index,:)=archive(index,:);
        	objF(index,:)=archiveobjF(index,:);
        	conV(index,:)=archiveconV(index,:);
        end
    end

    feasibleNum = size(find(conV==0),1);
    feasibleRatio=feasibleNum/size(p,1);
    
    if (feasibleRatio> 0.5 || feasibleRatioInitial>0)
    	p=archive;
    	objF=archiveobjF;
    	conV=archiveconV;
    end
    
    if (feasibleRatio==0 && rand<DivIndex)
        p=archive;
    	objF=archiveobjF;
    	conV=archiveconV;
    end
end