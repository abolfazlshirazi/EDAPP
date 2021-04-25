function [parent,smart,num_of_correctedpoints,overall_iter2correct,curFE] = mapping(parent,smart,constFcn,mappingtype,mappingmethod)
    num_of_correctedpoints = 0;
    overall_iter2correct = 0;
    curFE = 0;

    for ii=1:numel(parent)

        x_sam = parent{ii}.x_sam;
        center = parent{ii}.c;
            C = evalFcn(constFcn,x_sam);
                curFE = curFE + size(x_sam,1);
            curxfeas = x_sam(C<=0,:);
            curxinfeas = x_sam(C>0,:);
            
            num_of_correctedpoints = num_of_correctedpoints + size(curxinfeas,1);
            
            if size(curxfeas,1)==size(x_sam,1)
                x_sam = curxfeas;
            else
                for jj=1:size(curxinfeas,1)
                    curpoint = curxinfeas(jj,:);
                    
                    
                    [curpoint,this_iter2correct] = mapmechanism(curpoint,center,constFcn,mappingtype,mappingmethod);
                    overall_iter2correct = overall_iter2correct + this_iter2correct;
                    
                    curxinfeas(jj,:) = curpoint;
                end
                x_sam = [curxfeas;curxinfeas];
            end        
        
        parent{ii}.x_sam = x_sam;
    end
     
    for ii=1:numel(smart)

        x_sam = smart{ii}.x_sam;
        center = smart{ii}.c;
            C = evalFcn(constFcn,x_sam);
                curFE = curFE + size(x_sam,1);            
            curxfeas = x_sam(C<=0,:);
            curxinfeas = x_sam(C>0,:);
                        
            num_of_correctedpoints = num_of_correctedpoints + size(curxinfeas,1);
            
            if size(curxfeas,1)==size(x_sam,1)
                x_sam = curxfeas;
            else
                for jj=1:size(curxinfeas,1)
                    curpoint = curxinfeas(jj,:);
                    
                    
                    [curpoint,this_iter2correct] = mapmechanism(curpoint,center,constFcn,mappingtype,mappingmethod);
                    overall_iter2correct = overall_iter2correct + this_iter2correct;
                    
                    
                    curxinfeas(jj,:) = curpoint;
                end
                x_sam = [curxfeas;curxinfeas];
            end        
        
        smart{ii}.x_sam = x_sam;
    end
end