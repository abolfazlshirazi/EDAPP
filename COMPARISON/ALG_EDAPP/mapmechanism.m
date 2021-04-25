function [newp,iter] = mapmechanism(p,c,conFcn,type,method)

    itermax = 50;
    
        
    curp = p;    
    iter = 0;
    isfound = 0;
    for ii=itermax:-1:1
        if strcmp(type,'lin') 
            delta = (c-curp)./ii;
        elseif strcmp(type,'bis') 
            delta = (c-curp)./2;
        end

        if strcmp(method,'det')
            r=ones(1,size(p,2));
        elseif strcmp(method,'sto')
            r=rand(1,size(p,2));
        end        
        curp = curp + delta.*(r);

        iter = iter+1;
        
        CC = evalFcn(conFcn,curp);
        
        if isempty(find(CC>0,1))
            isfound = true;
            break;
        end
    end

    if ~isfound
        curp = c;
    end
    
    newp = curp;    
end