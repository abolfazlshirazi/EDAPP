function [x,c,FE] = seeding(conFcn,lb,ub,npop,FEMax)
    FE = 0; 
    lambda = 0.5; 
    nsel = round(lambda*npop);
    
    lambda2 = 0.2; 
    ninc = max(1,round(lambda2*npop)); 
    
    itermax = 100; 
    stallIter = 0;
    
    iter = 0;
    trymax = false;
    while FE <= FEMax - npop
        iter = iter + 1;
        
        if iter == 1
            x = unifrnd(repmat(lb,npop,1),repmat(ub,npop,1));
            c = evalFcn(conFcn,x);
            FE = FE + npop;            
    
            [c,ind] = sort(c);
            x = x(ind,:);
            
        elseif trymax
            trymax = false; 
            stallIter = 0;
            
            
            x = [x(1:ninc,:) ; unifrnd(repmat(lb,npop-ninc,1),repmat(ub,npop-ninc,1)) ];
            c = evalFcn(conFcn,x);
            FE = FE + npop;            
    
            [c,ind] = sort(c);
            x = x(ind,:);
            
        else
            
            x_sel = x(1:nsel,:);

            myCen = x_sel(1,:);
            myCov = cov(x_sel);

            x_new = mvnrnd(myCen,myCov,npop);

            x_new = max(x_new,repmat(lb,npop,1));
            x_new = min(x_new,repmat(ub,npop,1));

            c_new = evalFcn(conFcn,x_new);
            FE = FE + npop;

            if min(c_new) >= c(1)
                stallIter = stallIter + 1;
            else
                stallIter = 0;                
            end
            
            x_new = [x_new;x];
            c_new = [c_new;c];
            [c_new,ind] = sort(c_new);
            x_new = x_new(ind,:);

            x = x_new(1:npop,:);
            c = c_new(1:npop,:);
            
        end
        
        if nnz(c<=0)==npop
            disp('---------------------------------------------------------------')
            disp('    SUFFICIENT NUMBER OF FEASIBLE SOLUTIONS HAVE BEEN FOUND    ')
            disp('---------------------------------------------------------------')
            break;
        elseif nnz(c<=0)>0 
            continue;
        else
            if rem(iter,itermax)==0 
                trymax = true;
            end
        end
    end
    
end