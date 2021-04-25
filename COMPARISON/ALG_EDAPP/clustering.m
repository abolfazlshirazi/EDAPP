function [parentclusters,smartclusters,curFE] = clustering(x_sel,f_sel,constFcn)
        nsel = size(x_sel,1);
        curFE = 0;
        
        flag = false;
        
        for jj=1:nsel
            [idx,mu] = kmeans(x_sel,jj);
            
            nc = nnz(evalFcn(constFcn,mu)<=0);
            curFE = curFE + jj;
            if nc==jj 
                flag = true;
                break;
            end
        end        
        if ~flag
            error('Err!')
        end
        parentclusters = cell(1,nc);

        for ii=1:nc
            parentclusters{ii}.x = x_sel(idx==ii,:);
            parentclusters{ii}.f = f_sel(idx==ii);
            parentclusters{ii}.c = mu(ii,:);
        end
        
        
        lambda = 0.2;
        psi = 1;
        smartclusters = cell(1,nc);
        counter = 0;
        for ii=1:nc
            curx = x_sel(idx==ii,:);
            curf = f_sel(idx==ii,:);
            curmu = mu(ii,:);
            cursigma = std(curx);
            curn = size(curx,1);

            [curf,curind] = sort(curf);
            curx = curx(curind,:);
            
            nbest = ceil(lambda*curn);
            bestx = curx(1:nbest,:);
            bestf = curf(1:nbest,:);

            outliersindx = (sqrt(sum((bestx-repmat(curmu,nbest,1)).^2,2)) > sqrt(sum((psi.*cursigma).^2,2)));
            if isempty(find(outliersindx,1))
                % No outliers
            else            
                counter = counter+1;
                
                cursmart = struct;
                cursmart.x = bestx(outliersindx,:);
                cursmart.f = bestf(outliersindx);
                cursmart.parent = curmu;
                
                smartclusters{counter} = cursmart;
            end
        end
        
        smartclusters = smartclusters(1:counter);
        
end