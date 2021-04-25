function [Sol_arr,FE_arr] = EDAPP_run(prbNum)

    global initial_flag
    initial_flag = 0;    
    
    global fcnNum const_num;
    fcnNum = prbNum;
    
%%%%%%%%%%%%%%%%%%%%%%%%%% A. Shirazi: FE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%    
    global FCounter;
    FCounter = 0;
    Sol_arr = [];
    FE_arr = [];
    
    setup = Introd_Par(fcnNum);
        nvar  = setup.n;
        lb    = setup.xmin;
        ub    = setup.xmax;
        FEMax = setup.Max_FES;
        const_num = setup.gn(fcnNum) + setup.hn(fcnNum);

    npop = min(200,max(10*nvar,50));
    
    mappingtype = 'lin';
    mappingmethod = 'det';    

    objFcn = @objFcnCEC2020;
    conFcn = @conFcnCEC2020;
    
    disp('---------------------')
    disp('------ SEEDING ------')
    disp('---------------------')
    
    [x,c,FE] = seeding(conFcn,lb,ub,npop,FEMax); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% A. Shirazi: FE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%    
        Sol_arr(end+1,:) = x(1,:);
        FE_arr(end+1,1) = FCounter;  
    
    if nnz(c)~=0
        disp('-----------------------------------------------------------------------');
        disp('UNABLE TO FIND SUFFICIENT NUMBER OF FEASIBLE SOLUTIONS FOR OPTIMIZATION');
        disp('-----------------------------------------------------------------------');
        return;
    end

    f = evalFcn(objFcn,x);

    disp('--------------------------')
    disp('------ OPTIMIZATION ------')
    disp('--------------------------')
    
    iter = 0;
    while FE < FEMax
        iter=iter+1;
        
            [x_sel,f_sel]=selection(x,f);

            [parentclusters,smartclusters,curFE] = clustering(x_sel,f_sel,conFcn);    
            FE = FE + curFE;

            parentclusters = refinement(parentclusters,x_sel,f_sel);

            [parentclusters,smartclusters] = learning(parentclusters,smartclusters);
        
            [parentclusters,smartclusters] = sampling(parentclusters,smartclusters,npop);
        
            [parentclusters,smartclusters] = repairing(parentclusters,smartclusters,lb,ub);
        
            [parentclusters,smartclusters,~,overall_iter2correct,curFE] = mapping(parentclusters,smartclusters,conFcn,mappingtype,mappingmethod);
            FE = FE + overall_iter2correct + curFE;
        
            [parentclusters,smartclusters] = evaluation(parentclusters,smartclusters,objFcn);

            [x_new,f_new,x_best,f_best] = replacement(parentclusters,smartclusters,x,f);

            % Refine for cec2020 (FE limit)
            FE = FCounter;            
            % fprintf('Iter: %3.0f    Best Obj: %10.5g   FE/FEMax: %% %3.2f \n',[iter f_best min(100,100*FE/FEMax)])    
            x = x_new;
            f = f_new;
                    
            %%%%%%%%%%%%%%%%%%%%%%%%%% A. Shirazi: FE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%
                Sol_arr(end+1,:) = x_best;
                FE_arr(end+1,1) = FCounter;            
    end
   
end

function varout=evalFcn(func,varin)
    n=size(varin,1); % rows are for each input, columns for variables        
    varout=zeros(n,1);    
    
    for ii=1:n
        curvarout = func(varin(ii,:));
        varout(ii) = curvarout;
    end
end


function [parent,smart] = learning(parent,smart)
    for ii=1:numel(parent)
        parent{ii}.sigma2 = var(parent{ii}.x);
    end
    for ii=1:numel(smart)
        [~,bestind] = min(smart{ii}.f);
        
        %centeroid
        smart{ii}.c = smart{ii}.x(bestind,:); % we choose best as the center, can be mean
        
        % variance
            % distance can be: "from mean to parent" or "from best to parent"
            % distance can be considered as variance or std        
        smart{ii}.sigma2 = 2.*abs(smart{ii}.c - smart{ii}.parent);
    end
end

function [parent,smart] = sampling(parent,smart,npop)
    percluster = floor(npop/(numel(parent)+numel(smart))); % just a quick way for now
    for ii=1:numel(parent)
        parent{ii}.x_sam = mvnrnd(parent{ii}.c,parent{ii}.sigma2,percluster);
    end
    for ii=1:numel(smart)
        smart{ii}.x_sam = mvnrnd(smart{ii}.c,smart{ii}.sigma2,percluster);
    end
end

function [parent,smart] = repairing(parent,smart,lb,ub)
    for ii=1:numel(parent)
        n = size(parent{ii}.x_sam,1);
        parent{ii}.x_sam = max(repmat(lb,n,1),parent{ii}.x_sam);
        parent{ii}.x_sam = min(repmat(ub,n,1),parent{ii}.x_sam);        
    end
    for ii=1:numel(smart)
        n = size(smart{ii}.x_sam,1);
        smart{ii}.x_sam = max(repmat(lb,n,1),smart{ii}.x_sam);
        smart{ii}.x_sam = min(repmat(ub,n,1),smart{ii}.x_sam);        
    end
end

function [parent,smart,num_of_correctedpoints,overall_iter2correct,curFE] = mapping(parent,smart,constFcn,mappingtype,mappingmethod)
    num_of_correctedpoints = 0;
    overall_iter2correct = 0;
    curFE = 0;

    for ii=1:numel(parent)

        x_sam = parent{ii}.x_sam;
        center = parent{ii}.c;
            % separate feasible and infeasible points
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
            % separate feasible and infeasible points
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
                    
                    
%                     disp([num2str(kk) ' iterations for mapping this point'])
                    curxinfeas(jj,:) = curpoint;
                end
                x_sam = [curxfeas;curxinfeas];
            end        
        
        smart{ii}.x_sam = x_sam;
    end
end

function [newp,iter] = mapmechanism(p,c,conFcn,type,method)

%     c,p, conFcn
%     type   = 'lin'; % lin , bis
%     method = 'sto'; % det sto
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

function [parent,smart] = evaluation(parent,smart,func)
    for ii=1:numel(parent)      
        parent{ii}.f_sam = evalFcn(func,parent{ii}.x_sam); 
    end
    for ii=1:numel(smart)
        smart{ii}.f_sam = evalFcn(func,smart{ii}.x_sam);
    end
end

function [x,c,FE] = seeding(conFcn,lb,ub,npop,FEMax)
    FE = 0; % function evaluation
    lambda = 0.5; % selection
    nsel = round(lambda*npop);
    
    lambda2 = 0.2; 
    ninc = max(1,round(lambda2*npop)); % include in maxtry iterations
    
    itermax = 100; % every itermax... it refines population to start from the beginning
    stallIter = 0;
    
    iter = 0;
    trymax = false;
    while FE <= FEMax - npop
        iter = iter + 1;
        
        if iter == 1
            
            x = unifrnd(repmat(lb,npop,1),repmat(ub,npop,1));
            c = evalFcn(conFcn,x);
            FE = FE + npop;            
    
            % sort
            [c,ind] = sort(c);
            x = x(ind,:);
            
        elseif trymax
            trymax = false; 
            stallIter = 0;
            
            x = [x(1:ninc,:) ; unifrnd(repmat(lb,npop-ninc,1),repmat(ub,npop-ninc,1)) ];
            c = evalFcn(conFcn,x);
            FE = FE + npop;            
    
            % sort
            [c,ind] = sort(c);
            x = x(ind,:);
            
        else
            
           % select        
            x_sel = x(1:nsel,:);

            % learning
            myCen = x_sel(1,:);
            myCov = cov(x_sel);

            % sampling
            x_new = mvnrnd(myCen,myCov,npop);

            % repairing
            x_new = max(x_new,repmat(lb,npop,1));
            x_new = min(x_new,repmat(ub,npop,1));

            % evaluation
            c_new = evalFcn(conFcn,x_new);
            FE = FE + npop;

            % stall check
            if min(c_new) >= c(1) % if no improvement
                stallIter = stallIter + 1;
            else
                stallIter = 0;                
            end
            
            % replacement
            x_new = [x_new;x];
            c_new = [c_new;c];
            [c_new,ind] = sort(c_new);
            x_new = x_new(ind,:);

            x = x_new(1:npop,:);
            c = c_new(1:npop,:);
            
        end        
        
        % report
        % fprintf('Seeding Iter: %6.0f   Violation: %10.5g   FE/FEMax: %% %5.2f   Pop Feas: %% %3.2f   Stall: %i \n',[iter c(1) min(100,100*FE/FEMax) 100*nnz(c<=0)/npop stallIter])
        
        % check for stop
        if nnz(c<=0)==npop
            disp('---------------------------------------------------------------')
            disp('    SUFFICIENT NUMBER OF FEASIBLE SOLUTIONS HAVE BEEN FOUND    ')
            disp('---------------------------------------------------------------')
            break;
        elseif nnz(c<=0)>0 % at least one feasible solution has been found: continue then
            continue;
        else
            if rem(iter,itermax)==0 % iterations are passed for traymax
                % disp('    RETRYING...')
                trymax = true;
            end
        end
    end
    
end

function [x_sel,f_sel]=selection(x,f)
    trnc_factor = 0.5;
    npop = numel(f);
    nsel = round(trnc_factor*npop);

    if nsel < 2
        nsel = 2;
    end
    if nsel >= npop
        nsel = npop-1;
    end
        
    [f_sorted,ind_sorted] = sort(f);
    x_sorted = x(ind_sorted,:);
    
    x_sel = x_sorted(1:nsel,:);
    f_sel = f_sorted(1:nsel);
    
end

function [parentclusters,smartclusters,curFE] = clustering(x_sel,f_sel,constFcn)
        nsel = size(x_sel,1);
        curFE = 0;
        
%         disp('Searching for best number of clusters, satisfying the constraint')        
        flag = false; % found best cluster
        
        for jj=1:nsel
            [idx,mu] = kmeans(x_sel,jj);
            
            nc = nnz(evalFcn(constFcn,mu)<=0); % number of good clusters
            curFE = curFE + jj;
            if nc==jj % if all clusters are good
                flag = true;
                break;
            end
        end        
        if ~flag
            error('mage mishe?')
        end
        parentclusters = cell(1,nc);
        % saving parentclusters
        for ii=1:nc
            parentclusters{ii}.x = x_sel(idx==ii,:);
            parentclusters{ii}.f = f_sel(idx==ii);
            parentclusters{ii}.c = mu(ii,:);
        end
        

        lambda = 0.2; % percentage of top good points to check in each cluster
        psi = 1; % outlier distance for check (1sigma, 2sigma, 3sigma)
        smartclusters = cell(1,nc);
        counter = 0;
        for ii=1:nc
            curx = x_sel(idx==ii,:);
            curf = f_sel(idx==ii,:);
            curmu = mu(ii,:);
            cursigma = std(curx);
            curn = size(curx,1);

            % sort this cluster
            [curf,curind] = sort(curf);
            curx = curx(curind,:);
            
            %extract 10 of good points to check
            nbest = ceil(lambda*curn); % if it is 0 it will be disabled, if it is 1, all points will be checked            
            bestx = curx(1:nbest,:);
            bestf = curf(1:nbest,:);

            outliersindx = (sqrt(sum((bestx-repmat(curmu,nbest,1)).^2,2)) > sqrt(sum((psi.*cursigma).^2,2)));
            if isempty(find(outliersindx,1))
                % no outliers, then no smartclusters
            else            
                counter = counter+1;
                % create a new smart cluster
                cursmart = struct;
                cursmart.x = bestx(outliersindx,:);
                cursmart.f = bestf(outliersindx);
                cursmart.parent = curmu;
                
                smartclusters{counter} = cursmart;
            end
        end
        
        % trash the empty at the tail
        smartclusters = smartclusters(1:counter);
        
        %%%%%%%%%%% Now combine two clusters
end


function [parentclusters,num_of_ref] = refinement(parentclusters,x_sel,f_sel)
    num_of_ref = 0;
    for ii=1:numel(parentclusters)
        if size(parentclusters{ii}.x,1)==1    
            num_of_ref = num_of_ref+1;

            % first, lets remove the current point from the search space
            tocheck = (x_sel == repmat(parentclusters{ii}.x,size(x_sel,1),1));

            thisind = 0;
            for kk=1:size(tocheck,1)
                if tocheck(kk,:) == ones(1,size(x_sel,2))
                    thisind = kk;
                    break;
                end
            end
            if thisind~=0
                x_rem = x_sel((1:size(x_sel,1))'~=thisind,:);
                f_rem = f_sel((1:size(x_sel,1))'~=thisind);
            else

                % this should never happen. the cluster has one point and it should be center
                % if it is happenning, it maybe due to kmeans performance
                % if happens, no need to remove from search space
                x_rem = x_sel; % previous command line also works anyways
                f_rem = f_sel; % previous command line also works anyways
            end

            d = x_rem - repmat(parentclusters{ii}.x,size(x_rem,1),1);
            [~,q2] = min(sum(d.^2,2));
            thepointtoaddtothiscluster = x_rem(q2,:);
            thepointtoaddtothisclusterf = f_rem(q2);
            % update cluster
            parentclusters{ii}.x = [parentclusters{ii}.x; thepointtoaddtothiscluster];
            parentclusters{ii}.f = [parentclusters{ii}.f; thepointtoaddtothisclusterf];
            
            %it does NOT update center point
            
        end
    end
end

function [x_new,f_new,x_best,f_best] = replacement(parentclusters,smartclusters,x,f)
    % Best Elitism    
    npop = numel(f);

        % gathering
        x_sam = [];
        f_sam = [];
        for jj=1:numel(parentclusters)
            x_sam = [x_sam;parentclusters{jj}.x_sam];
            f_sam = [f_sam;parentclusters{jj}.f_sam];
        end
        for jj=1:numel(smartclusters)
            x_sam = [x_sam;smartclusters{jj}.x_sam];
            f_sam = [f_sam;smartclusters{jj}.f_sam];
        end
        
    x = [x;x_sam];
    f = [f;f_sam];
    
    [f_sorted,ind]= sort(f);
    x_sorted = x(ind,:);
    
    x_new = x_sorted(1:npop,:);
    f_new = f_sorted(1:npop,:);   
    
    % just for best elitism (be careful)
    x_best = x_new(1,:);
    f_best = f_new(1);

end

function f = objFcnCEC2020(x)
    global fcnNum;    
    f = cec20_func(x,fcnNum);
end

function out = conFcnCEC2020(x)
    global fcnNum const_num;
    [~,g,h] = cec20_func(x,fcnNum);    
    out = (sum([sum(g.*(g>0),1); sum(abs(h).*(abs(h)>1e-4),1)])./const_num);
end

