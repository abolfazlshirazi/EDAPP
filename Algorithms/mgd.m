function [x_best,f_best,x_,f_,ex_time] = mgd(objFcn,constFcn,nvar,lb,ub,x0,mappingtype,mappingmethod)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% CONSTRAINED ESTIMATION OF DISTRIBUTION ALGORITHM %%%%%%%%
%%%%%%%%%%%   MGD (Mixture of Gaussian Distribution)   %%%%%%%%%%%
%%%%%%%%    By: Abolfazl Shirazi (ashirazi@bcamath.org)   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Date: 2020-04-02   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This algorithm is a special kind of EDA that can be used to solve
% continuous optimization problems with nonlinear constraints. The
% algorithm works based on probabilistic models regarding the mixture
% of multivariate Gaussian distributions. It benefits from an intelligent
% learning process, a mapping method and an outlier detection technique.
% Feasibility of the solutions is guaranteed in this algorithm.

% IMPORTANT NOTE: 
%       THE CURRENT CODE IS THE INITIAL (UNCLEAN) VERSION OF THE
%       ALGORITHM. THE CODE IS STILL UNDER DEVELOPMENT. THE PURPOSE
%       OF PROVIDING THIS CODE IS FOR EVALUATION ONLY. PLEASE MODIFY
%       THE CODE ONLY IF YOU REALLY KNOW WHAT YOU ARE DOING.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------
% ------ Begin ------
% ------------------

    npop = 20*nvar; % must npop>=10
    ngen = 30*nvar;

    x = x0;  
    f = evalFcn(objFcn,x,'objective');

% ------------------
% ------ Loop ------
% ------------------
    x_ = zeros(ngen,nvar);
    f_ = zeros(ngen,1);
    
    thetimes = zeros(9,1);    
    for ii=1:ngen % Generations
    
        timeval = tic;
            [x_sel,f_sel]=selection(x,f);
        curdt = toc(timeval);
        thetimes(1) = thetimes(1)+curdt;
        
        timeval = tic;
            [parentclusters,smartclusters] = clustering(x_sel,f_sel,constFcn);    
        curdt = toc(timeval);
        thetimes(2) = thetimes(2)+curdt;
        
        timeval = tic;
            parentclusters = refinement(parentclusters,x_sel,f_sel);
        curdt = toc(timeval);
        thetimes(3) = thetimes(3)+curdt;
        
        timeval = tic;
            [parentclusters,smartclusters] = learning(parentclusters,smartclusters);
        curdt = toc(timeval);
        thetimes(4) = thetimes(4)+curdt;
        
        timeval = tic;
            [parentclusters,smartclusters] = sampling(parentclusters,smartclusters,npop);
        curdt = toc(timeval);
        thetimes(5) = thetimes(5)+curdt;
        
        timeval = tic;
            [parentclusters,smartclusters] = repairing(parentclusters,smartclusters,lb,ub);
        curdt = toc(timeval);
        thetimes(6) = thetimes(6)+curdt;
        
        timeval = tic;
            [parentclusters,smartclusters,~,~,omittime] = mapping(parentclusters,smartclusters,constFcn,mappingtype,mappingmethod);
                      
        curdt = toc(timeval)-omittime;
        thetimes(7) = thetimes(7)+curdt;
                
        timeval = tic;
            [parentclusters,smartclusters] = evaluation(parentclusters,smartclusters,objFcn);
        curdt = toc(timeval);
        thetimes(8) = thetimes(8)+curdt;

        timeval = tic;
            [x_new,f_new,x_best,f_best] = replacement(parentclusters,smartclusters,x,f);
        curdt = toc(timeval);
        thetimes(9) = thetimes(9)+curdt;
    
        x = x_new;
        f = f_new;
        
        x_(ii,:) = x_best;
        f_(ii,1) = f_best;
        
    end
    
    ex_time = sum(thetimes);
    
end

function varout=evalFcn(func,varin,fcnType,ExtraParams)
    n=size(varin,1); % Rows are for each input, columns for variables    
    varout=zeros(n,1);    
    if nargin == 3 || isempty(ExtraParams)
        for ii=1:n
            curvarout = func(varin(ii,:));
            if strcmp(fcnType,'objective')
                varout(ii) = curvarout;
            elseif strcmp(fcnType,'constraint')
                violations = sum(curvarout(curvarout>0)); % Sum of all constraint violations (c>0)
                if violations==0 % If all constraints are satisfied, return the sum of their negative values (c<=0)
                    violations = sum(curvarout);
                end
                varout(ii)=violations;
            else
                error('fcnType not recognized.')
            end
        end
    else
        for ii=1:n
            curvarout = func(varin(ii,:),ExtraParams);
            if strcmp(fcnType,'objective')
                varout(ii) = curvarout;
            elseif strcmp(fcnType,'constraint')
                violations = sum(curvarout(curvarout>0)); % Sum of all constraint violations (c>0)
                if violations==0 % If all constraints are satisfied, return the sum of their negative values (c<=0)
                    violations = sum(curvarout);
                end
                varout(ii)=violations;
            else
                error('fcnType not recognized.')
            end
        end
    end    
end

function [parent,smart] = learning(parent,smart)
    for ii=1:numel(parent)
        parent{ii}.sigma2 = var(parent{ii}.x);
    end
    for ii=1:numel(smart)
        [~,bestind] = min(smart{ii}.f);
        
        % Centeroid
        smart{ii}.c = smart{ii}.x(bestind,:); % We choose best as the center, can be mean
        
        % Variance
            % distance can be: "from mean to parent" or "from best to parent"
            % distance can be considered as variance or std        
        smart{ii}.sigma2 = 2.*abs(smart{ii}.c - smart{ii}.parent);
    end
end

function [parent,smart] = sampling(parent,smart,npop)
    percluster = floor(npop/(numel(parent)+numel(smart))); % Just a quick way for now
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

function [parent,smart,num_of_correctedpoints,overall_iter2correct,omittime] = mapping(parent,smart,constFcn,mappingtype,mappingmethod)
    num_of_correctedpoints = 0;
    overall_iter2correct = 0;
    omittime = 0;

    for ii=1:numel(parent)

        x_sam = parent{ii}.x_sam;
        center = parent{ii}.c;
            % Separate feasible and infeasible points
            C = evalFcn(constFcn,x_sam,'constraint');
            curxfeas = x_sam(C<=0,:);
            curxinfeas = x_sam(C>0,:);
            
            num_of_correctedpoints = num_of_correctedpoints + size(curxinfeas,1);
            
            if size(curxfeas,1)==size(x_sam,1)
                x_sam = curxfeas;
            else
                for jj=1:size(curxinfeas,1)
                    curpoint = curxinfeas(jj,:);
                    [curpoint,this_iter2correct,this_omittime] = mapmechanism(curpoint,center,constFcn,mappingtype,mappingmethod);
                    omittime = omittime + this_omittime;
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
            % Separate feasible and infeasible points
            C = evalFcn(constFcn,x_sam,'constraint');
            curxfeas = x_sam(C<=0,:);
            curxinfeas = x_sam(C>0,:);
                        
            num_of_correctedpoints = num_of_correctedpoints + size(curxinfeas,1);
            
            if size(curxfeas,1)==size(x_sam,1)
                x_sam = curxfeas;
            else
                for jj=1:size(curxinfeas,1)
                    curpoint = curxinfeas(jj,:);
                    [curpoint,this_iter2correct,this_omittime] = mapmechanism(curpoint,center,constFcn,mappingtype,mappingmethod);
                    omittime = omittime + this_omittime;
                    overall_iter2correct = overall_iter2correct + this_iter2correct;
                    curxinfeas(jj,:) = curpoint;
                end
                x_sam = [curxfeas;curxinfeas];
            end        
        smart{ii}.x_sam = x_sam;
    end
end

function [newp,iter,omittime] = mapmechanism(p,c,conFcn,type,method)
    itermax = 10;    
    omittime = 0;    
        
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
        
        timeval = tic;
            CC = evalFcn(conFcn,curp,'constraint');
        omittime = omittime + toc(timeval);
        
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

function [parent,smart] = evaluation(parent,smart,objFcn)
    for ii=1:numel(parent)
        parent{ii}.f_sam = evalFcn(objFcn,parent{ii}.x_sam,'objective');
    end
    for ii=1:numel(smart)
        smart{ii}.f_sam = evalFcn(objFcn,smart{ii}.x_sam,'objective');
    end
end

function [x_sel,f_sel]=selection(x,f)
    trnc_factor = 0.8;
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

function [parentclusters,smartclusters] = clustering(x_sel,f_sel,constFcn)
        nsel = size(x_sel,1);
        
        flag = false; % Found best cluster
        
        for jj=1:nsel
            [idx,mu] = kmeans(x_sel,jj);
            nc = nnz(evalFcn(constFcn,mu,'constraint')<=0); % Number of good clusters
            
            if nc==jj % If all clusters are good
                flag = true;
                break;
            end
        end        
        if ~flag
            error('mage mishe?')
        end
        parentclusters = cell(1,nc);
        % Saving parentclusters
        for ii=1:nc
            parentclusters{ii}.x = x_sel(idx==ii,:);
            parentclusters{ii}.f = f_sel(idx==ii);
            parentclusters{ii}.c = mu(ii,:);
        end
        
        % Until here, everything went good.
        % Now, check for good outliers
        
        lambda = 0.1; % Percentage of top good points to check in each cluster
        psi = 1; % Outlier distance for check (1sigma, 2sigma, 3sigma)
        smartclusters = cell(1,nc);
        counter = 0;
        for ii=1:nc
            curx = x_sel(idx==ii,:);
            curf = f_sel(idx==ii,:);
            curmu = mu(ii,:);
            cursigma = std(curx);
            curn = size(curx,1);

            % Sort this cluster
            [curf,curind] = sort(curf);
            curx = curx(curind,:);
            
            % Extract 10 of good points to check
            nbest = ceil(lambda*curn); % if it is 0 it will be disabled, if it is 1, all points will be checked
            bestx = curx(1:nbest,:);
            bestf = curf(1:nbest,:);

            outliersindx = (sqrt(sum((bestx-repmat(curmu,nbest,1)).^2,2)) > sqrt(sum((psi.*cursigma).^2,2)));
            if isempty(find(outliersindx,1))
                % No outliers, then no smartclusters
            else            
                counter = counter+1;
                % Create a new smart cluster
                cursmart = struct;
                cursmart.x = bestx(outliersindx,:);
                cursmart.f = bestf(outliersindx);
                cursmart.parent = curmu;
                
                smartclusters{counter} = cursmart;
            end
        end
        
        % Trash the empty at the tail
        smartclusters = smartclusters(1:counter);

        %%%%%%%%%%% Now combine two clusters
end

function [parentclusters,num_of_ref] = refinement(parentclusters,x_sel,f_sel)
    num_of_ref = 0;
    for ii=1:numel(parentclusters)
        if size(parentclusters{ii}.x,1)==1    
            num_of_ref = num_of_ref+1;

            % First, lets remove the current point from the search space
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
                % error('Something is wrong here!')
                % This should never happen!! The cluster has one point and it should be the centroid
                % If it is happenning, it maybe due to kmeans++ performance
                % If happened, no need to remove from search space
                x_rem = x_sel; % previous command line also works anyways
                f_rem = f_sel; % previous command line also works anyways
            end

            d = x_rem - repmat(parentclusters{ii}.x,size(x_rem,1),1);
            [~,q2] = min(sum(d.^2,2));
            thepointtoaddtothiscluster = x_rem(q2,:);
            thepointtoaddtothisclusterf = f_rem(q2);

            parentclusters{ii}.x = [parentclusters{ii}.x; thepointtoaddtothiscluster];
            parentclusters{ii}.f = [parentclusters{ii}.f; thepointtoaddtothisclusterf];
            
            %% it does NOT update center point
            
        end
    end
end

function [x_new,f_new,x_best,f_best] = replacement(parentclusters,smartclusters,x,f)
    % Best elitism    
    npop = numel(f);

        % Gathering
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
    
    % Just for best elitism (be careful here)
    x_best = x_new(1,:);
    f_best = f_new(1);

end

