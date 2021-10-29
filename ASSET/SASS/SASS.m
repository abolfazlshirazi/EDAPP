function [x,fx,gx,hx, convx,best_sol,best_of,best_conv, archive,hist_pos,memory_size, archive_c,archive_rp,current_eval] = ...
    SASS( x, fx,gx, hx, convx,best_sol,best_of,best_conv, archive,hist_pos,memory_size, archive_c,archive_rp, xmin, xmax,  n,...
    PopSize,  current_eval, I_fno, gg, const_num, Eg, Eh, Par)

ui=zeros(PopSize,n);

%% calc CR and F
mem_rand_index = ceil(memory_size * rand(PopSize, 1));
mu_c = archive_c(mem_rand_index);
mu_rp = archive_rp(mem_rand_index);
%% ========================= generate rp ==================================
rp = (mu_rp + 0.1*sqrt(pi)*(asin(-rand(1,PopSize))+asin(rand(1,PopSize))))';
rp(mu_rp == -1) = 0;
rp = min(rp, 1);
rp = max(rp, 0);

%% ========================= generate c ===================================
    c = mu_c + 0.1 * tan(pi * (rand(1,PopSize) - 0.5));
    pos = find(c <= 0);
    while ~ isempty(pos)
        c(pos) = mu_c(pos) + 0.1 * tan(pi * (rand(1,length(pos)) - 0.5));
        pos = find(c <= 0);
    end
    c = min(c, 1)';
%% ======================== generate new x =================================
popAll = [x;archive.pop]; %% set archive
r0 = 1 : PopSize;
%% generate random integer numbers
[r1, r2] = gnR1R2(PopSize, size(popAll, 1), r0);

%% calculation of pbest solution
pNP = max(round(0.1 * PopSize), 2); %% choose at least two best solutions
randindex = ceil(rand(1, PopSize) .* pNP); %% select from [1, 2, 3, ..., pNP]
randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
phix = x(randindex, :); %% randomly choose one of the top 10% solutions

%% calculation binary diagonal matrix
mask = rand(PopSize, n) < rp(:, ones(1, n)); % mask is used to indicate which elements of ui comes from the parent
rows = (1 : PopSize)'; cols = floor(rand(PopSize, 1) * n)+1; % choose one position where the element of ui doesn't come from the parent
jrand = sub2ind([PopSize n], rows, cols); mask(jrand) = true;
%% spherical search
zi = phix - x + x(r1, :) - popAll(r2, :);
%% Binary Orthogonal Matrix
A  = RandOrthMat(n,0);
for i = 1:PopSize
    B = diag(mask(i, :));
    ui(i, :) = x(i, :) + c(i).*zi(i, :)*A*B*A';
end
ui = han_boun(ui, xmax, xmin, x, PopSize,1);
%% evaluate
[fx_new,gx_new,hx_new] = cec20_func(ui,I_fno);fx_new = fx_new';
convx_new = (sum([sum(gx_new.*(gx_new>0),1); sum(abs(hx_new).*(abs(hx_new)>1e-4),1)],1)./const_num);
%% update FITNESS EVALUATIONS
current_eval =current_eval+PopSize;
%% Repair Solution
problem.constr_fun_name = @cec20_func;
problem.gn              = Par.gn;
problem.hn              = Par.hn;        
problem.I_fno           = I_fno;
problem.xmin            = Par.xmin;
problem.xmax            = Par.xmax;
problem.n               = n;
if (problem.hn(I_fno) > 0 && Par.n < 20) || problem.hn(I_fno) == 0
    problem.maxiter        = 4;
else
    problem.maxiter        = 3000;
end
for i = 1:PopSize
    if rand < 0.2 && convx_new(i) ~= 0 && mod(gg,n) == 0 && convx_new(i) < 1e2  
          [new_mutant,fes]     = QNrepair(problem, ui(i,:)', gx_new(:,i), hx_new(:,i));
          new_mutant          = han_boun(new_mutant', xmax, xmin, new_mutant',1,3);
          current_eval        = current_eval+fes;
          [fval, gv, hv]      = cec20_func(new_mutant,I_fno);
          fx_new(i)           = fval;                                                               
          convx_new(i)        = (sum([sum(gv.*(gv>0),1); sum(abs(hv).*(abs(hv)>1e-4),1)])./const_num);
          gx_new(:,i)         = gv;
          hx_new(:,i)         = hv;
          ui(i,:)             = new_mutant;
    end       
end

%% calc. imprv. for RP and C
diff = abs(fx - fx_new);
I =(fx_new < fx);
goodRP = rp(I == 1);
goodC = c(I == 1);

%% ========================= update archive ===============================
archive = updateArchive(archive, x(I == 1, :));
%% =================== update memory c and rp =============================
num_success_params = numel(goodRP);
if num_success_params > 0
    weightsSS = diff(I == 1)./ sum(diff(I == 1));
    %% for updating the memory of C
    archive_c(hist_pos) = (weightsSS * (goodC .^ 2))./ (weightsSS * goodC);
    
    %% for updating the memory of RP
    if max(goodRP) == 0 || archive_rp(hist_pos)  == -1
        archive_rp(hist_pos)  = -1;
    else
        archive_rp(hist_pos) = (weightsSS * (goodRP .^ 2)) / (weightsSS * goodRP);
    end
    
    hist_pos= hist_pos+1;
    if hist_pos > memory_size;  hist_pos = 1; end
end
%% ==================== update x and fitx =================================
I = comp_sol(fx_new,fx,gx_new,gx,hx_new,hx,Eg,Eh);
x(I == 1, :) = ui(I == 1, :);
fx( I == 1) = fx_new(I == 1);
gx(:, I == 1) = gx_new(:, I == 1);
hx(:, I == 1) = hx_new(:, I == 1);
convx(I == 1) = convx_new(I == 1);
%% sort new x, fitness
ind = eps_sort(fx,convx,0);
x  = x(ind,:);
fx = fx(ind);
gx = gx(:, ind);
hx = hx(:, ind);
convx = convx(ind);
%% update best solution found so far   
        if (convx(ind(1))==0 && best_conv==0 && fx(ind(1)) <= best_of) ||...
           (convx(ind(1)) == best_conv && fx(ind(1)) <= best_of) || convx(ind(1)) < best_conv
            best_sol   = x(ind(1), :); 				
            best_of =  fx(ind(1));
            best_conv = convx(ind(1));
        end

%
end
