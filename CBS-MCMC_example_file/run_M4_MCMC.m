clear all; close all
%% Settings of constrained based parameter search
<<<<<<< HEAD
NMC = 5e5; %5e5; % number of parameters in initial sampling
N_new = 1000;%2e3; %newly generated samples
N_add = 5e3;%5e3; %added samples if no samples fulfill the rules
N_tar = 10e4;%10e4; %finally samples following target distribution
=======
NMC = 5e5; % number of parameters in initial sampling
N_new = 1000; %newly generated samples
N_add = 5e3; %added samples if no samples fulfill the rules
N_tar = 1e5; %finally samples following target distribution
>>>>>>> 601373ae4ef77a4f6a20f9018b7ad66e4045ff9b
n_chain = 40; % number of Markov chains
n_min = n_chain; % minimum number of parameter sets that fulfill the rules for generating new proposals
N_tot = N_tar*5; % total number of function evaluations
% f=0; %factor to define how many less process constrains need to be accepted
% n_min=1; % minimum of accepted parameters that need to be accepted

%% cluster settings
cluster=1; % specifiy if cluster is used (1) or not (0)

if cluster
    % create a local cluster object
    pc = parcluster('local');
    
    % get the number of dedicated cores from environment
    num_workers = str2num(getenv('SLURM_NPROCS'));
    
    % explicitly set the JobStorageLocation to the tmp directory that is unique to each cluster job (and is on local, fast scratch)
    parpool_tmpdir = [getenv('TMP'),'/.matlab/local_cluster_jobs/slurm_jobID_',getenv('SLURM_JOB_ID')];
    mkdir(parpool_tmpdir);
    
    pc.JobStorageLocation = parpool_tmpdir;
    % start the parallel pool
    parpool(pc,num_workers,'AttachedFiles',{'run_M4_n.m','check_process_constrains_chem.m',...
        'M4.m','AT_init.m','BoundaryHandling.m','check_parameter_constrains_chem.m',...
        'LHS.m','stopevent.m'})
    
elseif ~cluster
    if isempty(gcp('nocreate'))% check if parallel pool is already runnung
        % if not start with two workers
        pc=parcluster('local');
        parpool(pc,'AttachedFiles',{'run_M4_n.m','check_process_constrains_chem.m',...
            'M4.m','AT_init.m','BoundaryHandling.m','check_parameter_constrains_chem.m'...
            'LHS.m','stopevent.m'})
    else
        % if yes attach and update files
        addAttachedFiles(gcp,{'run_M4_n.m','check_process_constrains_chem.m',...
            'M4.m','AT_init.m','BoundaryHandling.m','check_parameter_constrains_chem.m',...
            'LHS.m','stopevent.m'});
        updateAttachedFiles(gcp);
    end
end

%% initial parameter sampling
% read parameter bounds
pb=csvread('par_bounds_M4.csv');

% transform parameter space such that all parameters are>=1
p_trans=abs(pb(:,1))+1;
parmin =	(pb(:,1)+p_trans)';
parmax =	(pb(:,2)+p_trans)';

ParRange.minn = parmin; ParRange.maxn = parmax;

%% implement a constraint-based search algorithm %%
C=7; % number of rules to be fullfilled
t_points=401; % number of timepoints rejected from model function
n_mo=15; % number of model outputs;
n_par=length(parmin); % number of parameters

% A constraint-based search algorithm for parameter identification of environmental models (2014)
% Hydrology and Earth System Sciences 18:4861?4870 http://dx.doi.org/10.5194/hess-18-4861-2014
% Gharari S, Shafiei M, Hrachowitz M, Kumar R, Fenicia F, Gupta HV, Savenije HHG

while 1
    % sample parameters using latin hypercube sampling
    x = LHS(parmin,parmax,NMC);
    % step 1 : check parameter constrains
    % select parameters:
    % growth rates higher than dead rates, and dead rate of active higher than
    % dormant bacteria (line 1-3);
    sel=x(:,2)>=x(:,9) & x(:,9)>=x(:,10)... %bacteria A
        & x(:,20)>=x(:,23) & x(:,23)>=x(:,24)...%bacteria B
        & x(:,42)>=x(:,44) & x(:,44)>=x(:,45)...%bacteria D
        & x(:,59)>max(x(:,[58,60:64]),[],2) & x(:,58)> max(x(:,[60,62,64]),[],2); %order of sorption coefficients: HY highest; AT>DEA,NE,CY
    
    % define selected parameter sets that fullfill parameter constrains
    x_sel=x(sel,:);
    
    % Use function from DREAM for boundary handling and reflect parameter
    % values outside of ranges
    [x_sel] = BoundaryHandling(x_sel,ParRange,'Reflect');
    
    if ~isempty(x_sel)
        % estimate number of paramter sets that fullfill the parameter constrains
        NMC_sel=size(x_sel,1);
        idx=1:NMC_sel;
        c_sim = NMC_sel;
        
        % initialize variable that stores the number of accepted constrains
        n_sel = zeros(NMC_sel,1);
        n_idx_sel = zeros(NMC_sel,C);
        fx_sel = zeros(NMC_sel,t_points,n_mo);
        
        % retransform parameters
        x_run=x_sel-p_trans';
        
        parfor ii=1:NMC_sel
            [n_sel(ii),n_idx_sel(ii,:),fx_sel(ii,:,:)]=run_M4_n(x_run(ii,:)); % run models and estimated how many process rules were fullfilled
        end
        
        if length(n_sel(n_sel>=2))>=n_min % make sure we get n_min parameter set that fulfill at least two rules at the beginning with Monte Carlo simulation
            break;
        end
    end
    
    NMC=NMC+NMC;
end

% initialize variables that store initial sets, factor 2 to ensure that
% they are large enough
%x_fin=zeros(2*N_fin, n_par); % initialize final parameter set
%n_fin=zeros(2*N_fin, C); % initilize final rules accepted
%fx_fin=zeros(2*N_fin,t_points,n_mo);% initialize model outputs
% add current time to generate unique output files
f_t=datestr(datetime('now'),'yymmdd_HHMMSS');
f_out=['M4_out_',f_t,'.h5'];

% create data sets and hdf5 file to store output
% h5create(f_out,'/x_fin',[Inf, n_par],'ChunkSize',[1 n_par]);
% h5create(f_out,'/n_fin',[Inf, C],'ChunkSize',[1 C]);
% h5create(f_out,'/fx_fin',[Inf, t_points, n_mo],'ChunkSize',[1 t_points n_mo]);
% br=0; % set break boolean to decide if simulations are stopped
N_MCMC=N_new; % initialize N_MCMC
[idx_seed,x_seed,sumd,D,ind_seed] = kmedoids(x_sel,n_chain);
sigma_p = std(x_sel,0,1);

for i_C=2:C % iterate number of accepted rules
    
    suc = 0; % set boolean indicating success of current iteration
    n_acc = 0; % number of accepted samples
    n_new_run = 0; % number of model runs for the new proposed samples in each iteration
    
    while ~suc && c_sim<=N_tot
        n_add_run = 0; % number of model runs for addtionally added samples
        % estimate indices
        sel_1 = idx(n_sel>=i_C); % parameter sets that fullfill i_C or more rules
        n_unique = size(unique(x_sel(sel_1,:),'row'),1); % determine size of unique parameter sets
        
        if n_unique>=n_min % at least n_min parameter set fullfull i_C or more rules
            
            suc=1;
            
            num_sel = length(sel_1); % number of parameter sets that currently fullfill i_C or more rules, number of markov chains
            
            x_1=x_sel(sel_1,:); % parameter sets that fullfil C or more rules
            n_1 = n_sel(sel_1);
            n_idx_1 = n_idx_sel(sel_1,:);
            fx_1 = fx_sel(sel_1,:,:);
            
            n_acc = n_acc + num_sel; % currently accepted parameters
            
            N_MCMC = N_new;
            
            j_r=1;% set jump rate

            while c_sim<=N_tot
                [idx,x_seed,sumd,D,ind_seed] = kmedoids(unique(x_1,'rows'),n_chain);
                n_seed = n_1(ind_seed);
                n_idx_seed = n_idx_1(ind_seed,:);
                fx_seed = fx_1(ind_seed,:,:);
                
                % number of samples simulated from each Markov chain
                n_sample = ceil(N_MCMC/n_chain);
                
                N_MCMC = n_sample*n_chain;
                n_new_run = N_MCMC;
                
                % initialize sample set and outputs
                x_mcmc = zeros(n_sample,n_par,n_chain);
                fx_mcmc = zeros(n_sample,t_points,n_mo,n_chain);
                n_mcmc = zeros(n_sample,n_chain);
                n_idx_mcmc = zeros(n_sample,C,n_chain);
                
                sigma_p = j_r*std(unique(x_1,'rows'),0,1);
                
                % initialize length of markov chain
                len_each_chain = 1;
                while len_each_chain <= n_sample
                    % generate candidates
                    x_can_1 = x_seed;
                    for i_chain = 1:n_chain
                        % candidate sample
                        while 1
                            x_can = zeros(1,n_par);
                            for i_par = 1:n_par
                                % propose an candidate
                                
                                if i_C == C
                                    %x_can(i_par) = x_seed(i_chain,i_par) + 0.3*sigma_p(i_par)*randn;
                                    x_can(i_par) = x_seed(i_chain,i_par) + sigma_p(i_par)*randn;
                                else
                                    %x_can(i_par) = x_seed(i_chain,i_par) + 0.5*sigma_p(i_par)*randn;
                                    x_can(i_par) = x_seed(i_chain,i_par) + sigma_p(i_par)*randn;
                                end
                                % Use function from DREAM for boundary handling and reflect parameter
                                % values outside of ranges
                            end
                            [x_can] = BoundaryHandling(x_can,ParRange,'Reflect');
                            if check_parameter_constrains_chem(x_can)
                                x_can_1(i_chain,:)= x_can;
                                break; % make sure the new sample fulfill the parameter constraints
                            end
                        end
                    end
                    
                    parfor i_chain = 1:n_chain
                        x_mcmc(len_each_chain,:,i_chain) = x_can_1(i_chain,:);
                        
                        % retransform parameters
                        x_mcmc_run = x_mcmc(len_each_chain,:,i_chain)-p_trans';
                        % run the model for the proposed new sample
                        [n_mcmc(len_each_chain,i_chain),n_idx_mcmc(len_each_chain,:,i_chain),fx_mcmc(len_each_chain,:,:,i_chain)] = run_M4_n(x_mcmc_run);
                        if n_mcmc(len_each_chain,i_chain) >= i_C % accept, update the seeds
                            x_seed(i_chain,:) = x_mcmc(len_each_chain,:,i_chain);
                            n_seed(i_chain) = n_mcmc(len_each_chain,i_chain);
                            n_idx_seed(i_chain,:) = n_idx_mcmc(len_each_chain,:,i_chain);
                            fx_seed(i_chain,:,:) = fx_mcmc(len_each_chain,:,:,i_chain);
                            
                            % len_each_chain(i_chain) = len_each_chain(i_chain) + 1;
                            
                        else % reject, keep the current state
                            x_mcmc(len_each_chain,:,i_chain) = x_seed(i_chain,:);
                            n_mcmc(len_each_chain,i_chain) = n_seed(i_chain);
                            n_idx_mcmc(len_each_chain,:,i_chain) = n_idx_seed(i_chain,:);
                            fx_mcmc(len_each_chain,:,:,i_chain) = fx_seed(i_chain,:,:);
                            
                            % n_new_run = n_new_run + 1; % record the model runs
                        end
                        
                        
                    end
                    len_each_chain = len_each_chain + 1;
                end
                
                % reshape these sample sets
                x_sel = zeros(N_MCMC+num_sel,n_par);
                fx_sel = zeros(N_MCMC+num_sel,t_points,n_mo);
                n_sel = zeros(N_MCMC+num_sel,1);
                n_idx_sel = zeros(N_MCMC+num_sel,C);
                for i = 1:n_chain
                    x_sel((i-1)*(n_sample)+1:i*n_sample,:) = x_mcmc(:,:,i);
                    fx_sel((i-1)*(n_sample)+1:i*n_sample,:,:) = fx_mcmc(:,:,:,i);
                    n_sel((i-1)*(n_sample)+1:i*n_sample) = n_mcmc(:,i);
                    n_idx_sel((i-1)*(n_sample)+1:i*n_sample,:) = n_idx_mcmc(:,:,i);
                end
                x_sel(N_MCMC+1:end,:) = x_1;
                fx_sel(N_MCMC+1:end,:,:) = fx_1;
                n_sel(N_MCMC+1:end) = n_1;
                n_idx_sel(N_MCMC+1:end,:) = n_idx_1;
                
                
                c_sim = c_sim + n_new_run;
                
                if i_C < C
                    % estimate number of paramter sets that fullfill the parameter constrains
                    NMC_sel=size(x_sel,1);
                    % update idx vector
                    idx=1:NMC_sel;
                end
                
                a_r=size(unique(x_sel,'row'),1)/N_MCMC; % calculate proportion of uniquely accepted parameters
                
                prog=sprintf('i_C %7.0f | a_r %3.0f | N_MCMC %7.0f | n_new_run %7.0f | n_add_run %7.0f | c_sim %7.0f',[i_C, a_r*100,N_MCMC,n_new_run,n_add_run,c_sim]);
                disp(prog)

                if i_C<C
                    break;
                else % adapt jump rate to optimize acceptance rate with short MCMCs
                    if N_MCMC == N_tar
                        break;
                    elseif a_r <0.2
                        j_r=j_r/1.3;
                    elseif a_r > 0.5
                        j_r=j_r*1.3;
                    else
                        N_MCMC = N_tar;
                    end
                end
            end % tuning of acceptance rate in last MCMC
        else % less then n_min parameter set fulfills i_C or more rules, try to get more parameter sets from the previous iteration with i_C-1
            
            % number of samples simulated from each Markov chain
            n_sample = ceil((N_MCMC+N_add)/n_chain);
            
            N_MCMC = n_sample*n_chain;
            n_add_run = N_MCMC;
            
            % initialize sample set and outputs
            x_mcmc = zeros(n_sample,n_par,n_chain);
            fx_mcmc = zeros(n_sample,t_points,n_mo,n_chain);
            n_mcmc = zeros(n_sample,n_chain);
            n_idx_mcmc = zeros(n_sample,C,n_chain);
            
            % initialize length of markov chain
            len_each_chain = 1;
            while len_each_chain <= n_sample
                % generate candidates
                x_can_1 = x_seed;
                for i_chain = 1:n_chain
                    % candidate sample
                    while 1
                        x_can = zeros(1,n_par);
                        for i_par = 1:n_par
                            % propose an candidate
                            x_can(i_par) = x_seed(i_chain,i_par) + sigma_p(i_par)*randn;
                            
                        end
                        % Use function from DREAM for boundary handling and reflect parameter
                        % values outside of ranges
                        [x_can] = BoundaryHandling(x_can,ParRange,'Reflect');
                        if check_parameter_constrains_chem(x_can)
                            x_can_1(i_chain,:)= x_can;
                            break; % make sure the new sample fulfill the parameter constraints
                        end
                    end
                end
                
                parfor i_chain = 1:n_chain
                    
                    x_mcmc(len_each_chain,:,i_chain) = x_can_1(i_chain,:);
                    
                    % retransform parameters
                    x_mcmc_run = x_mcmc(len_each_chain,:,i_chain)-p_trans';
                    % run the model for the proposed new sample
                    [n_mcmc(len_each_chain,i_chain),n_idx_mcmc(len_each_chain,:,i_chain),fx_mcmc(len_each_chain,:,:,i_chain)] = run_M4_n(x_mcmc_run);
                    if n_mcmc(len_each_chain,i_chain) >= i_C-1 % accept, update the seeds
                        x_seed(i_chain,:) = x_mcmc(len_each_chain,:,i_chain);
                        n_seed(i_chain) = n_mcmc(len_each_chain,i_chain);
                        n_idx_seed(i_chain,:) = n_idx_mcmc(len_each_chain,:,i_chain);
                        fx_seed(i_chain,:,:) = fx_mcmc(len_each_chain,:,:,i_chain);
                        
                        % len_each_chain(i_chain) = len_each_chain(i_chain) + 1;
                        
                    else % reject, keep the current state
                        x_mcmc(len_each_chain,:,i_chain) = x_seed(i_chain,:);
                        n_mcmc(len_each_chain,i_chain) = n_seed(i_chain);
                        n_idx_mcmc(len_each_chain,:,i_chain) = n_idx_seed(i_chain,:);
                        fx_mcmc(len_each_chain,:,:,i_chain) = fx_seed(i_chain,:,:);
                        
                        % n_add_run = n_add_run + 1; % record the model runs
                    end
                    
                    
                end
                len_each_chain = len_each_chain + 1;
            end
            
            % reshape these sample sets
            x_sel = zeros(N_MCMC,n_par);
            fx_sel = zeros(N_MCMC,t_points,n_mo);
            n_sel = zeros(N_MCMC,1);
            n_idx_sel = zeros(N_MCMC,C);
            for i = 1:n_chain
                x_sel((i-1)*n_sample+1:i*n_sample,:) = x_mcmc(:,:,i);
                fx_sel((i-1)*n_sample+1:i*n_sample,:,:) = fx_mcmc(:,:,:,i);
                n_sel((i-1)*n_sample+1:i*n_sample) = n_mcmc(:,i);
                n_idx_sel((i-1)*n_sample+1:i*n_sample,:) = n_idx_mcmc(:,:,i);
            end
            
            c_sim = c_sim + n_add_run;
            
            % estimate number of paramter sets that fullfill the parameter constrains
            NMC_sel=size(x_sel,1);
            % update idx vector
            idx=1:NMC_sel;
            
            a_r=size(unique(x_sel,'row'),1)/N_MCMC; % calculate proportion of uniquely accepted parameters
            
            prog=sprintf('i_C %7.0f | a_r %3.0f | N_MCMC %7.0f | n_new_run %7.0f | n_add_run %7.0f | c_sim %7.0f',[i_C, a_r*100,N_MCMC,n_new_run,n_add_run,c_sim]);
            disp(prog)
            
        end
        
    end
end

h5create(f_out,'/x_fin',size(x_sel));
h5create(f_out,'/fx_fin',size(fx_sel),'ChunkSize',[1 t_points n_mo]);
h5create(f_out,'/n_fin',size(n_sel));

h5write(f_out, '/x_fin', x_sel-p_trans');
h5write(f_out, '/fx_fin', fx_sel);
h5write(f_out,'/n_fin',n_sel);