% %% CODE
% % Initial input data
% % Add paths
% % Input "ctrl" variable
% % Read input files
% % Compute idx for the different alternatives of time
% % Sampling combinations (or designs)
% % Apply PreDIA
% % Compute BME diagonal
% % Save results
% % Plot results
% 
% %% Initial input data
clc
clear
format compact
% multiplied by 4, because new outputs are every 0.25 days
t1 = 25*4 + 1;                 % Option 1: one sample every day 
t2 = 50*4 + 1;                 % Option 2: one sample every two-days
t4 = 100*4 + 1;                % Option 4: one sample every four-days

n_obs = 6%6                    % 6 types of observations =[Atrazine, CO2, Cyanuric, DEA, DIA, Hydroxyatrazine]
n_mod = 6%6                    % number of models
n_mc  = 100%50000                % NMC realizations for each model

meas_err = [0.10 0.10 0.10 0.10 0.10 0.10]';    % relative error = 10% for each observation type

t_max = 3                      % number of times series considered for the designs
flag_t = 1;                    % if =0 n_time=1,     (x,x,     1,x)
                               % if =1 n_time=n_time (x,x,n_time,x), times are separated
n_print = 5;                   % save results every n_print designs                     

% 
% %% ADD PATHS
if isunix
    main_path = pwd;
%    addpath([main_path ':' main_path '/aaux']); %ag addpath([main_path ':' main_path '/aux']);
%     addpath([main_path ':' main_path '/../outputs']);
elseif ispc
%    addpath ('aaux\') %ag addpath ('aux\')
%     addpath ('..\outputs\')
end


 

%% Leonardo READ INPUT DATA
% Needs to be modified for new input from Luciana M1_1, M1_2, M1_3, etc.
% if in windows
% cd  'D:\PROJECTS\Hohenheim_Atrazine_ODE_Luciana\OUTPUTS_Luciana\2019-03-06_Outputs_Atrazine_final'
if isunix
%     cd '/nfs/home_simtech/gonzalez/001_BME/Outputs/1_M1-M6'
    cd 'C:/Doctorado/Programa/Nuevo_PECCAD/NEW ATRAZINE/ENERGY_DISTANCE/Outputs/1_M1M6'
elseif ispc
    cd  'C:/Doctorado/Programa/Nuevo_PECCAD/NEW ATRAZINE/ENERGY_DISTANCE/Outputs/1_M1M6'
end


tic
list = dir('*.mat');
for n = 1:size(list)
    load(list(n).name);
end
toc
tic
gg = list;
gg = who('M*');  % need to do this
% from M1* to M6 ag(:,i)=who('Mi*')
ag(:,1) = who('M1*');
ag(:,2) = who('M2*');
ag(:,3) = who('M3*');
ag(:,4) = who('M4*');
ag(:,5) = who('M5*');
ag(:,6) = who('M6*');
toc

%% Compute idx for the different alternatives of time
% daily, every two days, every three days, every four days
time = 1:1:t1-1; 
idx1 = rem(time,1*4)== 0; s = t4-t1; idx = zeros(1,s); idx1 = [idx1 idx]; idx_t(1,:) = logical(idx1);
time = 1:1:t2-1;
idx2 = rem(time,2*4)== 0; s = t4-t2; idx = zeros(1,s); idx2 = [idx2 idx]; idx_t(2,:) = logical(idx2);
time = 1:1:t4-1;
idx_t(3,:) = rem(time,4*4)== 0; %s = t4-t2; idx = zeros(1,s); idx2 = [idx2 idx];

%% SAMPLING COMBINATIONS (OR DESIGNS)
sum_d =0;
for k = 1:n_obs
    k;
    temp=nchoosek(1:n_obs,k);
    size_d(k)=size(temp,1);
    sum_d = sum_d + size(temp,1);
    namef = ['design_' num2str(k)];
    save(namef,'temp')
end
n_d = sum_d*t_max

% Write down all possible designs
n_d_test = 0;
model_design = zeros(n_d,n_obs+1);
for t=1:t_max
    for k = 1:n_obs
        combinations = nchoosek(1:n_obs,k);
        for l=1:size(combinations,1)
            n_d_test = n_d_test+1;
            for m=1:size(combinations,2)
                model_design(n_d_test,combinations(l,m)) = 1;
                model_design(n_d_test,n_obs+1) = t;
            end
        end
    end
end
n_obs_d = sum(model_design(:,1:n_obs),2); % will tell how many types of observations are "measured" in "d"
%% Energy distances

%%% change this pathway
if isunix
%     cd '/nfs/home_simtech/gonzalez/001_BME/1_M1_M6'
    cd 'C:/Doctorado/Programa/Nuevo_PECCAD/NEW ATRAZINE/ENERGY_DISTANCE'    
    elseif ispc
    cd  'C:/Doctorado/Programa/Nuevo_PECCAD/NEW ATRAZINE/ENERGY_DISTANCE'
end

% Initialize matrices
n_d = n_d_test;
n_time = sum(idx1);                     % # observations along time
if flag_t ==1   
    E_dist = nan(n_mod,n_mod,n_time,n_d);   % =1 times are separated
else
    E_dist = nan(n_mod,n_mod,n_d);          % =0 there is only one time
end
E_dist2 = nan(n_mod,n_mod,n_d);
Max_Ed = nan(1,n_d);
Max_Ed_mod = nan(n_mod,n_d);

tic
for d = 1:7%n_d %parfor
    d
    ttime = cputime;

    % Define intervals to read for each design
    interval = '(idx_t(model_design(d,n_obs+1),:),1:n_mc)';
    
    % Load data for input and hypothetical observations
    n_comb = 0;
    n_comb = n_obs_d(d);    % will tell how many types of observations are "measured" in "d"

    if flag_t == 1
       comb_input = nan(n_mc,n_comb,n_time,n_mod);      % =1 times are separated
    else
       comb_input = nan(n_mc,n_comb*n_time,1,n_mod);    % =0 there is only one time 
    end
    
    for m=1:n_mod
        nn = 0; n1 = 0; n2 = 0;
        for k=1:n_obs
            if (model_design(d,k)==1)
                  nn = nn + 1;           % counter for # of observations type
                  n1 = n_time*(nn-1)+1;  % ini
                  n2 = n_time*nn;        % end
                  clear temp
                  temp = eval(strcat(ag{k,m},interval));
                  temp = temp';
                  if flag_t == 1          % each time step is separated
                      for t=1:n_time
                          comb_input(:,nn,t,m) = temp(:,t);
                      end
                  else                    % there is only one time 
                      comb_input(:,n1:n2,1,m) = temp;
                  end
                  % Normalize by measurement error
                  comb_input(:,nn,:,m) = comb_input(:,nn,:,m)./meas_err(nn);
            end
        end
    end
 
    %%% CALL Energy Distance HERE

    if flag_t == 1  % =1 each time step is separated
        E_dist(:,:,:,d)= EnergyDistance(comb_input,n_mc,n_comb,n_time,n_mod,flag_t);
    else            % =0  only one time step
        E_dist(:,:,d)= EnergyDistance(comb_input,n_mc,n_comb,n_time,n_mod,flag_t);
    end

  
    if flag_t ==1   
        E_dist2(:,:,d) = sum(E_dist(:,:,:,d),3)/n_time;    % average for all times
    else
        E_dist2(:,:,d) = E_dist(:,:,d);
    end
    
    % Look for maximums
    Max_Ed(d) = sum(sum(E_dist2(:,:,d)));
    Max_Ed(d);
    temp2 = tril(E_dist2(:,:,d).',-1) + triu(E_dist2(:,:,d));   %mirror matrix to have all components of the matrix
    Max_Ed_mod(:,d) =sum(temp2,2);

    time_d(d,1) = cputime - ttime;
    % Generate temporary results (in case of Re-start)
    if rem(d,n_print)==0 || d==1 || d==n_d;
        d
        save Results_time.mat time_d
        save Results_Max_Ed.mat Max_Ed
        save Results_Max_Ed_models.mat Max_Ed_mod
    end
    
    clear temp temp2
end              
toc 
% Look for global maximums
[Max_Ed_d,d1] = max(Max_Ed);
[Max_Ed_mod_d,d2] = max(Max_Ed_mod');
Best_d = [d1 Max_Ed_d];
Best_d_model = [d2' Max_Ed_mod_d'];

save Best_d.mat Best_d
save Best_d_model.mat Best_d_model

