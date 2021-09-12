function [n,n_idx,fx] = run_M4_n(p,~)
% p parameter vector
% output: number of fullfilled rules
%% Define fixed parameters %%

q(1)        = 216;     % Molar Weight of AT (mg AT mol^-1)
q(2)        = 197;     % Molar Weight of HY (mg HY mol^-1)
q(3)        = 156;     % Molar Weight of NE (mg NE mol^-1)
q(4)        = 170;     % Molar Weight of NI (mg NI mol^-1)
q(5)        = 188;     % Molar Weight of DEA (mg DEA mol^-1)
q(6)        = 174;     % Molar Weight of DIA (mg DIA mol^-1)
q(7)        = 129;     % Molar Weight of CYA (mg CYA mol^-1)
q(8)        = 12;      % Molar Weight of Carbon (mg C mol^-1)
q(9)        = 8;       % Number of C atoms in AT
q(10)       = 8;       % Number of C atoms in HY
q(11)       = 5;       % Number of C atoms in NE
q(12)       = 6;       % Number of C atoms in NI
q(13)       = 6;       % Number of C atoms in DEA
q(14)       = 5;       % Number of C atoms in DIA
q(15)       = 3;       % Number of C atoms in CYA
q(16)       = 0.1;     % Switch function parameter
q(17)       = 0.35;    % Average volumetric soil water content (cm^3 cm^-3)
q(18)       = 1.1;     % Bulk density of soils (g cm^-3)
q(19)       = 3/8;     % Stoichiometric formation of CYA from AT
q(20)       = 3/8;     % Stoichiometric formation of CYA from HY
q(21)       = 3/5;     % Stoichiometric formation of CYA from NE
q(22)       = 3/6;     % Stoichiometric formation of CYA from NI
q(23)       = 3/6;     % Stoichiometric formation of CYA from DEA
q(24)       = 3/5;     % Stoichiometric formation of CYA from DIA
q(25)       = 5/8;     % Stoichiometric formation of NE from HY
q(26)       = 6/8;     % Stoichiometric formation of NI from HY

%% set ODE solver options %%

abstol = 1e-9;
reltol = 1e-7;
o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent,'NonNegative',1:14);
t = linspace(0,100,401); % simulation period (d)

% INITIAL CONDITION

K_F     = p(58);
n_F     = p(65);
C_T     = 1e-01;   % mg AT g-1 soil
C_L0    = 0;       % mg AT cm3
th_V    = 0.35;    % Average volumetric soil water content (cm^3 cm^-3)
rho_B   = 1.1;     % Bulk density of soils (g cm^-3)
Bac     = 1e-03;   % Initial concentration of bacteria (mg C g^-1 soil)
N_C2    = 3;       % Number of bacteria involved

fun=@(C_L0) AT_init(C_L0,C_T,K_F,n_F,th_V,rho_B);
%  switch off solver progress information on solver progress of fsolve
options = optimoptions('fsolve','Display','none');
res = fsolve(fun,C_L0,options); % AT in solution

% INITIAL CONDITIONS

c(1)  = 0;           % Initial conc active Bacteria A (mg C g^-1 soil)
c(2)  = Bac/N_C2;    % Initial conc inactive Bacteria A (mg C g^-1 soil)
c(3)  = 0;           % Initial conc active Bacteria B (mg C g^-1 soil)
c(4)  = Bac/N_C2;    % Initial conc inactive Bacteria B (mg C g^-1 soil)
c(5)  = 0;           % Initial conc active Bacteria D (mg C g^-1 soil)
c(6)  = Bac/N_C2;    % Initial conc inactive Bacteria D (mg C g^-1 soil)
c(7)  = res;         % Initial Atrazine conc in Solution (mg AT cm^-3)
c(8)  = eps;         % Initial Hydroxyatrazine conc in Solution (mg HY cm^-3)
c(9)  = eps;         % Initial N-Ethylammelide conc in Solution (mg NE cm^-3)
c(10) = eps;         % Initial Deethylatrazine conc in Solution (mg DEA cm^-3)
c(11) = eps;         % Initial Deisopropylatrazine conc in Solution (mg DIA cm^-3)
c(12) = eps;         % Initial Cyanuric Acid conc in Solution (mg CYA cm^-3)
c(13) = 0;           % Initial CO2 conc (mg CO2 cm^-3)
c(14) = 0;           % DOC pool (mg C g^(-1)  soil)

ty=[]; cu=[];
try
    warning off
    tic
    [ty, cu] = ode45(@M4,t,c,o_opts,p',q); % ode15s
catch ME
    warning off
end

% set all outputs to zero and only fill them if model run was successful
fx=zeros(size(t,2),size(c,2)+1);
n=0; % total number of accepted rules
n_idx=zeros(7,1); % array indicating which rule was fullfilled

if (length(cu) == length(t) && isreal(cu))% check if model run was successful without imaginary results
    fx=[ty,cu];% store model output in fx
    [n,n_idx] = check_process_constrains_chem(ty,cu,C_T);    
end
