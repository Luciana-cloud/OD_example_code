function vf_ = M2(t,x,p,q)

%% STOCKS

Bac_B      = x(1); % Active Bacteria B [mg C g-1]
Bin_B      = x(2); % Inactive Bacteria B [mg C g-1]
Bac_C      = x(3); % Active Bacteria C [mg C g-1]
Bin_C      = x(4); % Inactive Bacteria C [mg C g-1]
Bac_D      = x(5); % Active Bacteria D [mg C g-1]
Bin_D      = x(6); % Inactive Bacteria D [mg C g-1]
CL_AT      = x(7); % ATrazine in solution [mg C cm-3]
CL_HY      = x(8); % ATrazine in solution [mg C cm-3]
CL_NE      = x(9); % ATrazine in solution [mg C cm-3]
CL_NI      = x(10); % ATrazine in solution [mg C cm-3]
CL_DEA     = x(11); % ATrazine in solution [mg C cm-3]
CL_DIA     = x(12); % ATrazine in solution [mg C cm-3]
CL_CA      = x(13); % ATrazine in solution [mg C cm-3]
% CO2        = x(14); % ATrazine in solution [mg C g-1]
% DOC        = x(15); % ATrazine in solution [mg C g-1]

%% PARAMETERS %%

%BACTERIA B

CT_AT_B   = 10^p(19); % Threshold AT conc. for growth Bacteria B (mg C cm^-3 water)
mu_max_B  = 10^p(20); % Max. specific growth rateof Bacteria B (d^-1)
k_HY_B    = 10^p(21); % HY growth substrate affinity coeff Bacteria B (g soil (mg C d)^-1)
k_NI_B    = 10^p(22); % NI growth substrate affinity coeff Bacteria B (g soil (mg C d)^-1)
a_Bac     = 10^p(23); % Specific deat rate of act. Bacteria B (d^-1)
a_Bin     = 10^p(24); % Specific deat rate of inact. Bacteria B (d^-1)
K_reac_B  = 10^p(25); % Coeff rate of reactivation Bacteria B (d^-1)
K_deac_B  = 10^p(26); % Coeff rate of deactivation Bacteria B (d^-1)
Y_HY_B    = 10^p(27); % Substrate HY uptake efficiency of Bacteria B (1)
Y_NI_B    = 10^p(28); % Substrate NI uptake efficiency of Bacteria B (1)

% BACTERIA C

CT_HY_C   = 10^p(29); % Threshold HY conc. for growth Bacteria C (mg C cm^-3 water)
mu_max_C  = 10^p(30); % Max. specific growth rateof Bacteria C (d^-1)
k_HY_C    = 10^p(31); % HY growth substrate affinity coeff Bacteria C (g soil (mg C d)^-1)
k_NE_C    = 10^p(32); % NE growth substrate affinity coeff Bacteria C (g soil (mg C d)^-1)
k_NI_C    = 10^p(33); % NI growth substrate affinity coeff Bacteria C (g soil (mg C d)^-1)
a_Cac     = 10^p(34); % Specific deat rate of act. Bacteria C (d^-1)
a_Cin     = 10^p(35); % Specific deat rate of inact. Bacteria C (d^-1)
K_reac_C  = 10^p(36); % Coeff rate of reactivation Bacteria C (d^-1)
K_deac_C  = 10^p(37); % Coeff rate of deactivation Bacteria C (d^-1)
Y_HY_C    = 10^p(38); % Substrate HY uptake efficiency of Bacteria C (1)
Y_NE_C    = 10^p(39); % Substrate NE uptake efficiency of Bacteria C (1)
Y_NI_C    = 10^p(40); % Substrate NI uptake efficiency of Bacteria C (1)

% BACTERIA D

CT_AT_D   = 10^p(41); % Threshold AT conc. for growth Bacteria D (mg C cm^-3 water)
mu_max_D  = 10^p(42); % Max. specific growth rateof Bacteria D (d^-1)
k_AT_D    = 10^p(43); % AT growth substrate affinity coeff Bacteria D (g soil (mg C d)^-1)
a_Dac     = 10^p(44); % Specific deat rate of act. Bacteria D (d^-1)
a_Din     = 10^p(45); % Specific deat rate of inact. Bacteria D (d^-1)
K_reac_D  = 10^p(46); % Coeff rate of reactivation Bacteria D (d^-1)
K_deac_D  = 10^p(47); % Coeff rate of deactivation Bacteria D (d^-1)
Y_AT_D    = 10^p(48); % Substrate AT uptake efficiency of Bacteria D (1)

% OTHER PARAMETERS

k_AT_HY   = 10^ p(49); % Max. rate of AT degradation to HY (mg C cm^-3 d^-1)
K_AT_HY   = 10^p(50); % Michaelis-Menten Const. AT-HY (mg C cm^-3)
f1        = p(51); % Fraction of C assimilation by Bacteria D (1)
f3        = p(53); % Fraction of AT used for DEA formation by bacteria D (1)
f4        = p(72); % Fraction of dead bacteria which goes to DOC in soils
k0        = 10^p(54); % Coeff of abiotic decomposition of AT to HY (d^-1)
d_CYA_CO2 = 10^p(55); % Coeff of degradation of CYA (d^-1)
Ki        = 10^p(56); % Coeff of inhibition of CYA degrad (mg NO3 g^-1)
c_NO3     = 10^p(57); % Conc of NO3 in soil (mg NO3 g^-1)

% SORPTION PARAMETERS

KF_AT     = p(58); % Freundlich coeff of AT sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_AT)
KF_HY     = p(59); % Freundlich coeff of HY sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_HY)
KF_NE     = p(60); % Freundlich coeff of NE sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_NE)
KF_NI     = p(61); % Freundlich coeff of NI sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_NI)
KF_DEA    = p(62); % Freundlich coeff of DEA sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_DEA)
KF_DIA    = p(63); % Freundlich coeff of DIA sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_DIA)
KF_CYA    = p(64); % Freundlich coeff of CYA sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_CYA)
nF_AT     = p(65); % Freundlich exponent of AT sorption isotherm (1)
nF_HY     = p(66); % Freundlich exponent of HY sorption isotherm (1)
nF_NE     = p(67); % Freundlich exponent of NE sorption isotherm (1)
nF_NI     = p(68); % Freundlich exponent of NI sorption isotherm (1)
nF_DEA    = p(69); % Freundlich exponent of DEA sorption isotherm (1)
nF_DIA    = p(70); % Freundlich exponent of DIA sorption isotherm (1)
nF_CYA    = p(71); % Freundlich exponent of CYA sorption isotherm (1)

%% VALUES OF CONSTANT %%

n         = q(16);    % Switch function parameter
th_V      = q(17);   % Average volumetric soil water content (cm^3 cm^-3)
rho_B     = q(18);    % Bulk density of soils (g cm^-3)
f_NE_CYA  = q(21);    % Stoichiometric formation of CYA from NE
f_NI_CYA  = q(22);    % Stoichiometric formation of CYA from NI
f_HY_NE   = q(25);    % Stoichiometric formation of NE from HY
f_HY_NI   = q(26);    % Stoichiometric formation of NI from HY

%% BIOKINETIC FUNCTIONS %%

% switch-functions Bacteria B - D:
tau_B = 1/(exp((CT_AT_B-(CL_AT+CL_HY+CL_NI))/(n*CT_AT_B))+1);
tau_C = 1/(exp((CT_HY_C-(CL_HY+CL_NI+CL_NE))/(n*CT_HY_C))+1);
tau_D = 1/(exp((CT_AT_D-CL_AT)/(n*CT_AT_D))+1);

% Degrad. Coeff of Atrazine by Bacteria B (mg C d^-1 cm^-3 water):
K_dAT_HY = (k_AT_HY*CL_AT)/(K_AT_HY+CL_AT);

% Substrate (HY, NI) dependante growth rate for Bacteria B (d^-1):
mu_B_HY = (mu_max_B*CL_HY*k_HY_B)/(mu_max_B + CL_HY*k_HY_B ...
    + CL_NI*k_NI_B);
mu_B_NI = (mu_max_B*CL_NI*k_NI_B)/(mu_max_B + CL_HY*k_HY_B ...
    + CL_NI*k_NI_B);

% Substrate (HY, NE, NI) dependante growth rate for Bacteria C (d^-1):
mu_C_HY = (mu_max_C*CL_HY*k_HY_C)/(mu_max_C + CL_HY*k_HY_C ...
    + CL_NE*k_NE_C + CL_NI*k_NI_C);
mu_C_NE = (mu_max_C*CL_NE*k_NE_C)/(mu_max_C + CL_HY*k_HY_C ...
    + CL_NE*k_NE_C + CL_NI*k_NI_C);
mu_C_NI = (mu_max_C*CL_NI*k_NI_C)/(mu_max_C + CL_HY*k_HY_C ...
    + CL_NE*k_NE_C + CL_NI*k_NI_C);

% Substrate (AT) dependante growth rate for Bacteria D (d^-1):
mu_D_AT = (mu_max_D*CL_AT*k_AT_D)/(mu_max_D + CL_AT*k_AT_D);

%% STOCKS %%

vf_ = zeros(15,1);

% Active degrading bacteria B (mg C g^-1 soil) Bac_B:
vf_(1) = Bac_B*f_HY_NE*mu_B_HY + Bac_B*mu_B_NI*f_NI_CYA ...
      - (1 - tau_B)*K_deac_B*Bac_B + tau_B*K_reac_B*Bin_B - Bac_B*a_Bac;

% Inactive degrad. bacteria B (mg C g^-1 soil) Bin_B
vf_(2) = (1 - tau_B)*K_deac_B*Bac_B - tau_B*K_reac_B*Bin_B - Bin_B*a_Bin;

% Active degrading bacteria C (mg C g^-1 soil) Bac_C:
vf_(3) = Bac_C*f_HY_NI*mu_C_HY + Bac_C*mu_C_NE*f_NE_CYA ...
      + Bac_C*mu_C_NI*f_NI_CYA - (1 - tau_C)*K_deac_C*Bac_C ...
      + tau_C*K_reac_C*Bin_C - Bac_C*a_Cac;

% Inactive degrad. bacteria C (mg C g^-1 soil) Bin_C:
vf_(4) = (1 - tau_C)*K_deac_C*Bac_C - tau_C*K_reac_C*Bin_C - Bin_C*a_Cin;

% Active degrading bacteria D (mg C g^-1 soil) Bac_D:
vf_(5) = Bac_D*f1*mu_D_AT - (1 - tau_D)*K_deac_D*Bac_D ...
      + tau_D*K_reac_D*Bin_D - Bac_D*a_Dac;

% Inactive degrad. bacteria D (mg C g^-1 soil) Bin_D:
vf_(6) = (1 - tau_D)*K_deac_D*Bac_D - tau_D*K_reac_D*Bin_D - Bin_D*a_Din;

% Atrazine conc in Solution (mg C cm^-3):
vf_(7) = (- (K_dAT_HY*Bac_B)- k0*(CL_AT) ...
      + (-(Bac_D/th_V)*(rho_B)*mu_D_AT*((f1/Y_AT_D) + (1 - f1)))) ...
      / (1 + (rho_B/th_V)*KF_AT * nF_AT*(CL_AT^(nF_AT-1)));

% Hydroxyatrazine conc in Solution (mg C cm^-3)
vf_(8) = (+ K_dAT_HY*Bac_B + k0*(CL_AT)...
      - ((Bac_B/th_V)*(rho_B)*mu_B_HY*((f_HY_NE/Y_HY_B)+(1-f_HY_NE)))...
      - ((Bac_C/th_V)*(rho_B)*mu_C_HY*((f_HY_NI/Y_HY_C)+(1-f_HY_NI))))...
      / ( 1 + (rho_B/th_V)*KF_HY*nF_HY*(CL_HY^(nF_HY-1)));

% N-Ethylammelide conc in Solution (mg C cm^-3)
vf_(9) = ((Bac_B/th_V)*(rho_B)*mu_B_HY*(1 - f_HY_NE)...
      - ((Bac_C/th_V)*(rho_B)*mu_C_NE*((f_NE_CYA/Y_NE_C)+(1-f_NE_CYA))))...
      / ( 1 + (rho_B/th_V)*KF_NE*nF_NE*(CL_NE^(nF_NE-1)) );
    
% N-Isopropylammelide conc in Solution (mg C cm^-3)
vf_(10) = ((Bac_C/th_V)*(rho_B)*mu_C_HY*(1 - f_HY_NI)...
      - ((Bac_B/th_V)*(rho_B)*mu_B_NI*((f_NI_CYA/Y_NI_B)+(1-f_NI_CYA)))...
      - ((Bac_C/th_V)*(rho_B)*mu_C_NI*((f_NI_CYA/Y_NI_C)+(1-f_NI_CYA))))...
      / ( 1 + (rho_B/th_V)*KF_NI*nF_NI*(CL_NI^(nF_NI-1)));

% Deethylatrazine conc in Solution (mg C cm^-3)
vf_(11) = ((Bac_D/th_V)*(rho_B)*mu_D_AT*(1-f1)*f3)...
       / (1+(rho_B/th_V)*KF_DEA*nF_DEA*(CL_DEA^(nF_DEA-1)));
    
% Deisopropylatrazine conc in Solution (mg C cm^-3)
vf_(12) = ((Bac_D/th_V)*(rho_B)*mu_D_AT*(1-f1)*(1-f3))...
       / (1+(rho_B/th_V)*KF_DIA*nF_DIA*(CL_DIA^(nF_DIA-1)));
    
% Cyanuric Acid conc in Solution (mg C cm^-3)
vf_(13) = ((Bac_C/th_V)*(rho_B)*mu_C_NE*(1-f_NE_CYA)...   
       + (Bac_B/th_V)*(rho_B)*mu_B_NI*(1-f_NI_CYA)...
       + (Bac_C/th_V)*(rho_B)*mu_C_NI*(1-f_NI_CYA)...      
       - (CL_CA*d_CYA_CO2*(Ki/(c_NO3 + Ki))))...
       / (1+(rho_B/th_V)*KF_CYA*nF_CYA*(CL_CA^(nF_CYA-1)));
   
% CO2 (mg C g^-1 soil)
vf_(14) = (CL_CA*d_CYA_CO2*(Ki/(c_NO3 + Ki)))*(th_V/rho_B)...
    +  Bac_D*mu_D_AT*f1*((1 - Y_AT_D)/Y_AT_D)...
    +  Bac_C*f_NE_CYA*mu_C_NE*((1 - Y_NE_C)/Y_NE_C) ...
    +  Bac_B*f_NI_CYA*mu_B_NI*((1 - Y_NI_B)/Y_NI_B) ...
    +  Bac_C*f_NI_CYA*mu_C_NI*((1 - Y_NI_C)/Y_NI_C) ...
    +  Bac_B*mu_B_HY*f_HY_NE*((1 - Y_HY_B)/Y_HY_B) ...
    +  Bac_C*mu_C_HY*f_HY_NI*((1 - Y_HY_C)/Y_HY_C)...
    +  (Bac_B*a_Bac + Bin_B*a_Bin + Bac_C*a_Cac + Bin_C*a_Cin + Bac_D*a_Dac ...
    +  Bin_D*a_Din)*(1-f4);

% DOC pool (mg C g^(-1)  soil)
vf_(15) = (Bac_B*a_Bac + Bin_B*a_Bin + Bac_C*a_Cac + Bin_C*a_Cin + Bac_D*a_Dac ...
    +  Bin_D*a_Din)*(1-f4);

end
