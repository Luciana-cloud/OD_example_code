function vf_= M6(t,x,p,q)

%% STOCKS

CL_AT      = x(1); % ATrazine in solution [mg C cm-3]
CL_HY      = x(2); % ATrazine in solution [mg C cm-3]
CL_DEA     = x(3); % ATrazine in solution [mg C cm-3]
CL_DIA     = x(4); % ATrazine in solution [mg C cm-3]
CL_CA      = x(5); % ATrazine in solution [mg C cm-3]

%% PARAMETERS %%

% FIRST ORDER DECAY PARAMETERS %

dAT_HY    = 10^p(14);          % Decay rate of AT to HY transformation [d^-1]
dAT_DD    = 10^p(15);          % Decay rate of AT to DD transformation [d^-1]
dHY_CYA   = 10^p(16);          % Decay rate of HY to CYA transformation [d^-1]
dDEA_CYA  = 10^p(17);          % Decay rate of DEA to CYA transformation [d^-1]
dDIA_CYA  = 10^p(18);         % Decay rate of DIA to CYA transformation [d^-1]
dCYA_CO2  = 10^p(19);         % Decay rate of CYA to CO2 transformation [d^-1]

% SORPTION PARAMETERS

KF_AT     = p(4); % Freundlich coeff of AT sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_AT)
KF_HY     = p(5); % Freundlich coeff of HY sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_HY)
KF_DEA    = p(6); % Freundlich coeff of DEA sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_DEA)
KF_DIA    = p(7); % Freundlich coeff of DIA sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_DIA)
KF_CYA    = p(8); % Freundlich coeff of CYA sorption isotherm (mg C g^-1 soil / (mg C cm^-3)^nF_CYA)
nF_AT     = p(9); % Freundlich exponent of AT sorption isotherm (1)
nF_HY     = p(10); % Freundlich exponent of HY sorption isotherm (1)
nF_DEA    = p(11); % Freundlich exponent of DEA sorption isotherm (1)
nF_DIA    = p(12); % Freundlich exponent of DIA sorption isotherm (1)
nF_CYA    = p(13); % Freundlich exponent of CYA sorption isotherm (1)

% OTHER PARAMETERS

Ki        = 10^p(2); % Coeff of inhibition of CYA degrad (mg NO3 g^-1)
c_NO3     = 10^p(3); % Conc of NO3 in soil (mg NO3 g^-1)
f3        = p(1); % Fraction of AT used for DEA formation by bacteria D (1)

%% VALUES OF CONSTANT %%

th_V      = q(17);   % Average volumetric soil water content (cm^3 cm^-3)
rho_B     = q(18);    % Bulk density of soils (g cm^-3)
f_HY_CYA  = q(20);    % Stoichiometric formation of CYA from HY
f_DEA_CYA = q(23);    % Stoichiometric formation of CYA from DEA
f_DIA_CYA = q(24);    % Stoichiometric formation of CYA from DIA

%% STOCKS %%

vf_ = zeros(6,1);

% Atrazine conc in Solution (mg C cm^-3):
vf_(1) = (- CL_AT * (dAT_HY+dAT_DD))...
      / (1 + (rho_B/th_V)*KF_AT * nF_AT*(CL_AT^(nF_AT-1)));

  % Hydroxyatrazine conc in Solution (mg C cm^-3)
vf_(2) = ((CL_AT*dAT_HY)-CL_HY*(dHY_CYA))...
      / ( 1 + (rho_B/th_V)*KF_HY*nF_HY*(CL_HY^(nF_HY-1)));

  % Deethylatrazine conc in Solution (mg C cm^-3)
vf_(3) = ((f3*CL_AT*dAT_DD)-CL_DEA*dDEA_CYA)...
       / (1+(rho_B/th_V)*KF_DEA*nF_DEA*(CL_DEA^(nF_DEA-1)));
    
% Deisopropylatrazine conc in Solution (mg C cm^-3)
vf_(4) = (((1-f3)*CL_AT*dAT_DD)-CL_DIA*dDIA_CYA)...
       / (1+(rho_B/th_V)*KF_DIA*nF_DIA*(CL_DIA^(nF_DIA-1)));

% Cyanuric Acid conc in Solution (mg C cm^-3)

vf_(5) = (CL_HY*dHY_CYA*f_HY_CYA... 
       +  CL_DEA*dDEA_CYA*f_DEA_CYA+CL_DIA*dDIA_CYA*f_DIA_CYA...
       - (CL_CA*dCYA_CO2*(Ki/(c_NO3 + Ki))))...
       / (1+(rho_B/th_V)*KF_CYA*nF_CYA*(CL_CA^(nF_CYA-1)));
   
% CO2 (mg C g^-1 soil)
vf_(6) = CL_HY*dHY_CYA*(1-f_HY_CYA)*(th_V/rho_B) ...
    + CL_DEA*dDEA_CYA*(1-f_DEA_CYA)*(th_V/rho_B)... *(th_V/rho_B)
    + CL_DIA*dDIA_CYA*(1-f_DIA_CYA)*(th_V/rho_B)... *(th_V/rho_B)
    + CL_CA*dCYA_CO2*((Ki/(c_NO3 + Ki)))*(th_V/rho_B);

end
