%% SOLUTION ATRAZINE - INITIAL CONDITIONS %%

% C_L = Atrazine in the solution phase, K_F*C_L^n_F = Atrazine in the sorbed phase
% C_T = Total Atrazine = C_L*(psy/roh) + K_F*C_L^n_F

function out = AT_init(C_L,C_T,K_F,n_F,th_V,rho_B)
 out = K_F*C_L^n_F + C_L*(th_V/rho_B) - C_T;
end