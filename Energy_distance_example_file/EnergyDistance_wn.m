function [D,temp2,Eij] = EnergyDistance_wn(data,n_mc,n_types,n_time,n_mod);

% n_mod         % # models
% n_types       % # observation types
% n_time        % # of observations along time
% n_mc          % # of realizations per model


tic

% computing fragments "Eij" of D^2

% figure which terms Eij are actually needed
% this is copied from the commands
% that compute D(i,j) based on Eij(i,j)
% and use only certain combinations of i,j
usecount=zeros(n_mod,n_mod);
for i=1:n_mod
    for j=i:n_mod
        usecount(i,j)=1;
        usecount(i,i)=1;
        usecount(j,j)=1;
    end
end

Eij = zeros(n_mod,n_mod);
for i=1:n_mod % enumerate all models (rows)
    for j=1:n_mod % enumerate all models, upper triangular (columns)
        if usecount(i,j)==1
            XY_SSE = zeros(n_mc,n_mc); % sumatoria, average (matrix)
            for k=1:n_types
                for l=1:n_time
                    % extract current "X" and "Y" (i,j) for data "k,l"
                    current_Xkl   = data(:,k,l,i); % vector realization
                    current_Ykl   = data(:,k,l,j);
                    % compute all SE between X and Y and aggregate to SSE
                    XY_SSE    = XY_SSE+(current_Xkl-current_Ykl').^2;
                end
            end
        end
        % take square root and average to get RMSE -> E_ij
        Eij(i,j)=mean(sqrt(XY_SSE(:))); % realization 
    end
end

% compute Dsquare and D matrices
D2 = zeros(n_mod,n_mod);
for i=1:n_mod
  for j=i:n_mod
    D2(i,j)=2*Eij(i,j)-Eij(i,i)-Eij(j,j); % 
  end
end
% copy symmetric elements
% temp2 = tril(E_dist2(:,:,d).',-1) + triu(E_dist2(:,:,d));   %mirror matrix to have all components of the matrix

for i=1:n_mod
  for j=1:i-1
    D2(i,j)=D2(j,i);
  end
end
temp2 = sum(D2(:)<0); % 
D  = sqrt(D2);
toc