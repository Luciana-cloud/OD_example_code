function [n,n_idx] = check_process_constrains_chem(t,mo,extra)
%checks process constrains of chemostat-retentostat models
% input:
% t   time points
% mo  model outputs
% vector of extra parameters needed for setting the constrains
% output:
% n      total number of accepted constrains
% n_idx  index vector indicating which constrain was fullfilled

C_T=extra(1); % total added Atrazine concentration
n=0; % total number of accepted rules
n_idx=zeros(7,1); % array indicating which rule was fullfilled

% Rule 1: Half life AT time between 25 and 5 days.
if mo(end,7)<mo(1,7)/2 % check if half of initial AT was degraded at all
    t_half = min(t(mo(:,7)<=mo(1,7)/2));% determine approximately the half life
    if t_half>=5 && t_half<=25
        n=n+1;
        n_idx(n)=1;
    end
end

% Rule 2: Minimum AT concentration is at least 1 micrograms/l.
if mo(end,7)>=1e-5
    n=n+1;
    n_idx(n)=1;
end

% Rule 3: mineralization of initially added AT is between 20-80%
if (mo(end,13)/C_T>=0.2 && mo(end,13)/C_T<=0.8)
    n=n+1;
    n_idx(n)=1;
end

% Rule 4: The mimimum concentrationof AT+DEA+DIA+HY at the end is 1 micrograms/L.
if sum([mo(end,7),mo(end,8),mo(end,10),mo(end,11)])>1e-5
    n=n+1;
    n_idx(n)=1;
end

% Rule 5,6,7: half lifes of HY, DIA and DEA between 2-30 days
%prob=1;
for j=[8,10,11]
    prob=1;
    if max(mo(:,j))>1e-8 % check if metabolite was built at all
        if (mo(end,j)>=max(mo(:,j))/2) % accumulation until end?
            prob=0;
        else
            t_max =min(t(mo(:,j)==max(mo(:,j)))); % find time of max peak
            idx_max=t>t_max; % store index of times > t_max
            t_idx=t(idx_max); % store new time vector
            t_half = min(t_idx(mo(idx_max,j)<=max(mo(:,j))/2))-t_max;% determine approximately the half life
            if ~(t_half>=2 && t_half<=30)
                prob=0;
            end
        end
    end
    if prob==1
        n=n+1;
        n_idx(n)=1;
    end
    
end
end


