close all;
clear all;
file1='H:/M4_out_210424_082240.h5';
file2='H:/M4_out_210424_152014.h5';
file3='H:/M4_out_210428_161139.h5';


% read parameter bounds
pb=csvread('par_bounds_M4.csv');

parmin =	pb(:,1);
parmax =	pb(:,2);

h5disp(file2)
n_fin1=h5read(file1,'/n_fin');
x_fin1=h5read(file1,'/x_fin');
n_fin2=h5read(file2,'/n_fin');
x_fin2=h5read(file2,'/x_fin');
n_fin3=h5read(file3,'/n_fin');
x_fin3=h5read(file3,'/x_fin');

C=7;
t_points=401; % number of timepoints rejected from model function
n_mo=15; % number of model outputs
n_par=72; % number of parameters

n1=sum(n_fin1,2);
n2=sum(n_fin2,2);
n3=sum(n_fin3,2);

n_sel1=find(n1==C);
n_sel2=find(n2==C);
n_sel3=find(n3==C);


tol=1e-12;
[x_sel1,IA,IC]=uniquetol(x_fin1(n_sel1,:),tol,'Byrows',true);
[x_sel2,IA,IC]=uniquetol(x_fin2(n_sel2,:),tol,'Byrows',true);
[x_sel3,IA,IC]=uniquetol(x_fin3(n_sel3,:),tol,'Byrows',true);


n_p=size(x_sel1,2); % number of parameters
n_s=size(x_sel1,1); % number of sets

% K-S test
h=ones(n_p,1);
for iter_p=1:ceil(n_p/2)
h(iter_p) = kstest2(x_sel2(:,iter_p),x_sel3(:,iter_p));
end


figure('Name','Accepted parameters 1')
for iter_p=1:ceil(n_p/2)
    prior=unifrnd(parmin(iter_p),parmax(iter_p),1,n_s);
    subplot(ceil(n_p/2/5),5,iter_p)
    h1=histogram(prior);
    hold on;
    h2=histogram(x_sel1(:,iter_p));
    hold on;
    h3=histogram(x_sel2(:,iter_p));
    hold on;
    h4=histogram(x_sel3(:,iter_p));
    h1.Normalization = 'probability';
    h1.BinWidth = abs(parmax(iter_p)-parmin(iter_p))/25;
    h2.Normalization = 'probability';
    h2.BinWidth = h1.BinWidth;
    h3.Normalization = 'probability';
    h3.BinWidth = h1.BinWidth;
    h4.Normalization = 'probability';
    h4.BinWidth = h1.BinWidth; 
    title(num2str(iter_p))
    hold off;
end
legend({'prior','run1','run2','run3'},'Location','eastoutside');

figure('Name','Accepted parameters 2')
for iter_p=n_p/2+1:ceil(n_p)
    prior=unifrnd(parmin(iter_p),parmax(iter_p),1,n_s);
    subplot(ceil(n_p/2/5),5,iter_p-n_p/2)
    h1=histogram(prior);
    hold on;
    h2=histogram(x_sel1(:,iter_p));
    hold on;
    h3=histogram(x_sel2(:,iter_p));
    hold on;
    h4=histogram(x_sel3(:,iter_p));
    h1.Normalization = 'probability';
    h1.BinWidth = abs(parmax(iter_p)-parmin(iter_p))/25;
    h2.Normalization = 'probability';
    h2.BinWidth = h1.BinWidth;
    h3.Normalization = 'probability';
    h3.BinWidth = h1.BinWidth; 
    h4.Normalization = 'probability';
    h4.BinWidth = h1.BinWidth; 
    title(num2str(iter_p))
    hold off;
end
legend({'prior','run1','run2','run3'},'Location','eastoutside');




% read only selected data sets
fx_sel=zeros(n_mod,t_points,n_mo);
%for i_sel=1:n_mod
%i_sel
%   fx_sel(i_sel,:,:)=h5read(file1,'/fx_fin',[idx(i_sel) 1 1],[1 t_points n_mo]);
%end
fx_all=h5read(file1,'/fx_fin');
% randomly sample n_mod model runs
n_mod=100;
idx=datasample(n_sel1,n_mod,'Replace',false);
idx_compounds = [8,9,11,12,14];
names_compounds= {'Atrazine','Hydroxyatratzin','DIA','DEA','CO2'};
fx_sel=zeros(n_mod,t_points,n_mo);
fx_sel=fx_all(idx,:,:);

figure('Name','Compounds')
for i_f=1:length(idx_compounds)
   subplot(ceil(length(idx_compounds)/2),2,i_f) 
   for i_s=1:n_mod
        fx=squeeze(fx_sel(i_s,:,:));
        plot(fx(:,1),fx(:,idx_compounds(i_f)));
        hold on;
   end
   hold off;
   title(names_compounds{i_f});
   xlabel('Time (d)');
   ylabel('concentration (mg L^-1)');
    
end

% % plot parallel oordinate plot for randomly selected parameter sets
% figure('Name','Parallel coordinate plot')
% n_mod=10;
% x_p=datasample(x_sel,n_mod,'Replace',false);
% n_p=size(x_p,2);
% for iter_p=1:ceil(n_p/15)
%      subplot(ceil(n_p/15),1,iter_p)
%      idx_par = ((iter_p-1)*15+1):(min(iter_p*15,n_p));
%      parallelplot(x_p(:,idx_par));
%      title (num2str(idx_par));
% end