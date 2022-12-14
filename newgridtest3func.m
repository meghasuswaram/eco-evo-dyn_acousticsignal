close all
clear all
clc
rmax = 5;  %Max of reproductions with signaling
rmin = 1; %Min no. of reproductions without signaling 
dmax = 1;  %Max mortality due to starvation
dmin = 0.01; %Min mortality due to starvation 
ini_pop = 1000;
alpha2 = 0.5; %beta
k = 10^5; % Carrying capacity 
sigmag = 1; % genetic variation 
tmax = 500;  % no of time steps to run 
alphap = 200;

rmini = [1:1:10];
rmaxi = [1:1:10];
rminoutzbar = zeros(10,10);
rminoutn = zeros(10,10);
rminout2zbar = zeros(10,10);
rminout2n = zeros(10,10);

for i = 1: length(rmini)
    rmin = rmini(i)
for j = 1: length(rmaxi)
    rmax = rmaxi(j)

RSn  = reprostarvpop( rmin, rmax, dmax, dmin, ini_pop , alphap, alpha2, k,tmax,sigmag);
RSzbar  = reprostarvzbar( rmin, rmax, dmax, dmin, ini_pop , alphap, alpha2, k,tmax,sigmag);
popvec = RSn;
zbarvec = RSzbar;
% bimodal(j,i) = bimodalcheck;
stepinfo(RSn); 
S = stepinfo(RSn);
S1 = sum(popvec(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
rminoutn(j,i) = S1;%N*
rminoutzbar(j,i) = S2;%zbar*

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
contourf(rmini, rmaxi, rminoutzbar)
colorbar
xlabel({'rmin';'Increase in reproductive reward with low z (as rmin gets closer to rmax) ==>'})
ylabel({'rmax';'Increase in reproductive incentive ==>'})
title({'Mean signal Steady State';'Zbar*'})

figure(3)
contourf(rmini, rmaxi, rminoutn,'ShowText','on')
colormap(summer)
colorbar
xlabel({'rmin';'Increase in reproductive reward with low z (as rmin gets closer to rmax) ==>'})
ylabel({'rmax';'Increase in reproductive incentive ==>'})
title({'Population Steady State';'N*'})
