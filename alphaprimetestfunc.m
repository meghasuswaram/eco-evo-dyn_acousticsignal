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
tmax = 100;  % no of time steps to run 


alphapm = [0 0.2 0.4 0.8 1 2 4 8 10 20 40 60 80 100 200 400 600 800 1000];
alphapmoutzbar = zeros(1,19);
alphapmoutn = zeros(1,19);
alphaouttime = zeros(1,19);

for m = 1: length(alphapm)
    alphap = alphapm(m) 
RSn  = reprostarvpop( rmin, rmax, dmax, dmin, ini_pop , alphap, alpha2, k,tmax,sigmag);
RSzbar  = reprostarvzbar( rmin, rmax, dmax, dmin, ini_pop , alphap, alpha2, k,tmax,sigmag);
popvec = RSn;
zbarvec = RSzbar;
S1 = sum(popvec(:,end));
S2 = zbarvec(:,end);
alphapmoutzbar(m) = S2;
alphapmoutn(m) = S1;
alphaouttime(m) = length(zbarvec);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
semilogx(alphapm, alphapmoutzbar,'MarkerFaceColor',[0 0.447 0.741])
title('Syllable rate steady state for different alphaprime')
xlabel('alpha prime')
ylabel('zbar*')


figure(3)
semilogx(alphapm, alphapmoutn,'MarkerFaceColor',[0 0.447 0.741])
title('Population steady state for different alphaprime')
xlabel('alpha prime')
ylabel('N*')

figure(4)
scatter(alphapm, alphapmoutzbar,'filled')
title('Syllable rate steady state for different alphaprime')
xlabel('alpha prime')
ylabel('zbar*')


figure(5)
scatter(alphapm, alphapmoutn,'filled')
title('Population steady state for different alphaprime')
xlabel('alpha prime')
ylabel('N*')
