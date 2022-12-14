close all 
clear all
clc
%Initial conditions
dmax = 1;
dmin = 0.01;
alphap = 0.2;
alpha2 = 0.5;
mut = 0.05;
k = 10^5;
rmini = [0.1:0.1:2];
rmaxi = [0.1:0.1:2];
rminout2 = zeros(20,20);
rmaxout = zeros(20,20);
rmaxout2 = zeros(20,20);
for i = 1: length(rmini)
    rmin = rmini(i)
for j = 1: length(rmaxi)
    rmax = rmaxi(j)
% Initital Distribution of traits
zsize = 100;
tmax = 1000;
x = (1: 1: zsize);
y = makedist('Normal', 50, 10);
z = pdf(y,x);
% calculate zbar (this will inform R(zbar)) 
zbar = y.mu;
% N(z) i.e what is the no. of individuals with that trait (z)?
nz = 1000*z;
figure(1)
plot(x,nz)
title('initial distribution of traits')
xlabel('z')
ylabel('no of individuals')
nzt = zeros(zsize, tmax);
nzt(:,1) = nz;
%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
R = repro(x, y.mu, rmin, rmax, alphap,z,nz);
S = starv(x,dmin, dmax, alpha2);

winsize = 10;
win = randi(100,[1,winsize]);
threshold = 0.1;
tdel = 1;
while mean(abs(diff(win))) > threshold%checking for convergence
for   t = 2:1000 %checking for convergence
    nznew = nz + (nz.* R* (1 - (sum(nz)/k)) - S.* nz)*tdel; %.*(max(0,R*(1-sum(nz)/k)))
    nz = nznew;
    nzt(:,t)=nznew;

y.mu = dot(nz/sum(nz),x);
ymu = y.mu;
ymutime(:,t) = ymu;
win(1:(winsize-1)) = win(2:winsize);
win(winsize) = sum(nz);
% varp = var(nz/sum(nz),x);
% varptime(:,t) = varp;
end
end
stepinfo(nz);
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
rminout(j,i) = S1;%N*
rminout2(j,i) = S2;%zbar*
end
end
figure(2)
contourf(rmini-dmin, rmaxi-dmax, rminout, 'ShowText','on')
xlabel('rmin-dmin')
ylabel('rmax-dmax')
title('N*')

figure(3)
contourf(rmini-dmin, rmaxi-dmax, rminout2, 'ShowText','on')
colorbar
xlabel('rmin-dmin')
ylabel('rmax-dmax')
title('zbar*')


figure(6)
contourf(rmini, rmaxi, rminout2)
colorbar
xlabel({'rmin';'Increase in reproductive reward with low z (as rmin gets closer to rmax) ==>'})
ylabel({'rmax';'Increase in reproductive incentive ==>'})
title({'Mean signal Steady State';'Zbar*'})

figure(7)
contourf(rmini, rmaxi, rminout,'ShowText','on')
colormap(summer)
colorbar
xlabel({'rmin';'Increase in reproductive reward with low z (as rmin gets closer to rmax) ==>'})
ylabel({'rmax';'Increase in reproductive incentive ==>'})
title({'Population Steady State';'N*'})





rminout2(rminout2 > 1) = 0 ;
figure(8)
h = heatmap(rmini, rmaxi, rminout2);
h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
colorbar
xlabel({'rmin';'Increase in reproductive reward with low z -->'})
ylabel({'rmax';'Increase in reproductive incentive ==>'})
title('Bimodality')


