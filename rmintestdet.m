close all
clear all
clc
%Initial conditions
rmax = 2;
dmax = 1;
dmin = 0.01;
alpha = 0.2;
alpha2 = 0.5;
mut = 0.05;
k = 10^5;
rmini = 0.1:0.1:2.0;
rminout = zeros(1,20);
rminout2 = zeros(1,20);

for i = 1: length(rmini)
    rmin = rmini(i)
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
R = repro(x, y.mu, rmin, rmax, alpha,z,nz);
S = starv(x,dmin, dmax, alpha2);
% update N(z) + babyz - deadz
%nz = nz   +  { (1- mut-rate) *repro * nz + (mut-rate/2 * repro * nz(i+1)) +
%(mut-rate/2 * repro * nz(i-1))* ( 1 - sum(nz)/K )}   - starv**nz
%timeloop 
%zloop - for rdistribution of traits 
winsize = 10;
win = randi(100,[1,winsize]);
threshold = 0.1;
tdel = 1;
while mean(abs(diff(win))) > threshold%checking for convergence
for   t = 2:1000 %checking for convergence
    nznew = nz + (nz.* R* (1 - (sum(nz)/k)) - S.* nz)*tdel; %.*(max(0,R*(1-sum(nz)/k)))
    nz = nznew;
    nzt(:,t)=nznew;
%nz = nz + (((1-mut)*(R.*nz)) + (((mut/2)* (R.*nz(t+1)) )+ ((mut/2)*(R.*nz(t-1))))* (1 - (sum(nz)/k)))- S.*nz;
%if sum(nz) < 0  % if population falls below zero, then stop updating 
    %print("Negative crickets!")
   % break
%end
y.mu = dot(nz/sum(nz),x);
ymu = y.mu;
ymutime(:,t) = ymu;
win(1:(winsize-1)) = win(2:winsize);
win(winsize) = sum(nz);
% figure(2)
% plot(x,nz)
% title('new distribution of traits')
% xlabel('z')
% ylabel('no of individuals')
%pause(0.005)
end
end
stepinfo(nz);
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
rminout(i) = S1;%N*
rminout2(i) = S2;%zbar*
end

figure(3)
plot(sum(nzt(:,:)))
title('N vs time')
xlabel('time')
ylabel('N')
figure(4)
plot(ymutime)
title('Zbar vs time')
xlabel('time')
ylabel('Zbar')
figure(5)
scatter(rmini, rminout, 'filled')
title('N* vs rmin')
xlabel('rmin')
ylabel('N*')
figure(6)
scatter(rmini, rminout2, 'filled')

xlabel({'rmin';'Increase in reproductive reward with low z -->'})
ylabel({'Steady state of chorus mean syllable rate*'; 'zbar*'})
%endzloop
%endtimeloop