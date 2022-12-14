clear all
close all
clc
%Initial conditions
rmax = 2;
rmin = 0.1;
dmax = 1;
dmin = 0.01;
alpha = 0.2;
alpha2 = 0.5;
mut = 0.05;
k = 10^5;
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
nzt = zeros(zsize, tmax);
nzt(:,1) = nz;
%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
R = repro(x, y.mu, rmin, rmax, alpha,z,nz);
S = starv(x,dmin, dmax, alpha2);
winsize = 10;
win = randi(100,[1,winsize]);
threshold = 0.1;
tdel = 1;
t = 0:200:1000; 
while mean(abs(diff(win))) > threshold  %checking for convergence
for   i = 1:length(t) %checking for convergence
    nznew = nz + (nz.* R* (1 - (sum(nz)/k)) - S.* nz)*tdel; %.*(max(0,R*(1-sum(nz)/k)))
    nz = nznew;
    nzt(:,i)=nznew;
y.mu = dot(nz/sum(nz),x);
ymu = y.mu;
ymutime(:,i) = ymu;
win(1:(winsize-1)) = win(2:winsize);
win(winsize) = sum(nz);

end
end
hold on
figure(1)
hold on
plot(x,nz,'Color', [0.9 0.9 0.9],'DisplayName','rmin = 0.1')
title('new distribution of traits')
xlabel('z')
ylabel('no of individuals')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
%Initial conditions
rmax = 2;
rmin = 0.5;
dmax = 1;
dmin = 0.01;
alpha = 0.2;
alpha2 = 0.5;
mut = 0.05;
k = 10^5;
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
nzt = zeros(zsize, tmax);
nzt(:,1) = nz;
%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
R = repro(x, y.mu, rmin, rmax, alpha,z,nz);
S = starv(x,dmin, dmax, alpha2);
winsize = 10;
win = randi(100,[1,winsize]);
threshold = 0.1;
tdel = 1;
t = 0:200:1000; 
while mean(abs(diff(win))) > threshold  %checking for convergence
for    i = 1:length(t)  %checking for convergence
    nznew = nz + (nz.* R* (1 - (sum(nz)/k)) - S.* nz)*tdel; %.*(max(0,R*(1-sum(nz)/k)))
    nz = nznew;
    nzt(:,i)=nznew;

y.mu = dot(nz/sum(nz),x);
ymu = y.mu;
ymutime(:,i) = ymu;
win(1:(winsize-1)) = win(2:winsize);
win(winsize) = sum(nz);

end
end
figure(1)
plot(x,nz,'Color', [0.8 0.8 0.8],'DisplayName','rmin = 0.5')
title('new distribution of traits')
xlabel('z')
ylabel('no of individuals')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
%Initial conditions
rmax = 2;
rmin = 0.9;
dmax = 1;
dmin = 0.01;
alpha = 0.2;
alpha2 = 0.5;
mut = 0.05;
k = 10^5;
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
nzt = zeros(zsize, tmax);
nzt(:,1) = nz;
%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
R = repro(x, y.mu, rmin, rmax, alpha,z,nz);
S = starv(x,dmin, dmax, alpha2);
winsize = 10;
win = randi(100,[1,winsize]);
threshold = 0.1;
tdel = 1;
t = 0:200:1000
while mean(abs(diff(win))) > threshold  %checking for convergence
for   i = 1:length(t) %checking for convergence
    nznew = nz + (nz.* R* (1 - (sum(nz)/k)) - S.* nz)*tdel; %.*(max(0,R*(1-sum(nz)/k)))
    nz = nznew;
    nzt(:,i)=nznew;

y.mu = dot(nz/sum(nz),x);
ymu = y.mu;
ymutime(:,i) = ymu;
win(1:(winsize-1)) = win(2:winsize);
win(winsize) = sum(nz);

end
end
 figure(1)
plot(x,nz,'Color', [0.6 0.6 0.6],'DisplayName','rmin = 0.9')
title('new distribution of traits')
xlabel('z')
ylabel('no of individuals')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
%Initial conditions
rmax = 2;
rmin = 1.5;
dmax = 1;
dmin = 0.01;
alpha = 0.2;
alpha2 = 0.5;
mut = 0.05;
k = 10^5;
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
nzt = zeros(zsize, tmax);
nzt(:,1) = nz;
%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
R = repro(x, y.mu, rmin, rmax, alpha,z,nz);
S = starv(x,dmin, dmax, alpha2);
winsize = 10;
win = randi(100,[1,winsize]);
threshold = 0.1;
tdel = 1;
t = 0:200:1000
while mean(abs(diff(win))) > threshold  %checking for convergence
for    i = 1:length(t)  %checking for convergence
    nznew = nz + (nz.* R* (1 - (sum(nz)/k)) - S.* nz)*tdel; %.*(max(0,R*(1-sum(nz)/k)))
    nz = nznew;
    nzt(:,i)=nznew;

y.mu = dot(nz/sum(nz),x);
ymu = y.mu;
ymutime(:,i) = ymu;
win(1:(winsize-1)) = win(2:winsize);
win(winsize) = sum(nz);

end
end
figure(1)
plot(x,nz,'Color', [0.4 0.4 0.4],'DisplayName','rmin = 1.5')
title('new distribution of traits')
xlabel('z')
ylabel('no of individuals')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
%Initial conditions
rmax = 2;
rmin = 2;
dmax = 1;
dmin = 0.01;
alpha = 0.2;
alpha2 = 0.5;
mut = 0.05;
k = 10^5;
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
nzt = zeros(zsize, tmax);
nzt(:,1) = nz;
%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
R = repro(x, y.mu, rmin, rmax, alpha,z,nz);
S = starv(x,dmin, dmax, alpha2);
winsize = 10;
win = randi(100,[1,winsize]);
threshold = 0.1;
tdel = 1;
t = 0:200:1000
while mean(abs(diff(win))) > threshold  %checking for convergence
for  i = 1:length(t)  %checking for convergence
    nznew = nz + (nz.* R* (1 - (sum(nz)/k)) - S.* nz)*tdel; %.*(max(0,R*(1-sum(nz)/k)))
    nz = nznew;
    nzt(:,i)=nznew;
y.mu = dot(nz/sum(nz),x);
ymu = y.mu;
ymutime(:,i) = ymu;
win(1:(winsize-1)) = win(2:winsize);
win(winsize) = sum(nz);

end
end
 figure(1)
plot(x,nz, 'k-*','DisplayName','rmin = 2')
title('new distribution of traits')
xlabel('z')
ylabel('no of individuals')

ylim([0 85000])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold off
legend('show')
breakyaxis([7000 70000]);