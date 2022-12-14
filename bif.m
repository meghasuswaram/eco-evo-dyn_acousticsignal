close all 
clear all
clc
%Initial conditions
dmax = 1;
dmin = 0.01;
ini_pop = 1000;
alphap = 1;
alpha2 = 0.5;
k = 10^5;
sigmag = 1; % genetic variation 
tmax = 1000;  % no of time steps to run 
rmin = 1;
rmaxi = [2:0.01:10];
rminoutzbar = zeros(801,801);
rminoutn = zeros(801,801);
rminout2zbar = zeros(801,801);
rminout2n = zeros(801,801);
%  peak1d = zeros(10,10);
% peak1 = zeros(100,10);
% valley1 = zeros(10,10);
% uniquex = zeros(10,10);
% uniquey = zeros(10,10);
% uniquez = zeros(10,10);


for j = 1: length(rmaxi)
    rmax = rmaxi(j)


% Initital Distribution of traits
zsize = 100;

x = (1: 1: zsize);
y = makedist('Normal', 50, 10);
z = pdf(y,x);
% calculate zbar (this will inform R(zbar)) 
zbar = y.mu;
% N(z) i.e what is the no. of individuals with that trait (z)?
nz = ini_pop*z;
% figure(1)
% plot(x,nz)
% title('initial distribution of traits')
% xlabel('z')
% ylabel('density')
nzt = zeros(zsize, tmax);
nzt(:,1) = nz;
%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
tdel = 1;
  %run loop only if population sum is greater than 0
for   t = 1:tmax %for t = 1 to tmax
    if sum(nz) > 0.0001
    R = repro(x, y.mu, rmin, rmax, alphap,nz); %Reproductive function
    S = starv(x,dmin, dmax, alpha2); %Starvation function 
    fitness = (R*(1-(sum(nz)/k)))-S;
    nznew = nz.*fitness; % New population until carrying capacity
    nz = nznew; % replace old population with new population
    nz(nz < 0) = 0; % if population falls below zero, make population zero 
%     nzmut = zeros(1,zsize); % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %introdue variance in next generation 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for l = 1:zsize % for iteration from 1 to 100 of syllable rate
%        indz = nz(l); % picking number of individuals at each trait
%        mutdist = pdf(makedist('Normal', l, sigmag),x); %creating normal pdf at each trait
%        nzmut = nzmut + (indz * mutdist);  % making a new mutation population       
%     end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     nz = nzmut; %replacing new population with individuals who have variation 
    nzt(:,t)=nz; % storing new population vector as dataframe at every time step
%
if nz >=0 & sum(nz) ~= 0
y.mu = dot(nz/sum(nz),x);
else if sum(nz) == 0
y.mu = 0;
    end
 end
ymu = y.mu;
ymutime(:,t) = ymu;


end
end


stepinfo(nz); 
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
rminoutn(:,j) = S1;%N*
rminoutzbar(:,j) = S2;%zbar*]
peak1a(:,j) = (sum(nzt(:,:)));
peak3a(:,j) = S.PeakTime;
peak1b(:,j) = peak1a(end-10:end);
if length(peak3a) > 10
peak3b{j} = peak3a(end-10:end);
else
    peak3b{j} = peak3a(end);
end
peak1c(:,j) = round(peak1b(:,j));
peak3c{j} = round(peak3b{j},2);
peak1d{j} = unique(peak1c(:,j));
peak3d{j} = unique(peak3c{j});

% peak1(:,:) = findpeaks(peak1b(:,j))
% valley1(:,:) = islocalmin(peak1a(end-2:end));
% uniquex(:,:) = unique(peak1);
% uniquey(:,:) = unique(valley1);

end

rows = cellfun(@numel,peak1d);
cols = size(peak1d,2);
peak1e = zeros(max(rows),cols);
for k = 1:cols
    peak1e(1:rows(k),k) = peak1d{k};
end
peak1e(peak1e==0) = nan;

s1 = scatter(rmaxi, peak1e,'r','.', 'DisplayName','alpha-prime = 1')
xlabel({'Increase in reproductive incentive ==>';'r_{max}'})
ylabel('Population size of the host')
title('Bifurcation')
hold on

clear all
clc

%Initial conditions
dmax = 1;
dmin = 0.01;
ini_pop = 1000;
alphap = 1000;
alpha2 = 0.5;
k = 10^5;
sigmag = 1; % genetic variation 
tmax = 1000;  % no of time steps to run 
rmin = 1;
rmaxi = [2:0.01:10];
rminoutzbar = zeros(801,801);
rminoutn = zeros(801,801);
rminout2zbar = zeros(801,801);
rminout2n = zeros(801,801);
%  peak1d = zeros(10,10);
% peak1 = zeros(100,10);
% valley1 = zeros(10,10);
% uniquex = zeros(10,10);
% uniquey = zeros(10,10);
% uniquez = zeros(10,10);


for j = 1: length(rmaxi)
    rmax = rmaxi(j)


% Initital Distribution of traits
zsize = 100;

x = (1: 1: zsize);
y = makedist('Normal', 50, 10);
z = pdf(y,x);
% calculate zbar (this will inform R(zbar)) 
zbar = y.mu;
% N(z) i.e what is the no. of individuals with that trait (z)?
nz = ini_pop*z;
% figure(1)
% plot(x,nz)
% title('initial distribution of traits')
% xlabel('z')
% ylabel('density')
nzt = zeros(zsize, tmax);
nzt(:,1) = nz;
%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
tdel = 1;
  %run loop only if population sum is greater than 0
for   t = 1:tmax %for t = 1 to tmax
    if sum(nz) > 0.0001
    R = repro(x, y.mu, rmin, rmax, alphap,nz); %Reproductive function
    S = starv(x,dmin, dmax, alpha2); %Starvation function 
    fitness = (R*(1-(sum(nz)/k)))-S;
    nznew = nz.*fitness; % New population until carrying capacity
    nz = nznew; % replace old population with new population
    nz(nz < 0) = 0; % if population falls below zero, make population zero 
%     nzmut = zeros(1,zsize); % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %introdue variance in next generation 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for l = 1:zsize % for iteration from 1 to 100 of syllable rate
%        indz = nz(l); % picking number of individuals at each trait
%        mutdist = pdf(makedist('Normal', l, sigmag),x); %creating normal pdf at each trait
%        nzmut = nzmut + (indz * mutdist);  % making a new mutation population       
%     end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     nz = nzmut; %replacing new population with individuals who have variation 
    nzt(:,t)=nz; % storing new population vector as dataframe at every time step
%
if nz >=0 & sum(nz) ~= 0
y.mu = dot(nz/sum(nz),x);
else if sum(nz) == 0
y.mu = 0;
    end
 end
ymu = y.mu;
ymutime(:,t) = ymu;


end
end


stepinfo(nz); 
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
rminoutn(:,j) = S1;%N*
rminoutzbar(:,j) = S2;%zbar*]
peak1a(:,j) = (sum(nzt(:,:)));
peak3a(:,j) = S.PeakTime;
peak1b(:,j) = peak1a(end-10:end);
if length(peak3a) > 10
peak3b{j} = peak3a(end-10:end);
else
    peak3b{j} = peak3a(end);
end
peak1c(:,j) = round(peak1b(:,j));
peak3c{j} = round(peak3b{j},2);
peak1d{j} = unique(peak1c(:,j));
peak3d{j} = unique(peak3c{j});

% peak1(:,:) = findpeaks(peak1b(:,j))
% valley1(:,:) = islocalmin(peak1a(end-2:end));
% uniquex(:,:) = unique(peak1);
% uniquey(:,:) = unique(valley1);

end

rows = cellfun(@numel,peak1d);
cols = size(peak1d,2);
peak1e = zeros(max(rows),cols);
for k = 1:cols
    peak1e(1:rows(k),k) = peak1d{k};
end
peak1e(peak1e==0) = nan;

s2 = scatter(rmaxi, peak1e,'b','.','DisplayName','alpha-prime = 1000')
xlabel({'Increase in reproductive incentive ==>';'r_{max}'})
ylabel('Population size of the host')
title('Bifurcation')


