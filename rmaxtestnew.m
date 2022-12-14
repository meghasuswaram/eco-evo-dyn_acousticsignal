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
rmaxi = [1:0.1:10];
rminoutzbarmean = zeros(1,91);
rminoutzbarmode = zeros(1,91);
rminoutn = zeros(1,91);
rminout2zbar = zeros(1,91);
rminout2n = zeros(1,91);
temp_index = 1;
for j = 1: length(rmaxi)
    rmax = rmaxi(j)

%diff_rmax_rmin(1,temp_index) = rmax-rmin;   %storing rmax-rmin in a list at a specific index
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
for   t = 2:tmax %for t = 1 to tmax
    if sum(nz) > 0.0001
    R = repro(x, y.mu, rmin, rmax, alphap,nz); %Reproductive function
    S = starv(x,dmin, dmax, alpha2); %Starvation function 
    fitness = (R*(1-(sum(nz)/k)))-S;
    nznew = nz.*fitness; % New population until carrying capacity
    nz = nznew; % replace old population with new population
    nz(nz < 0) = 0; % if population falls below zero, make population zero 
    nzmut = zeros(1,zsize); % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%introdue variance in next generation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
peaksplot = findpeaks(nz);
noofpeaks(:,t) =length(peaksplot); %storing no of peaks at each time step 
bimodalcheck = any(noofpeaks == 2);


end
end

bimodal(j) = bimodalcheck;
stepinfo(nz); 
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
S3 = ymutime(:,end-1);
rminoutn(j) = S1;%N*
rminoutzbarmode(j) = S2;%zbar*
rminoutzbarmean(j) = S3;%zbar*
end


figure(1)
p1 = plot(rmaxi, rminoutzbarmean,'Color' ,  [1 0 0],'DisplayName','alphaprime = 1','Linewidth',1.5)
xlabel({'rmax';'Increase in reproductive incentive ==>'})
title({'Mean signal Steady State';'Zbar*'})
ylim([0 100])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

clear all
%Initial conditions
dmax = 1;
dmin = 0.01;
ini_pop = 1000;
alphap = 10;
alpha2 = 0.5;
k = 10^5;
sigmag = 1; % genetic variation 
tmax = 1000;  % no of time steps to run 
rmin = 1;
rmaxi = [1:0.1:10];
rminoutzbarmean = zeros(1,91);
rminoutzbarmode = zeros(1,91);
rminoutn = zeros(1,91);
rminout2zbar = zeros(1,91);
rminout2n = zeros(1,91);
temp_index = 1;
for j = 1: length(rmaxi)
    rmax = rmaxi(j)

%diff_rmax_rmin(1,temp_index) = rmax-rmin;   %storing rmax-rmin in a list at a specific index
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
for   t = 2:tmax %for t = 1 to tmax
    if sum(nz) > 0.0001
    R = repro(x, y.mu, rmin, rmax, alphap,nz); %Reproductive function
    S = starv(x,dmin, dmax, alpha2); %Starvation function 
    fitness = (R*(1-(sum(nz)/k)))-S;
    nznew = nz.*fitness; % New population until carrying capacity
    nz = nznew; % replace old population with new population
    nz(nz < 0) = 0; % if population falls below zero, make population zero 
    nzmut = zeros(1,zsize); % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%introdue variance in next generation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
peaksplot = findpeaks(nz);
noofpeaks(:,t) =length(peaksplot); %storing no of peaks at each time step 
bimodalcheck = any(noofpeaks == 2);


end
end

bimodal(j) = bimodalcheck;
stepinfo(nz); 
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
S3 = ymutime(:,end);
rminoutn(j) = S1;%N*
rminoutzbarmode(j) = S2;%zbar*
rminoutzbarmean(:,j) = S3;%zbar*
end

figure(1)
p2 = plot(rmaxi, rminoutzbarmean,'Color' ,  [0.8 0.2 0],'DisplayName','alphaprime = 10','Linewidth',1.5)
xlabel({'rmax';'Increase in reproductive incentive ==>'})
ylabel({'Mean signal Steady State';'Zbar*'})
ylim([0 100])
title({'Low alpha vs High alpha'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%Initial conditions
dmax = 1;
dmin = 0.01;
ini_pop = 1000;
alphap = 100;
alpha2 = 0.5;
k = 10^5;
sigmag = 1; % genetic variation 
tmax = 1000;  % no of time steps to run 
rmin = 1;
rmaxi = [1:0.1:10];
rminoutzbarmean = zeros(1,91);
rminoutzbarmode = zeros(1,91);
rminoutn = zeros(1,91);
rminout2zbar = zeros(1,91);
rminout2n = zeros(1,91);
temp_index = 1;
for j = 1: length(rmaxi)
    rmax = rmaxi(j)

%diff_rmax_rmin(1,temp_index) = rmax-rmin;   %storing rmax-rmin in a list at a specific index
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
for   t = 2:tmax %for t = 1 to tmax
    if sum(nz) > 0.0001
    R = repro(x, y.mu, rmin, rmax, alphap,nz); %Reproductive function
    S = starv(x,dmin, dmax, alpha2); %Starvation function 
    fitness = (R*(1-(sum(nz)/k)))-S;
    nznew = nz.*fitness; % New population until carrying capacity
    nz = nznew; % replace old population with new population
    nz(nz < 0) = 0; % if population falls below zero, make population zero 
    nzmut = zeros(1,zsize); % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%introdue variance in next generation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
peaksplot = findpeaks(nz);
noofpeaks(:,t) =length(peaksplot); %storing no of peaks at each time step 
bimodalcheck = any(noofpeaks == 2);


end
end

bimodal(j) = bimodalcheck;
stepinfo(nz); 
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
S3 = ymutime(:,end);
rminoutn(j) = S1;%N*
rminoutzbarmode(j) = S2;%zbar*
rminoutzbarmean(:,j) = S3;%zbar*
end

figure(1)
p3 = plot(rmaxi, rminoutzbarmean, 'Color' ,  [0.8 0.4 0],'DisplayName','alphaprime = 100','Linewidth',1.5)
xlabel({'rmax';'Increase in reproductive incentive ==>'})
ylabel({'Mean signal Steady State';'Zbar*'})
ylim([0 100])
title({'Low alpha vs High alpha'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%Initial conditions
dmax = 1;
dmin = 0.01;
ini_pop = 1000;
alphap = 200;
alpha2 = 0.5;
k = 10^5;
sigmag = 1; % genetic variation 
tmax = 1000;  % no of time steps to run 
rmin = 1;
rmaxi = [1:0.1:10];
rminoutzbarmean = zeros(1,91);
rminoutzbarmode = zeros(1,91);
rminoutn = zeros(1,91);
rminout2zbar = zeros(1,91);
rminout2n = zeros(1,91);
temp_index = 1;
for j = 1: length(rmaxi)
    rmax = rmaxi(j)

%diff_rmax_rmin(1,temp_index) = rmax-rmin;   %storing rmax-rmin in a list at a specific index
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
for   t = 2:tmax %for t = 1 to tmax
    if sum(nz) > 0.0001
    R = repro(x, y.mu, rmin, rmax, alphap,nz); %Reproductive function
    S = starv(x,dmin, dmax, alpha2); %Starvation function 
    fitness = (R*(1-(sum(nz)/k)))-S;
    nznew = nz.*fitness; % New population until carrying capacity
    nz = nznew; % replace old population with new population
    nz(nz < 0) = 0; % if population falls below zero, make population zero 
    nzmut = zeros(1,zsize); % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%introdue variance in next generation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for l = 1:zsize % for iteration from 1 to 100 of syllable rate
       indz = nz(l); % picking number of individuals at each trait
       mutdist = pdf(makedist('Normal', l, sigmag),x); %creating normal pdf at each trait
       nzmut = nzmut + (indz * mutdist);  % making a new mutation population       
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nz = nzmut; %replacing new population with individuals who have variation 
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
peaksplot = findpeaks(nz);
noofpeaks(:,t) =length(peaksplot); %storing no of peaks at each time step 
bimodalcheck = any(noofpeaks == 2);


end
end

bimodal(j) = bimodalcheck;
stepinfo(nz); 
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
S3 = ymutime(:,end);
rminoutn(j) = S1;%N*
rminoutzbarmode(j) = S2;%zbar*
rminoutzbarmean(:,j) = S3;%zbar*
end

figure(1)
p4 = plot(rmaxi, rminoutzbarmean,'Color' , [0.02 0.5 0.8],'DisplayName','alphaprime = 200','Linewidth',1.5)
xlabel({'rmax';'Increase in reproductive incentive ==>'})
ylabel({'Mean signal Steady State';'Zbar*'})
ylim([0 100])
title({'Low alpha vs High alpha'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%Initial conditions
dmax = 1;
dmin = 0.01;
ini_pop = 1000;
alphap = 500;
alpha2 = 0.5;
k = 10^5;
sigmag = 1; % genetic variation 
tmax = 1000;  % no of time steps to run 
rmin = 1;
rmaxi = [1:0.1:10];
rminoutzbarmean = zeros(1,91);
rminoutzbarmode = zeros(1,91);
rminoutn = zeros(1,91);
rminout2zbar = zeros(1,91);
rminout2n = zeros(1,91);
temp_index = 1;
for j = 1: length(rmaxi)
    rmax = rmaxi(j)

%diff_rmax_rmin(1,temp_index) = rmax-rmin;   %storing rmax-rmin in a list at a specific index
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
for   t = 2:tmax %for t = 1 to tmax
    if sum(nz) > 0.0001
    R = repro(x, y.mu, rmin, rmax, alphap,nz); %Reproductive function
    S = starv(x,dmin, dmax, alpha2); %Starvation function 
    fitness = (R*(1-(sum(nz)/k)))-S;
    nznew = nz.*fitness; % New population until carrying capacity
    nz = nznew; % replace old population with new population
    nz(nz < 0) = 0; % if population falls below zero, make population zero 
    nzmut = zeros(1,zsize); % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%introdue variance in next generation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
peaksplot = findpeaks(nz);
noofpeaks(:,t) =length(peaksplot); %storing no of peaks at each time step 
bimodalcheck = any(noofpeaks == 2);


end
end

bimodal(j) = bimodalcheck;
stepinfo(nz); 
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
S3 = ymutime(:,end);
rminoutn(j) = S1;%N*
rminoutzbarmode(j) = S2;%zbar*
rminoutzbarmean(:,j) = S3;%zbar*
end

figure(1)
p5 = plot(rmaxi, rminoutzbarmean,'Color' , [0.02 0.4 0.8],'DisplayName','alphaprime = 500','Linewidth',1.5)
xlabel({'rmax';'Increase in reproductive incentive ==>'})
ylabel({'Mean signal Steady State';'Zbar*'})
ylim([0 100])
title({'Low alpha vs High alpha'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
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
rmaxi = [1:0.1:10];
rminoutzbarmean = zeros(1,91);
rminoutzbarmode = zeros(1,91);
rminoutn = zeros(1,91);
rminout2zbar = zeros(1,91);
rminout2n = zeros(1,91);
temp_index = 1;
for j = 1: length(rmaxi)
    rmax = rmaxi(j)

%diff_rmax_rmin(1,temp_index) = rmax-rmin;   %storing rmax-rmin in a list at a specific index
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
for   t = 2:tmax %for t = 1 to tmax
    if sum(nz) > 0.0001
    R = repro(x, y.mu, rmin, rmax, alphap,nz); %Reproductive function
    S = starv(x,dmin, dmax, alpha2); %Starvation function 
    fitness = (R*(1-(sum(nz)/k)))-S;
    nznew = nz.*fitness; % New population until carrying capacity
    nz = nznew; % replace old population with new population
    nz(nz < 0) = 0; % if population falls below zero, make population zero 
    nzmut = zeros(1,zsize); % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%introdue variance in next generation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
peaksplot = findpeaks(nz);
noofpeaks(:,t) =length(peaksplot); %storing no of peaks at each time step 
bimodalcheck = any(noofpeaks == 2);


end
end

bimodal(j) = bimodalcheck;
stepinfo(nz); 
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
S3 = ymutime(:,end);
rminoutn(j) = S1;%N*
rminoutzbarmode(j) = S2;%zbar*
rminoutzbarmean(:,j) = S3;%zbar*
end

figure(1)
p6 = plot(rmaxi, rminoutzbarmean,'Color' ,[0 0 1],'DisplayName','alphaprime = 1000', 'Linewidth',1.5)
xlabel({'rmax';'Increase in reproductive incentive ==>'})
ylabel({'Mean signal Steady State';'Zbar*'})
ylim([0 100])

hold off
% legend('show')
