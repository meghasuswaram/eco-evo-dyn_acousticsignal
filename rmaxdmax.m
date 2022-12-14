close all 
clear all
clc
%Initial conditions
dmaxi = [1:0.5:10];
dmin = 0.01;
ini_pop = 1000;
alphap = 10;
alpha2 = 0.5;
k = 10^5;
sigmag = 1; % genetic variation 
tmax = 1000;  % no of time steps to run 
rmin = 1;
rmaxi = [1:0.5:10];
rminoutzbarmean = zeros(19,19);
rminoutzbarmode = zeros(19,19);
rminoutn = zeros(19,19);
rminout2zbar = zeros(19,19);
rminout2n = zeros(19,19);
temp_index = 1;
for i = 1: length(dmaxi)
    dmax = dmaxi(i)
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
figure(1)
plot(x,nz)
title('initial distribution of traits')
xlabel('z')
ylabel('density')
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

bimodal(j,i) = bimodalcheck;
stepinfo(nz); 
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S2 = S.PeakTime; %zbar*
S3 = ymutime(:,end);
rminoutn(j,i) = S1;%N*
rminoutzbarmode(j,i) = S2;%zbar*
rminoutzbarmean(j,i) = S3;%zbar*
end
end



figure(6)
contourf(dmaxi, rmaxi, rminoutzbarmean)
colorbar
xlabel({'dmax';'Increase in cost of signaling ==>'})
ylabel({'delta_r';'Increase in incentive ==>'})
title({'Mean signal Steady State';'Zbar*'})

figure(7)
contourf(dmaxi, rmaxi, rminoutn,'ShowText','on')
colormap(summer)
colorbar
xlabel({'dmax';'Increase in cost of signaling ==>'})
ylabel({'delta_r';'Increase in incentive ==>'})
title({'Population Steady State';'N*'})

figure(8)
contourf(dmaxi, rmaxi, rminoutzbarmode)
colorbar('Ticks',[2,100],...
         'TickLabels',{'Silent','Conspicuous'})
xlabel({'dmax';'Increase in cost of signaling ==>'})
ylabel({'delta_r';'Increase in reproductive reward ==>'})
title({'Mode signal Steady State';'Zbar*'})

