close all 
clear all
clc
%Initial conditions
dmax = 1;
dmin = 0.01;
ini_pop = 1000;
alphap = 0.2*ini_pop;
alpha2 = 0.5;
k = 10^5;
sigmag = 1; % genetic variation 
tmax = 1000;  % no of time steps to run 
rmini = [1:0.5:10];
rmaxi = [1:0.5:10];
dyntype = zeros(19,19);



for i = 1: length(rmini)
    rmin = rmini(i)
for j = 1: length(rmaxi)
    rmax = rmaxi(j)
if rmax > rmin

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


 end
end
peak1a(j,i) = {sum(nzt(:,:))};

end
end
end
peak1aa = peak1a;
tf = cellfun('isempty',peak1aa);
peak1aa(tf) = {ones(1,1000)};
peak1b = cellfun(@(peak1aa) peak1aa(end-20:end),peak1aa,'UniformOutput',false);
peak1c = cellfun(@(peak1b) round(peak1b,-2),peak1b,'UniformOutput',false);
peak1d = cellfun(@(peak1c) unique(peak1c),peak1c,'UniformOutput',false);
peak1e = cellfun(@(peak1d) length(peak1d),peak1d,'UniformOutput',false);
peak1f = cell2mat(peak1e);

peak1f(peak1f==1) = 0;

peak1f(peak1f>=2 & peak1f <=7) = 1;

peak1f(peak1f>=6) = 2;



figure(2)
h = heatmap(rmini, rmaxi, peak1f);
h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
%colorbar
xlabel({'rmin';'Increase in reproductive reward with low z -->'})
ylabel({'rmax';'Increase in reproductive incentive ==>'})
title({'Population';'Stability = 0, Period doubling =1, Chaos = 2'})
grid off
