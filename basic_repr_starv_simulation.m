close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmax = 9.5;  %Max of reproductions with signaling
rmin = 1; %Min no. of reproductions without signaling 
dmax = 1;  %Max mortality due to starvation
dmin = 0.01; %Min mortality due to starvation 
ini_pop = 1000;
alphap = 100; % alphaprime
alpha2 = 0.5; %beta
k = 10^5; % Carrying capacity 
sigmag = 1; % genetic variation 
tmax = 1000;  % no of time steps to run 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up Initital Distribution of traits and popultaion vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zsize = 100; 
x = (1: 1: zsize); % Creating a vector from 0 to 100
y = makedist('Normal', 50, 10); %Making a normal distrubution with mean and std dev 
z = pdf(y,x); % Making a probability distribution 
%plot(x,z); % Plotting to check 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% calculate zbar (this will inform R(zbar)) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zbar = y.mu; %Setting zbar to y.mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N(z) i.e what is the no. of individuals with that trait (z)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nz = ini_pop*z; %Putting 1000 individuals across the pdf
% % figure(1)
%  plot(x,nz)
%  title('initial distribution of traits')
%  xlabel('z')
% ylabel('Density distribution of population (indv/m^2)')
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
nzt = zeros(zsize, tmax); % initializing a population dataframe across time to be filled 
nzt(:,1) = nz; % the first column for t=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tdel = 1; % if we need to slow down the function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if nz > 0   %run loop only if population sum is greater than 0
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
for   t = 2:tmax %for t = 1 to tmax
if sum(nz) > 0.0001 % run simulation if total pop size is greater than 0 
    R = repro(x, y.mu, rmin, rmax, alphap,nz); %Reproductive function
    S = starv(x,dmin, dmax, alpha2); %Starvation function
    fitness = (R*(1-(sum(nz)/k)))-S;
    nznew =  nz.*fitness;
    %nznew = ((nz.* (R* (1 - (sum(nz)/k)))) - (S.* nz))*tdel; % New population until carrying capacity
    nz = nznew ;% replace old population with new population
    nz(nz < 0) = 0; % if individuals in population falls below zero, make it zero 
    
    fitnesstime(:,t)= fitness;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%introdue variance in next generation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nzmut = zeros(1,zsize); % 
    for i = 1:zsize % for iteration from 1 to 100 of syllable rate
       indz = nz(i); % picking number of individuals at each trait
       mutdist = pdf(makedist('Normal', i, 3),x); %creating normal pdf at each trait
       nzmut = nzmut + (indz * mutdist);  % making a new mutation population       
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nz = nzmut; %replacing new population with individuals who have variation 
    nzt(:,t)=nz; % storing new population vector as dataframe at every time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculating zbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(nz) ~= 0
y.mu = dot(nz/sum(nz),x);
else if sum(nz) == 0
y.mu = 0;
    end
 end
%y.mu = sum (nz.*x)/sum(nz);
%y.mu = dot(nz/sum(nz),x); %dot product to find the mean 
zbar = y.mu; %replacing the old zbar with new pop zbar
zbartime(:,t) = zbar; % Storing zbar over all timesteps 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting new population after variation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
%     figure(3)
%     plot(x,fitness)
%     title('per-capita/Malthusian growth')
% xlabel('z - syllable rate')
% ylabel('individual fitness')
%  % ylim([2 4])
%      xlim([0 100])
%  % pause(0.5)
% figure(4) 
%  plot(x,nz)
% title('new distribution of traits')
% xlabel('z')
% ylabel('Density distribution of population (indv/m^2)')
%  %pause(.5)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Peak calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peaksplot = findpeaks(nz); %calculating the peak values
noofpeaks(:,t) =length(peaksplot); %storing no of peaks at each time step 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
plot(sum(nzt(:,:)))
% title('N vs time')
% xlabel('time')
% ylabel('N')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
plot(zbartime)
% title('Zbar vs time')
% xlabel('time')
% ylabel('Zbar')

% figure(7)
% plot(fitnesstime(:,(100:1000)))
% xlabel('z')
% ylabel('fitness')
% %ylim([0,2])




