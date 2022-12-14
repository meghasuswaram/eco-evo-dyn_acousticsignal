
function [RSzbar] = reprostarvzbar(rmin, rmax, dmax, dmin, ini_pop , alphap, alpha2, k,tmax,sigmag)

zsize = 100; 
x = (1: 1: zsize); % Creating a vector from 0 to 100
y = makedist('Normal', 50, 10); %Making a normal distrubution with mean and std dev 
z = pdf(y,x); % Making a probability distribution 
zbar = y.mu;

nz = ini_pop*z;  %Putting 1000 individuals across the pdf

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
nzt = zeros(zsize, tmax); % initializing a population dataframe across time to be filled 
nzt(:,1) = nz; % the first column for t=1

  
    
for   t = 1:tmax %for t = 1 to tmax
if sum(nz) > 0.0001 % run simulation if total pop size is greater than 0 
    R = repro(x, y.mu, rmin, rmax, alphap,nz); %Reproductive function
    S = starv(x,dmin, dmax, alpha2); %Starvation function
    fitness = (R*(1-(sum(nz)/k)))-S;
    nznew =  nz.*fitness;
    %nznew = ((nz.* (R* (1 - (sum(nz)/k)))) - (S.* nz))*tdel; % New population until carrying capacity
    nz = nznew ;% replace old population with new population
    nz(nz < 0) = 0; % if individuals in population falls below zero, make it zero 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%introdue variance in next generation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nzmut = zeros(1,zsize); % 
    for i = 1:zsize % for iteration from 1 to 100 of syllable rate
       indz = nz(i); % picking number of individuals at each trait
       mutdist = pdf(makedist('Normal', i, sigmag),x); %creating normal pdf at each trait
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
   

end
end

RSzbar = zbartime;
end

