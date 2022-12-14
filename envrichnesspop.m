%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Mortality percentage 
%100 inidviduals in a population with different starting energy reserves
%Calculate Survival probability 
%Calculate Mortality probability 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
%Initial conditions
b = 30; %resting metabolic energy
s = 50; %signal strength in energy units
x1 = 1:10:1000; % different initial reserves at beginning for hundred different individuals
tmax = 100;


for i = 1: length(x1)% for different starting reserves consumptions 
    x = x1(i);
for t = 1:tmax % run for so many timesteps
    %environmental richness distribution
g = random('Normal',100,10); %randomly pick environmental gain from a normally varying environment
x = x - b - s + g;
xnew = x; 
xnewt(:,t)=xnew; %save all fueltank values for different time steps
if x > 10000 % if fuel tank reach a max of 10000 energy units, stop, cant eat more
    x = 10000;
end
if x <= 1000 %if fuel tank deplenished, stop. Dead!
    xnew = 0; 
xnewt(:,t)=xnew; %save all fueltank values for different time steps
end
end
probsurv(i) = xnewt(end); % calculate and store prob of surv by last value at last time step (whether tank full or zero) for different values of S
end
%Calculating mortality percentage below
surv = sum(probsurv > 2000); %no. of individuals surviving among the 100
survperc = (surv/length(x1))*100; % percentage of individuals surviving among the 100
mortalitypercentage = (100 - survperc) % Mortality percentage among the population

%Plotting surviving and dead individuals 
figure
scatter(x1, probsurv) % plot (prob of surv for different individuals with different starting energy reserves)
title('Surviving Individuals')
xlabel('Different initial - energy reserve tanks')
ylabel('Survival')
yline(2000,'-','Survival threshold'); % draw a line for survival threshold at energy tank 
txt = {'Dead individuals'};
text(765,1900,txt);
txt1 = ['(','Population Mortality % =', num2str(mortalitypercentage),')'];
text(10,2900,txt1);


