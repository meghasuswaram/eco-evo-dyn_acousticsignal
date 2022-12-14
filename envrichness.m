close all
clear all
clc
%Initial conditions
b = 30; %resting metabolic energy
s1 = 50:1:80; %signal strength in energy units
x = 500; % energy tank at beginning 
tmax = 1000; %Max no. of time steps

for i = 1: length(s1)% for different signal consumptions 
    s = s1(i);
for t = 2:tmax % run for so many timesteps
%environmental richness distribution
g = random('Normal',100,10); %randomly pick environmental gain from a normally varying environment
x = x - b - s + g;
xnew = x; 
xnewt(:,t)=xnew; %save all fueltank values for different time steps
if x > 10000 % if fuel tank reach a max of 10000 energy units, stop, cant eat more
    x = 10000;
end
if x <= 500 %if fuel tank deplenished, stop. Dead!
    xnew = 0; 
xnewt(:,t)=xnew; %save all fueltank values for different time steps
end
end
probsurv(i) = xnewt(end); % calculate and store prob of surv by last value at last time step (whether tank full or zero) for different values of S
figure
plot(xnewt)
title('change in energy reserve for all elements in environmental variability distribution')
xlabel('time units')
ylabel('energy units')
end
figure
scatter(s1, probsurv) % plot (prob of surv vs different values of energy consumption)
title('change in survival for individuals')
xlabel('different signal consumption')
ylabel('energy units')
txt = {'Surviving individuals'};
text(80,10500,txt);
txt1 = {'Dead individuals'};
text(80,500,txt1)




 



