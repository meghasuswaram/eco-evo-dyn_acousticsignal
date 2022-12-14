close all
clear all
clc
%Initial conditions
rmax = 2;  %Max of reproductions with signaling
rmin = 0.5; %Min no. of reproductions without signaling 
dmax = 1;  %Max 
dmin = 0.01;
alphap = 0.2;
alpha2 = 0.5;
k = 10^5;
sigmag =1;

% Initital Distribution of traits
zsize = 100;
tmax = 1000;
x = (1: 1: zsize);
y = makedist('Normal', 50, 10);
mutdist = makedist('Normal',0, sigmag);
z = pdf(y,x);
 plot(x,z);
 

% calculate zbar (this will inform R(zbar)) 
zbar = y.mu;

% N(z) i.e what is the no. of individuals with that trait (z)?
nz = 1000*z;
figure(1)
plot(x,nz)
title('initial distribution of traits')
xlabel('z')
ylabel('Density distribution of population (indv/m^2)')

nzt = zeros(zsize, tmax);
nzt(:,1) = nz;

%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
R = repro(x, y.mu, rmin, rmax, alphap,z,nz);
S = starv(x,dmin, dmax, alpha2);
winsize = 10;
win = randi(100,[1,winsize]);
threshold = 0.1;
tdel = 1;

while mean(abs(diff(win))) > threshold  %checking for convergence
for   t = 1:100 %checking for convergence
    nznew = nz + (nz.* R* (1 - (sum(nz)/k)) - S.* nz)*tdel; 
    nznew1 = nznew - nz;
    nz = nznew1;
 win(1:(winsize-1)) = win(2:winsize);
 win(winsize) = sum(nz);
  
    nz(nz < 0) = 0;
%     nzmut = zeros(1,zsize)
%     %introdue variance in next generation 
%     for i = 1:zsize
%        indz = nz(i);
%       mutdist = pdf(makedist('Normal', i, sigmag),x);
%        
%        nzmut = nzmut + (indz * mutdist);
       
%        if indz > 0
%        draw = round(normrnd(0,sigmag,[1,indz]),0)
% 
%        mutz = max(min(i + (draw), zsize),1)
%        uniquez = unique(mutz)
%        for j = 1: length(uniquez)
%        nzmut(uniquez(j)) =  nzmut(uniquez(j))+ numel(find(mutz==uniquez(j)))
%        end
% %        end
%     end
%     nz = nzmut;
    nzt(:,t)=nz;
% if nz > 0
y.mu = dot(nz/sum(nz),x);
% else 
% y.mu = 0;
% end
ymu = y.mu;
ymutime(:,t) = ymu;


figure(2)
 plot(x,nz)
title('new distribution of traits')
xlabel('z')
ylabel('no of individuals')
% pause(0.005)
peaksplot = findpeaks(nz);
noofpeaks(:,t) =length(peaksplot);
% F(t) = getframe(gcf) ;
      drawnow
end
end
%  % create the video writer with 1 fps
%   writerObj = VideoWriter('myVideo.avi');
%   writerObj.FrameRate = 5;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for t=2:length(F)
%     % convert the image to a frame
%     frame = F(t) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

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
%stepinfo(nz)

%  noofpeaks(noofpeaks < 2) = 0;
%  noofpeaks(noofpeaks > 0) = 1;

%endzloop
%endtimeloop
