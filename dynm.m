function [dynamics] = dynm(traj, tstart, tend)

x = findpeaks(traj[tstart:tend])
y = %find valleys for the same trajectory
unique = %save the unique x and y together
z = unique

end

%across different values of rmax, run the above function save unique peaks
%and valleys

% plot the unique peaks and valleys across rmax

%flipped bifurcation (rmax* or rcritcal)

% plot where the bifurcation occurs (threshold) -  when fuzziness > some
% threshold 
% visually confirm it 



    
    
    