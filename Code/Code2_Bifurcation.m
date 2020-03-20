echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: Code2_Bifurcation.m
% - Author: XYZ
% - Created date: March 20, 2020
% - Modified date: March 20, 2020
% - Notes:
%       1.) You can change any function in line28.
%       2.)
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all
disp('Running...')
tic

%% Define parameters
r_s = 0.0;              % the starting parameter "r"
r_e = 5.0;              % the ending parameter "r"
r_step = 0.01;          % the step of parameter "r"
x_s = -1.0;               % the lower limit of guess parameter "x"
x_e = 1.0;                % the upper limit of guess parameter "x"
iters = 1000;
nSeeds = 100;

%% Preallocating variables & functions
r = (r_s:r_step:r_e);
f = @(r, x_n) r*x_n*(1-x_n);

%% Main code
for rCount = 1:length(r)
    for nSeed = 1:nSeeds
        x_n = x_e + (x_s-x_e)*rand(1);
        
        for iter = 1:iters
            x_n = f(r(rCount), x_n);
        end
        
        data_r(rCount, nSeed) = r(rCount);
        data_x(rCount, nSeed) = x_n;
    end
end

%% Plot and save results
figure(1), set(gcf, 'WindowStyle', 'docked')
plot(data_r, data_x, 'k.')
grid on, set(gca, 'fontsize', 14)
xlabel('r', 'fontsize', 16, 'fontweight', 'bold'), ylabel('x', 'fontsize', 16, 'fontweight', 'bold')
title(['f = ', strrep(char(f),'@(r,x_n)','')])

% save images and datasets
frame = getframe(gcf);
imwrite(frame.cdata, 'Bifurcation.png')

para.r_s = r_s;     
para.r_e = r_e;         
para.r_step = r_step;   
para.x_s = x_s;    
para.x_e = x_e; 
para.iters = iters;
para.nSeeds = nSeeds;
save('para.mat', 'para')

%%
toc
disp('Done.')

