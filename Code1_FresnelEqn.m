echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: Code1.m
% - Author: ²øµ¾ÒR (Xiang-Yu Zhuang)
% - Created date: March 13, 2020
% - Modified date: March 16, 2020
% - Notes:
%       1.) Using "surf()" with "shading interp" can make 3D surface pretty.
%       2.) abs(complex)^2 ~= (complex)^2
%       3.) (Fresnel equation) if (n_i < n_t), then the polarized light from 
%           the ground is parallel to the ground.
%       4.) The transmittance of white light LED (sapphire (n=1.76) stcked 
%           with GaN (n=2.2)) is about 33%. Using coating or scraped surface 
%           can enhance more transmit light.
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all
disp('Running...')
tic

%% Define units
global meter kg sec Hz deg rad
meter = 1;
kg = 1;
sec = 1;
Hz = 1/sec;
deg = 1;
rad = 1;

%% Define parameters
n_i = 2.2;                    % the incident refractive index
n_t = 1.0;                    % the transmit refractive index

theta_s = 0*(deg);            % the starting angle
theta_e = 90*(deg);           % the ending angle
theta_steps = 1000;

%% Preallocating variables
theta_i = linspace(theta_s, theta_e, theta_steps);
theta_t = asind(n_i*sind(theta_i) / n_t); % Snell's law

%% main code (Fresnel equation)
Fresnel = { @(n_i, n_t, theta_i, theta_t) (n_i*cosd(theta_i)-n_t*cosd(theta_t)) ./ (n_i*cosd(theta_i)+n_t*cosd(theta_t));
    @(n_i, n_t, theta_i, theta_t) (2*n_i*cosd(theta_i)) ./ (n_i*cosd(theta_i)+n_t*cosd(theta_t));
    @(n_i, n_t, theta_i, theta_t) (n_t*cosd(theta_i)-n_i*cosd(theta_t)) ./ (n_i*cosd(theta_t)+n_t*cosd(theta_i));
    @(n_i, n_t, theta_i, theta_t) (2*n_i*cosd(theta_i)) ./ (n_i*cosd(theta_t)+n_t*cosd(theta_i)) };

R_veri = abs(Fresnel{1}(n_i, n_t, theta_i, theta_t)).^2;
R_hori = abs(Fresnel{3}(n_i, n_t, theta_i, theta_t)).^2;
R_ratio = R_veri ./ R_hori;

%% plot and save results
figure(1), set(gcf, 'WindowStyle', 'docked')
subplot 211, plot(theta_i,R_veri,theta_i,R_hori)
grid on, set(gca, 'fontsize', 14)
xlabel('\theta_i (deg)', 'fontsize', 16, 'fontweight', 'bold'), ylabel('R', 'fontsize', 16, 'fontweight', 'bold')
legend('R_{\perp}','R_{||}','location','northwest', 'fontsize', 16, 'fontweight', 'bold')

subplot 212, semilogy(theta_i, R_ratio)
grid on, set(gca, 'fontsize', 14)
xlabel('\theta_i (deg)', 'fontsize', 16, 'fontweight', 'bold'), ylabel('R_{\perp}/R_{||}', 'fontsize', 16, 'fontweight', 'bold')

% save images and datasets
frame = getframe(gcf);
imwrite(frame.cdata, 'Fresnel eqn_case2.png')

para.n_i = n_i;
para.n_t = n_t;
para.theta_s = theta_s;
para.theta_e = theta_e;
para.theta_steps = theta_steps;
save('para2.mat', 'para')

%%
toc
disp('Done.')

