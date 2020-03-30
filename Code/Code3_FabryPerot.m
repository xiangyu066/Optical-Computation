echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: Code3_FabryPerot.m
% - Author: 莊翔淯 (Xiang-Yu Zhuang)
% - Created date: March 27, 2010
% - Modified date: April 2, 2020
% - Notes:
%       1.) T_max (相位差=0, n*pi)
%       2.) R1=R2, 反射全部變穿透, T_max = 1
%       3.) free spectral range = c0/(2*n*L*cos(theta))
%       4.) finesse = FSR/resolution, 要越大反射率要越高
%       5.) 入射角=0, FSR=c0/(2*n*L); delta(lamda_fsr)=lambda^2/(2*n*L),
%       可以檢驗是不是Fabry-perot造成的
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all
disp('Running...')
tic

%% Define units
global m nm mm deg rad
nm = 1;
m = 1E+9 *(nm);
mm = 1E-3 *(m);

deg = 1;
rad = 1;

%% Define parameters
n_i = 1.00; 
n_t = 1.46;
L = 1.0 *(mm);

wavlength_s = 500 *(nm);
wavelength_e = 501 *(nm);
dw = 0.0001 *(nm);
theta_i = 45 *(deg);

R1 = 0.8;
R2 = 0.6;

%% Preallocating variables and functions
theta_t = asind(n_i*sind(theta_i)/n_t);
lambda = (wavlength_s:dw:wavelength_e)';
phase_diff = (2*pi*n_t*L*cosd(theta_t) ./ lambda) *(rad); % phase difference

f = @(R1, R2, phi) (1-R1)*(1-R2)./((1-sqrt(R1*R2))^2+4*sqrt(R1*R2)*(sin(phi)).^2);

%% main code 
% calculate transmission
T = f(R1, R2, phase_diff);
[pks, locs, widths, proms] = findpeaks(T, lambda);

% plot wavelength versus transmission
figure(1), set(gcf, 'WindowState', 'maximized')
findpeaks(T, lambda,'Annotate','extents')
text(locs+0.005, pks+0.015, num2str(locs), 'color', [0,0.45,0.74])
text(locs+0.005, pks-0.5*proms-0.015, num2str(widths,4), 'color', [0.9,0.7,0.3])

% calculate free spectral range
fsr = zeros(length(pks)-1,1);
for nCount = 2:length(pks)
    fsr(nCount-1) = locs(nCount)-locs(nCount-1);
    
    % annotate FSR
    hold on, plot(locs(nCount-1:nCount), [1.2, 1.2], 'ko-')
    text(0.75*locs(nCount-1)+0.25*locs(nCount), 1.2-0.015, [num2str(fsr(nCount-1)), 'nm'])
end

% adjust graphic parameters
xlim([wavlength_s, wavelength_e]), ylim([0, 1.6]), yticks((0:0.2:1)), grid on
legend({'Fabry-Perot', 'peak', 'prominence', 'width (half-prominence)', 'FSR'}, 'location', 'Northeast')
set(gca, 'fontsize', 14)
xlabel('\lambda_i [nm]', 'fontsize', 16, 'fontweight', 'bold'), 
ylabel('Transmission', 'fontsize', 16, 'fontweight', 'bold')
title(['(n_i, n_t, R_1, R_2, L, \theta_i) = (', num2str(n_i, 3), ', ', num2str(n_t, 3), ', ', num2str(R1, 3), ', ', num2str(R2, 3), ', ', num2str(L/mm, 3), 'mm, ', num2str(theta_i, 2), '^\circ)'])

% save images and datasets
frame = getframe(gcf);
imwrite(frame.cdata, 'Fabry-Perot interferometer.png')

para.unit = {'nm','deg'};
para.n_i = n_i; 
para.n_t = n_t;
para.L = L;
para.wavlength_s = wavlength_s;
para.wavelength_e = wavelength_e;
para.dw = dw;
para.theta_i = theta_i;
para.R1 = R1;
para.R2 = R2;
save('para.mat', 'para')

%%
toc
disp('Done.')

