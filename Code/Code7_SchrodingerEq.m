echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: Code7_SchrodingerEq.m
% - Author: XYZ
% - Created date: May 8, 2020
% - Modified date: May 13, 2020
% - Notes:
%       1.) Keep units consistent.
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all
disp('Running...')
tic

%% Define units
global m nm kg sec J eV mass_e h_bar
m = 1;
kg = 1;
sec = 1;
J = 1;
nm = 1E-9 *(m);
eV = 1.60218E-19 *(J);

mass_e = 9.10938E-31 *(kg);
h_bar = (6.62607E-34)/(2*pi) *(J*sec);

%% Determinate simulation parameters
dz = 0.1 *(nm);
zRange = 24.0 *(nm);
qWell_wide = 4.0 *(nm);
V_upper = 0.18 *(eV);
V_lower = 0.00 *(eV);

%% Preallocating variables and functions
z = (-zRange/2:dz:zRange/2);

% constructure the quantum well
V = zeros(size(z));

V(z<-qWell_wide/2) = V_upper;
V(z>qWell_wide/2) = V_upper;
V(z>=0 & z<=qWell_wide/2) = V_lower;
V(z<=0 & z>=-qWell_wide/2) = V_lower;

figure(1), plot(z/nm, V/eV,'k'),grid on
xlim([min(z/nm),max(z/nm)]), ylim([min(V/eV),max(1.1*V/eV)])
xlabel('z [nm]','fontweight','bold','fontsize',16),ylabel('V (eV)','fontweight','bold','fontsize',16)

% construct the Hamiltonian operator
H0 = diag((h_bar^2/(mass_e*dz^2)) + V);
H1 = diag(ones(1,length(z)-1)*(-h_bar^2/(2*mass_e*dz^2)),1);
H2 = diag(ones(1,length(z)-1)*(-h_bar^2/(2*mass_e*dz^2)),-1);
H = H0+H1+H2;

%% Calculate field propagation
[fn,En] = eig(H);

%%
plotModes = sum(diag(En)>V_lower & diag(En)<V_upper);
figure(1),
for plotMode = 1:plotModes
    rescale = ((En(1,1)-V_lower))/(max(fn(:,plotMode))-min(fn(:,plotMode)));
    fn_rescal_plot = fn(:,plotMode)*rescale + En(plotMode,plotMode);
    hold on, plot(z/nm, fn_rescal_plot/eV, 'k--')
    hold on, plot(z/nm, En(plotMode, plotMode)*ones(size(z))/eV)
end
xlim([min(z/nm), max(z/nm)])
    
frame = getframe(gcf);
imwrite(frame.cdata, 'Code7_SchrodingerEq.png')

%%
toc
disp('Done.')
