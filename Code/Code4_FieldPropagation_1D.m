echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: Code4_FieldPropagation_1D.m
% - Author: XYZ
% - Created date: April 10, 2020
% - Modified date: April 11, 2020
% - Notes:
%       1.) Take single slit and double slits as examples.
%       2.) At far-field, the propagation equation is equal to the Fourier transform.
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all
disp('Running...')
tic

%% Define units
global um m cm mm nm
um = 1;
m = 1E+6 *(um);
cm = 1E-2 *(m);
mm = 1E-3 *(m);
nm = 1E-9 *(m);

%% Determinate simulation parameters
wavelength = 400 *(nm);
U0 = 1;
rectSize = 100 *(um); % the width of slit
nSlits = 2;
dist_Slit = 500 *(um); % the separation distance between slits
du = 1 *(um);
Z = 1 *(m); % propagation distance
dx = 0.5 *(um);
xRange = 20 *(mm);

%% Preallocating variables and functions
k = 2*pi/wavelength;
switch nSlits
    case 1
        u = (-rectSize/2:du:rectSize/2);
    case 2
        u1 = (-rectSize/2:du:rectSize/2)-(nSlits-1)*dist_Slit/2;
        u2 = (-rectSize/2:du:rectSize/2)+(nSlits-1)*dist_Slit/2;
        u = [u1, u2];
end
x = (-xRange:dx:xRange)';

FresConst = 1/(1i*wavelength);

% using vetorization to boost calculation speed
uu = repmat(u,[length(x),1]);
xx = repmat(x,[1,length(u)]);

%% Calculate field propagation
rr = sqrt( Z^2 + (xx-uu).^2 );
Ux = FresConst * Z * sum(( U0 * exp(1i*k.*rr) ./ rr.^2) * du ,2);
UUx = abs(Ux).^2;

%%
figure, set(gcf, 'WindowStyle', 'docked')
switch nSlits
    case 1
        subplot 211, fplot(U0*rectangularPulse(min(u), max(u), sym('u')), [2*min(u), 2*max(u)]), grid on
    case 2
        subplot 211, fplot(U0*(rectangularPulse(min(u1), max(u1), sym('u'))+rectangularPulse(min(u2), max(u2), sym('u'))), [2*min(u), 2*max(u)]), grid on
end
set(gca, 'fontsize', 16)
title('U(u)', 'fontsize', 16, 'fontweight', 'bold')
xlabel('x [\mum]', 'fontsize', 16, 'fontweight', 'bold'), ylabel('Intensity [a.u.]', 'fontsize', 16, 'fontweight', 'bold')

subplot 212, plot(x/mm, UUx), grid on
set(gca, 'fontsize', 16)
title(['U(x, Z = ', num2str(Z/m), ' m)'], 'fontsize', 16, 'fontweight', 'bold')
xlabel('x [mm]', 'fontsize', 16, 'fontweight', 'bold'), ylabel('Intensity [a.u.]', 'fontsize', 16, 'fontweight', 'bold')

% save images and datasets
frame = getframe(gcf);
imwrite(frame.cdata, 'Code4_FieldPropagation_1D_double slits.png')

para.unit = 'um';
para.wavelength = wavelength;
para.U0 = 1;
para.rectSize = rectSize;
para.nSlits = nSlits;
para.dist_Slit = dist_Slit;
para.du = du;
para.Z = Z;
para.dx = dx;
para.xRange = xRange;
save('Code4_FieldPropagation_1D_double_slits_para.mat', 'para')

%%
toc
disp('Done.')

