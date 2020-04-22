echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: Code5_FieldPropagation_1D_lens.m
% - Author: XYZ
% - Created date: April 17, 2020
% - Modified date: April 22, 2020
% - Notes:
%       1.) If the variable is complex, then converting the row vector into
%           the column vector replaces "'" with "reshape()".
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
dist_Slit = 200 *(um); % the separation distance between slits
du = 1 *(um);
p = 0.75 *(m); % propagation distance
dx1 = 1 *(um);
dx2 = 1 *(um);
focalL = 50 *(cm);
x1Range = 20 *(mm);
x2Range = 400 *(um);

%% Preallocating variables and functions
FresConst = 1/(1i*wavelength);
k = 2*pi/wavelength;
q = 1/(1/focalL-1/p);

% define the information of the object
switch nSlits
    case 1
        u = (-rectSize/2:du:rectSize/2);
    case 2
        u_lhs = (-rectSize/2:du:rectSize/2)-(nSlits-1)*dist_Slit/2;
        u_rhs = (-rectSize/2:du:rectSize/2)+(nSlits-1)*dist_Slit/2;
        u = [u_lhs, u_rhs];
end

% define the observation coordinate
x1 = (-x1Range:dx1:x1Range);
x2 = (-x2Range:dx2:x2Range);

% lens phase shift
t = exp(-1i*k/(2*focalL)*x1(:).^2); 

% using vetorization to boost calculation speed
uu1 = repmat(u,[length(x1(:)),1]);
xx1 = repmat(x1(:),[1,length(u)]);

uu2 = repmat(x1,[length(x2(:)),1]);
xx2 = repmat(x2(:),[1,length(x1)]);

%% Calculate field propagation
% the field propagates in the air
rr1 = sqrt( p^2 + (xx1-uu1).^2 );
Ux1 = FresConst * p * sum(( U0 * exp(1i*k.*rr1) ./ rr1.^2) * du ,2);
Ix1 = abs(Ux1).^2;

% the field goes through the lens
Ux_lens = t.*Ux1;

% the field propagates in the air, again
rr2 = sqrt( q^2 + (xx2-uu2).^2 );
UUx_lens = repmat(reshape(Ux_lens,1,[]),[length(x2),1]);
Ux2 = FresConst * q * sum(( UUx_lens .* exp(1i*k.*rr2) ./ rr2.^2) * du ,2);
Ix2 = abs(Ux2).^2;

%%
figure(1), set(gcf, 'WindowStyle', 'docked')
switch nSlits
    case 1
        subplot 311, fplot(U0*rectangularPulse(min(u), max(u), sym('u')), [-x2Range, x2Range]), grid on
    case 2
        subplot 311, fplot(U0*(rectangularPulse(min(u_lhs), max(u_lhs), sym('u'))+rectangularPulse(min(u_rhs), max(u_rhs), sym('u'))), [-x2Range, x2Range]), grid on
end
set(gca, 'fontsize', 16)
title('U(u)', 'fontsize', 16, 'fontweight', 'bold')
xlabel('u [\mum]', 'fontsize', 16, 'fontweight', 'bold'), ylabel('Intensity [a.u.]', 'fontsize', 16, 'fontweight', 'bold')

subplot 312, plot(x1/mm, Ix1, x1/mm,abs(Ux_lens).^2,'--','linewidth',2), grid on
set(gca, 'fontsize', 16)
legend({'Before lens','After lens'})
title(['U(x_1, p = ', num2str(p/m), ' m)'], 'fontsize', 16, 'fontweight', 'bold')
xlabel('x_1 [mm]', 'fontsize', 16, 'fontweight', 'bold'), ylabel('Intensity [a.u.]', 'fontsize', 16, 'fontweight', 'bold')

subplot 313, plot(x2,Ix2), grid on
set(gca, 'fontsize', 16)
title(['U(x_2, q = ', num2str(q/m), ' m)'], 'fontsize', 16, 'fontweight', 'bold')
xlabel('x_2 [\mum]', 'fontsize', 16, 'fontweight', 'bold'), ylabel('Intensity [a.u.]', 'fontsize', 16, 'fontweight', 'bold')

% save images and datasets
frame = getframe(gcf);
imwrite(frame.cdata, 'Code5_FieldPropagation_1D_lens_double_slits.png')

para.unit = 'um';
para.wavelength = wavelength;
para.U0 = 1;
para.rectSize = rectSize;
para.nSlits = nSlits;
para.dist_Slit = dist_Slit;
para.du = du;
para.p = p;
para.dx1 = dx1;
para.dx2 = dx2;
para.focalL = focalL;
para.x1Range = x1Range;
para.x2Range = x2Range;
save('Code5_FieldPropagation_1D_lens_double_slits_para.mat', 'para')

%%
toc
disp('Done.')
