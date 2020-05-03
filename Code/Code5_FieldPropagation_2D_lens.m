echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: Code5_FieldPropagation_2D_lens.m
% - Author: XYZ
% - Created date: April 10, 2020
% - Modified date: April 29, 2020
% - Notes:
%       1.) 
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all
disp('Running...')

%% Define units
global um m cm mm nm
um = 1;
m = 1E+6 *(um);
cm = 1E-2 *(m);
mm = 1E-3 *(m);
nm = 1E-9 *(m);

%% Determinate simulation parameters
wavelength = 400 *(nm);
PatternType = 2; % "1" is single slit; "2" is symbol-"+"
Int0 = 1;
uLeng = 20 *(um);
vLeng = 100 *(um);
p = 75 *(cm);
focalL = 50 *(cm);
du = 1 *(um);
dx_p = 250 *(um);
dx_q = 1 *(um);
xyRange_p = 100 *(mm);

%% Preallocating variables and functions
k = 2*pi/wavelength;
FresConst = 1/(1i*wavelength);
q = ((1/focalL)-(1/p))^-1;
M = q/p;

% calculate suitable range
if uLeng>=vLeng
    uvRange = 1.2*M*uLeng;
    xyRange_q = uvRange;
else
    uvRange = 1.2*M*vLeng;
    xyRange_q = uvRange;
end

% build the corordinate of the image plane
dv = du;
switch PatternType
    case 1
        [uu, vv]= meshgrid(-uvRange/2:du:uvRange/2, -uvRange/2:dv:uvRange/2);
        U0 = zeros(size(uu));
        U0(uu<=uLeng/2&uu>=-uLeng/2&vv<=vLeng/2&vv>=-vLeng/2) = Int0;
    case 2
        [uu, vv]= meshgrid(-uvRange/2:du:uvRange/2, -uvRange/2:dv:uvRange/2);
        U0 = zeros(size(uu));
        U0(uu<=uLeng/2&uu>=-uLeng/2&vv<=vLeng/2&vv>=-vLeng/2) = Int0; 
        U0(uu<=vLeng/2&uu>=-vLeng/2&vv<=uLeng/2&vv>=-uLeng/2) = Int0;
end

dy_p = dx_p;
[xx_p, yy_p] = meshgrid(-xyRange_p/2:dx_p:xyRange_p/2, -xyRange_p/2:dy_p:xyRange_p/2);

dy_q = dx_q;
[xx_q, yy_q]= meshgrid(-xyRange_q/2:dx_q:xyRange_q/2, -xyRange_q/2:dy_q:xyRange_q/2);

% 
U_p = zeros(size(yy_p,1), size(xx_p,2));
tt = exp(-1i*k*(xx_p.^2+yy_p.^2)/(2*focalL)); 
U_q = zeros(size(yy_q,1), size(xx_q,2));

%% Calculate field propagation
% propagating for a distance "p" before the lens
tic
for xi_p = 1:size(xx_p,2)
    for yi_p = 1:size(yy_p,1)
        U = 0;
        for ui = 1:size(uu,2)
            for vi = 1:size(vv,1)
                r = sqrt( p^2 + (xx_p(yi_p,xi_p)-uu(vi,ui))^2 + (yy_p(yi_p,xi_p)-vv(vi,ui))^2);
                U = U + FresConst * p * U0(vi,ui) * exp(1i*k*r) / r^2 * du * dv;
            end
        end
        U_p(yi_p,xi_p) = U;
    end
end
toc
UU_p = abs(U_p).^2;

% calculate phase through a lens
U_lens = tt.*U_p;
UU_lens = abs(U_lens).^2;

% propagating for a distance "q" after the lens
tic
for xi_q = 1:size(xx_q,2)
    for yi_q = 1:size(yy_q,1)
        U = 0;
        for xi_p = 1:size(xx_p,2)
            for yi_p = 1:size(yy_p,1)
                r = sqrt( q^2 + (xx_q(yi_q,xi_q)-xx_p(yi_p,xi_p))^2 + (yy_q(yi_q,xi_q)-yy_p(yi_p,xi_p))^2);
                U = U + FresConst * q * U_lens(yi_p,xi_p) * exp(1i*k*r) / r^2 * dx_p * dy_p;
            end
        end
        U_q(yi_q,xi_q) = U;
    end
end
toc
UU_q = abs(U_q).^2;

%% Show and save simulation results
x_p = mean(xx_p,1);
y_p = mean(yy_p,2);
U_p_x = UU_p(round(size(UU_p,1)/2),:);
U_p_y = UU_p(:,round(size(UU_p,2)/2));
U_lens_x = UU_lens(round(size(UU_lens,1)/2),:);
U_lens_y = UU_lens(:,round(size(UU_lens,1)/2));
x_q = mean(xx_q,1);
y_q = mean(yy_q,2);

figure(1), set(gcf,'windowstate','maximized')
subplot(2,2,[1 3]), imagesc(x_p/mm, y_p/mm, UU_p), daspect([1 1 1])
hold on, line(x_p,zeros(size(x_p)),'color','r')
hold on, line(zeros(size(y_p)),y_p,'color','r')
set(gca,'fontsize',16)
title('Propagating for a distance "p" before the lens','fontsize',16,'fontweight','bold')
xlabel('X_p [mm]','fontsize',24,'fontweight','bold'), ylabel('Y_p [mm]','fontsize',24,'fontweight','bold')
subplot(2,2,2), plot(x_p/mm, U_p_x), grid on
set(gca,'fontsize',16)
xlabel('X_p [mm]','fontsize',24,'fontweight','bold'), ylabel('Intensity [a.u]','fontsize',24,'fontweight','bold')
subplot(2,2,4), plot(y_p/mm, U_p_y), grid on
set(gca,'fontsize',16)
xlabel('Y_p [mm]','fontsize',24,'fontweight','bold'), ylabel('Intensity [a.u]','fontsize',24,'fontweight','bold')
frame = getframe(gcf);
imwrite(frame.cdata, 'Code5_FieldPropagation_before_Lens.png')

figure(2), set(gcf,'windowstate','maximized')
subplot(2,2,[1 3]), imagesc(x_p/mm, y_p/mm, UU_lens), daspect([1 1 1])
hold on, line(x_p,zeros(size(x_p)),'color','r')
hold on, line(zeros(size(y_p)),y_p,'color','r')
set(gca,'fontsize',16)
title('Propagating after the lens','fontsize',16,'fontweight','bold')
xlabel('X_q [mm]','fontsize',24,'fontweight','bold'), ylabel('Y_q [mm]','fontsize',24,'fontweight','bold')
subplot(2,2,2), plot(x_p/mm, U_lens_x), grid on
set(gca,'fontsize',16)
xlabel('X_{lens} [mm]','fontsize',24,'fontweight','bold'), ylabel('Intensity [a.u]','fontsize',24,'fontweight','bold')
subplot(2,2,4), plot(y_p/mm, U_lens_y), grid on
set(gca,'fontsize',16)
xlabel('Y_{lens} [mm]','fontsize',24,'fontweight','bold'), ylabel('Intensity [a.u]','fontsize',24,'fontweight','bold')
frame = getframe(gcf);
imwrite(frame.cdata, 'Code5_FieldPropagation_through_Lens.png')

figure(3), set(gcf,'windowstate','maximized')
subplot(1,2,1), imagesc(x_q/um, x_q/um, U0), daspect([1 1 1])
set(gca,'fontsize',16)
title('U_0','fontsize',16,'fontweight','bold')
xlabel('U [\mum]','fontsize',24,'fontweight','bold'), ylabel('V [\mum]','fontsize',24,'fontweight','bold')
subplot(1,2,2), imagesc(x_q/um, y_q/um, UU_q), daspect([1 1 1])
set(gca,'fontsize',16)
title('Propagating through a lens for a distance "q"','fontsize',16,'fontweight','bold')
xlabel('X_q [\mum]','fontsize',24,'fontweight','bold'), ylabel('Y_q [\mum]','fontsize',24,'fontweight','bold')
frame = getframe(gcf);
imwrite(frame.cdata, 'Code5_FieldPropagation_after_Lens.png')

%%
disp('Done.')

