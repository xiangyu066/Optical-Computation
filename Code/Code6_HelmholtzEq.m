echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: Code6_HelmholtzEq.m
% - Author: XYZ
% - Created date: May 1, 2020
% - Modified date: 
% - Notes:
%       1.) 
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all
disp('Running...')
tic

%% Define units
global um nm
um = 1;
nm = 1E-3 *(um);

%% Determinate simulation parameters
wavelength = 1 *(um);
dx = 100 *(nm);
xRange = 20 *(um);
n_air = 1;
n_medium = 1.5;
xRange_medium = 10 *(um);
type = 1; % triangle is "1"; rectangle is "2"
nModes = 4;

%% Preallocating variables and functions
k0 = 2*pi/wavelength;
x = (-xRange/2:dx:xRange/2);
nMode = 0;

% select the distribution of the refractive index
n = zeros(size(x));
switch type
    case 1
        % the triangular distribution of the refractive index
        n(x<-xRange_medium/2) = n_air;
        n(x>xRange_medium/2) = n_air;
        n(x>=0 & x<=xRange_medium/2) = -((n_medium-n_air)/(xRange_medium/2))*x(x>=0 & x<=xRange_medium/2)+n_medium;
        n(x<=0 & x>=-xRange_medium/2) = ((n_medium-n_air)/(xRange_medium/2))*x(x<=0 & x>=-xRange_medium/2)+n_medium;
    case 2
        n(x<-xRange_medium/2) = n_air;
        n(x>xRange_medium/2) = n_air;
        n(x>=0 & x<=xRange_medium/2) = n_medium;
        n(x<=0 & x>=-xRange_medium/2) =n_medium;
end
figure(1), subplot 211, plot(x,n), ylim([0,2]), grid on
xlabel('x [\mum]','fontweight','bold'),ylabel('refractive index','fontweight','bold')

% construct the Helmholtz matrix
H0 = diag((-2/(dx^2)) + (k0*n).^2);
H1 = diag(ones(1,length(x)-1)/dx^2,1);
H2 = diag(ones(1,length(x)-1)/dx^2,-1);
H = H0+H1+H2;

%% Calculate field propagation
[E,beta_squre] = eig(H);
n_eff = sqrt(diag(beta_squre))/k0;

%%
figure(1), 
for idx = length(n_eff):-1:1
    tf = isreal(n_eff(idx));
    if (tf)
        if (n_eff(idx)>min(n)) && (n_eff(idx)<max(n)) % mode condition
            nMode = nMode+1;
            subplot 212, hold on, plot(x,E(:,idx))
            if (nMode == nModes)
                grid on
                xlabel('x [\mum]','fontweight','bold'),ylabel('E(x)','fontweight','bold')
                legend()
                break;
            end
        end
    end
end
frame = getframe(gcf);
imwrite(frame.cdata, 'Code6_HelmholtzEq_triangle.png')

%%
toc
disp('Done.')

