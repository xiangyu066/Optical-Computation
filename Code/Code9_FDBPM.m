function Code9_FDBPM
echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: Code9_FDBPM.m
% - Author: XYZ
% - Created date: June 5, 2020
% - Modified date: June 16, 2020
% - Notes:
%       1.) A*x=b, x=inv(A)*b=A\b;
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
wavelength = 1.55 *(um);
dx = 0.5 *(um);
dz = 1 *(um);
X = 90 *(um);
Z = 510 *(um);
n_air = 1;
n_medium = 2;

%% Preallocating variables and functions
k0 = 2*pi/wavelength;
aw = 1/(dx)^2;
ax = -2/(dx)^2;
ae = 1/(dx)^2;

% Construct the map of the refractive index
[nxx,nzz] = meshgrid(0:dx:X,0:dz:Z); % the corordinate of the refractive index
nn = n_air*ones(size(nxx));
nn(nxx>=20*um & nxx<=25*um) = n_medium;
nn(nxx>25*um & nxx<=70*um & nzz>=10*um) = n_medium;

figure, mesh(nxx,nzz,nn), colormap('jet'), colorbar
xlim([0,X]), ylim([0,Z])
xlabel('X [\mum]'), ylabel('Z [\mum]')
title('The distribution of refractive index')
view([0,90])

% Given the initial field
phi_x = (0:dx:X); % the corordinate of the initial field
phi0 = zeros(size(phi_x))+1e-4;
phi0(phi_x>=20*um & phi_x<=25*um) = 1;

figure, plot(phi_x,phi0)
xlabel('X [\mum]'), ylabel('Amplitude')
title('The initial field')

phi = zeros(size(nxx));
phi(1,:) = phi0;
ixs = length(phi_x);

A = zeros(1,ixs);
B = zeros(1,ixs);
C = zeros(1,ixs-1);
D = zeros(1,ixs);

%% FD-BPM
for iz = 1:size(nzz,1)-1
    n_effs = HelmoholtzEq(nxx(iz+1,:),dx,k0,nn(iz+1,:),n_air,n_medium);
    n_eff = n_effs(1); % choose fundamental mode
    beta = k0*n_eff;
    for ix = 1:ixs
        if (ix == 1)
            gammaL = phi(iz,1)/phi(iz,2);
            B(ix) = -aw*gammaL+( -ax+4*(1i)*beta/dz-k0^2*(nn(iz+1,ix)^2-n_eff^2) );
            C(ix) = -ae;
            D(ix) = aw*gammaL*phi(iz,ix)+( ax+4*(1i)*beta/dz + k0^2*(nn(iz,ix)^2-n_eff^2) )*phi(iz,ix)+ae*phi(iz,ix+1);
        elseif (ix == ixs)
            gammaR = phi(iz,ix)/phi(iz,ix-1);
            A(ix) = -aw;
            B(ix) = -ae*gammaR + (-ax+4*(1i)*beta/dz-k0^2*(nn(iz+1,ix)^2-n_eff^2));
            D(ix) = aw*phi(iz,ix-1)+(ae*gammaR+ax+4*(1i)*beta/dz+k0^2*(nn(iz,ix)^2-n_eff^2))*phi(iz,ix);
        else
            A(ix) = -aw;
            B(ix) = -ax + 4*(1i)*beta/dz - k0^2*(nn(iz+1,ix)^2-n_eff^2);
            C(ix) = -ae;
            D(ix) = aw*phi(iz,ix-1) + (ax+4*(1i)*beta/dz+k0^2*(nn(iz,ix)^2-n_eff^2))*phi(iz,ix) + ae*phi(iz,ix+1);
        end
    end
    
    % solve next field
    HB = diag(B);
    HC = diag(C(1:ixs-1),1);
    HA = diag(A(2:ixs),-1);
    H = HA+HB+HC;
    
    new_phi = H\D.';
    phi(iz+1,:) = new_phi.';
end
Intensity = abs(phi).^2;

figure, mesh(nxx,nzz,Intensity), colormap('jet'), colorbar
xlim([0,X]), ylim([0,Z])
xlabel('X [\mum]'), ylabel('Z [\mum]')
title('The evolution of field in the medium')
view([0,90])

%%
toc, disp('Done...')
end

function n_effs = HelmoholtzEq(x,dx,k0,n,n_air,n_medium)
H0 = diag((-2/(dx^2)) + (k0*n).^2);
H1 = diag(ones(1,length(x)-1)/dx^2,1);
H2 = diag(ones(1,length(x)-1)/dx^2,-1);
H = H0+H1+H2;

[E,beta_squre] = eig(H);
n_effs = sqrt(diag(beta_squre))/k0;
n_effs(n_effs<=n_air | n_effs>=n_medium) = [];
end
