function TDFDAWFS2DSG

% TDFDAWFS2DSG
% This is a program of Time Domain Finite Difference Acoustic Wave Field Simulating with 2-Dimension Staggered Grid.
% Written by Tche.L. from USTC, 2016,6.

clc; clear; close all;
% format long;

%% Input parameters

nx = 101;               % the number of grid nodes in x-direction.
nz = 101;               % the number of grid nodes in z-direction.
npmlz = 20;            	% the number of grid nodes in top and bottom side of PML absorbing boundary.
npmlx = 20;             % the number of grid nodes in left and right side of PML absorbing boudary.
sx = 50;                % the grid node number of source position in x-direction.
sz = 50;                % the grid node number of source position in z-direction.
dx = 5;                 % the grid node interval in x-direction; Unit: m.
dz = 5;                 % the grid node interval in z-direction; Unit: m.
nt = 500;               % the number of time nodes for wave calculating.
dt = 1e-3;            	% the time node interval; Unit: s.
nppw = 12;              % the node point number per wavelength for dominant frequency of Ricker wavelet source.
ampl = 1.0e0;           % the amplitude of source wavelet.
xrcvr = 1:3:nx;         % the grid node number in x-direction of reciver position on ground.
nodr = 3;               % half of the order number for spatial difference.

%% Determine the difference coefficients

B = [1 zeros(1,nodr - 1)]';
A = NaN*ones(nodr,nodr);
for i = 1:1:nodr
    A(i,:) = (1:2:2*nodr - 1).^(2*i - 1);
end
C = A\B;

%% Model and source

Nz = nz + 2*npmlz;
Nx = nx + 2*npmlx;

vp = 2000*ones(Nz,Nx);                                                      % the velocity of acoustic wave of model; Unit: m/s.
rho = 1000*ones(Nz,Nx);                                                     % the density of model; Unit: kg/m^3.
rho(fix(Nz/3):end,fix(Nx/2):end) = 500;
vp(fix(Nz/3):end,fix(Nx/2):end) = 1000;

f0 = min(vp(:))/(min(dx,dz)*nppw);                                          % the dominant frequency of source Ricker wavelet; Unit: Hz.
t0 = 1/f0;                                                                  % the time shift of source Ricker wavelet; Unit: s; Suggest: 0.02 if fm = 50, or 0.05 if fm = 20.
t = dt*(1:1:nt);
src = (1 - 2*(pi*f0.*(t - t0)).^2).*exp( - (pi*f0*(t - t0)).^2);           	% the time series of source wavelet.
% The source wavelet formula refers to the equations (18) of Collino and Tsogka, 2001.

%% Perfectly matched layer absorbing factor

% R = 1e-6;                                                                   % Recommend: $R = 1e-2, if $nabsr = 5; $R = 1e-3, if $nabsr = 10; $R = 1e-4, if $absr = 20.
% dpml0z = log(1/R)*3*max(vp(:))/(2*npmlz);
dpml0z = 3*max(vp(:))/dz*(8/15 - 3/100*npmlz + 1/1500*npmlz^2);
dpmlz = zeros(Nz,Nx);
dpmlz(1:npmlz,:) = (dpml0z*((npmlz: - 1:1)./npmlz).^2)'*ones(1,Nx);
dpmlz(npmlz + nz + 1:Nz,:) = dpmlz(npmlz: - 1:1,:);
dpml0x = 3*max(vp(:))/dx*(8/15 - 3/100*npmlx + 1/1500*npmlx^2);
dpmlx = zeros(Nz,Nx);
dpmlx(:,1:npmlx) = ones(Nz,1)*(dpml0x*((npmlx: - 1:1)./npmlx).^2);
dpmlx(:,npmlx + nx + 1:Nx) = dpmlx(:,npmlx: - 1:1);
% The PLM formula refers to the equations (2) and (3) of Marcinkovich and Olsen, 2003.

%% Wavefield calculating

rho1 = rho;             % or = [(rho(:,1:end - 1) + rho(:,2:end))./2 (2*rho(:,end) - rho(:,end - 1))];
rho2 = rho;             % or = [(rho(1:end - 1,:) + rho(2:end,:))./2; (2*rho(end,:) - rho(end - 1,:))];

Coeffi1 = (2 - dt.*dpmlx)./(2 + dt.*dpmlx);
Coeffi2 = (2 - dt.*dpmlz)./(2 + dt.*dpmlz);
Coeffi3 = 1./rho1./dx.*(2*dt./(2 + dt.*dpmlx));
Coeffi4 = 1./rho2./dz.*(2*dt./(2 + dt.*dpmlz));
Coeffi5 = rho.*(vp.^2)./dx.*(2*dt./(2 + dt.*dpmlx));
Coeffi6 = rho.*(vp.^2)./dz.*(2*dt./(2 + dt.*dpmlz));

% +++++++++++++++++++++++++++++++++++++ approximate coeffient ++++++++++++++++++++++++++++++++++++++
% Coeffi1 = 1 - dt.*dpmlx;
% Coeffi2 = 1 - dt.*dpmlz;
% Coeffi3 = 1./rho./dx.*dt;
% Coeffi4 = 1./rho./dz.*dt;
% Coeffi5 = rho.*(vp.^2)./dx.*dt;
% Coeffi6 = rho.*(vp.^2)./dz.*dt;
% --------------------------------------------------------------------------------------------------

NZ = Nz + 2*nodr;                                                           % All values of the outermost some columns are set to zero to be a boundary condition: all of wavefield values beyond the left and right boundary are null.
NX = Nx + 2*nodr;                                                           % All values of the outermost some rows are set to zero to be a boundary condition: all of wavefield values beyond the top and bottom boundary are null.

Znodes = nodr + 1:NZ - nodr;
Xnodes = nodr + 1:NX - nodr;
znodes = nodr + npmlz + 1:nodr + npmlz + nz;
xnodes = nodr + npmlx + 1:nodr + npmlx + nx;
nsrcz = nodr + npmlz + sz;
nsrcx = nodr + npmlx + sx;

Ut = NaN*ones(NZ,NX,2);                                                     % the wavefield value preallocation.
Uz = zeros(NZ,NX,2);                                                        % The initial condition: all of wavefield values are null before source excitation.
Ux = zeros(NZ,NX,2);                                                        % The initial condition: all of wavefield values are null before source excitation.
Vz = zeros(NZ,NX,2);                                                        % The initial condition: all of wavefield values are null before source excitation.
Vx = zeros(NZ,NX,2);                                                        % The initial condition: all of wavefield values are null before source excitation.
Psum = NaN*ones(Nz,Nx);

U = NaN*ones(nz,nx,nt);

tic;
for it = 1:1:nt
    fprintf('The calculating time node is: it = %d\n',it);
    Ux(nsrcz,nsrcx,1) = Ux(nsrcz,nsrcx,1) + ampl*src(it)./2;
    Uz(nsrcz,nsrcx,1) = Uz(nsrcz,nsrcx,1) + ampl*src(it)./2;
    Ut(:,:,1) = Ux(:,:,1) + Uz(:,:,1);
    U(:,:,it) = Ut(znodes,xnodes,1);
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(Ut(Znodes,Xnodes + i,1) - Ut(Znodes,Xnodes + 1 - i,1));
    end
    Vx(Znodes,Xnodes,2) = Coeffi1.*Vx(Znodes,Xnodes,1) - Coeffi3.*Psum;
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(Ut(Znodes + i,Xnodes,1) - Ut(Znodes + 1 - i,Xnodes,1));
    end
    Vz(Znodes,Xnodes,2) = Coeffi2.*Vz(Znodes,Xnodes,1) - Coeffi4.*Psum;
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(Vx(Znodes,Xnodes - 1 + i,2) - Vx(Znodes,Xnodes - i,2));
    end
    Ux(Znodes,Xnodes,2) = Coeffi1.*Ux(Znodes,Xnodes,1) - Coeffi5.*Psum;
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(Vz(Znodes - 1 + i,Xnodes,2) - Vz(Znodes - i,Xnodes,2));
    end
    Uz(Znodes,Xnodes,2) = Coeffi2.*Uz(Znodes,Xnodes,1) - Coeffi6.*Psum;
    Ut(:,:,1) = Ut(:,:,2);
    Uz(:,:,1) = Uz(:,:,2);
    Ux(:,:,1) = Ux(:,:,2);
    Vz(:,:,1) = Vz(:,:,2);
    Vx(:,:,1) = Vx(:,:,2);
end
toc;

%% Plotting

% Wavefield Snapshot
figure;% colormap gray;
clims = [min(U(:)) max(U(:))]./5;
for it = 1:5:nt
    imagesc((0:nx - 1).*dx,(0:nz - 1).*dz,U(:,:,it));%,clims);
    set(gca,'xaxislocation','top'); axis equal; axis([0 (nx - 1)*dx 0 (nz - 1)*dz]);
    colorbar; xlabel('x distance (m)'); ylabel('z depth (m)');
    title(sprintf('the snapshot of %.1f ms',it*dt*1e3),'position',[(nx - 1)*dx/2,(nz - 1)*dz*(1 + 0.07)]);
    pause(0.01);
end

% Synthetic Seismogram
syngram(:,:) = U(1,xrcvr,:);
synmax = max(abs(syngram(:)));
rcvrintv = (xrcvr(2) - xrcvr(1))*dx;
syngram = syngram./synmax.*(rcvrintv/2);
[nsyn,~] = size(syngram);
figure; hold on;
for i = 1:1:nsyn
    plot(syngram(i,:) + (xrcvr(i) - 1)*dx,t.*1e3);
end
xlabel('x distance (m)'); ylabel('travel time (ms)');
title('Synthetic Seismogram','position',[(xrcvr(1) + xrcvr(nsyn))*dx/2,t(end)*1e3*(1 + 0.07)]);
set(gca,'xaxislocation','top');
set(gca,'YDir','reverse');
hold off;

end

%% References

% Collino and Tsogka, 2001. Geophysics, Application of the perfectly matched absorbing layer model to the linear elastodynamic problem in anisotropic heterogeneous media.
% Marcinkovich and Olsen, 2003. Journal of Geophysical Research, On the implementation of perfectly mathced layers in a three-dimensional fourth-order velocity-stress finite difference scheme.
