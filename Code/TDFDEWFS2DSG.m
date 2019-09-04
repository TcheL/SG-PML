function TDFDEWFS2DSG

% TDFDEWFS2DSG
% This is a program of Time Domain Finite Difference Elastic Wave Field Simulating with 2-Dimension Staggered Grid.
% Written by Tche.L. from USTC, 2016.7.

clc; clear; close all;
% format long;

%% Input parameters

nx = 159;               % the number of grid nodes in x-direction.
nz = 159;               % the number of grid nodes in z-direction.
npmlz = 20;            	% the number of grid nodes in top and bottom side of PML absorbing boundary.
npmlx = 20;             % the number of grid nodes in left and right side of PML absorbing boudary.
sx = 80;                % the grid node number of source position in x-direction.
sz = 80;                % the grid node number of source position in z-direction.
dx = 5;                 % the grid node interval in x-direction; Unit: m.
dz = 5;                 % the grid node interval in z-direction; Unit: m.
nt = 500;               % the number of time nodes for wave calculating.
dt = 1e-3;            	% the time node interval; Unit: s.
nppw = 12;              % the node point number per wavelength for dominant frequency of Ricker wavelet source.
ampl = 1.0e0;           % the amplitude of source wavelet.
xrcvr = 1:3:nx;         % the grid node number in x-direction of reciver position on ground.
nodr = 1;               % half of the order number for spatial difference.

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

vp = 2000*ones(Nz,Nx);                                                      % the velocity of P-wave of model; Unit: m/s.
vs = 1000*ones(Nz,Nx);                                                      % the velocity of S-wave of model; Unit: m/s.
rho = 1000*ones(Nz,Nx);                                                     % the density of model; Unit: kg/m^3.
% vp(fix(Nz/3):end,fix(Nx/2):end) = 1500;

lmd = rho.*(vp.^2 - 2.*vs.^2);                                              % the lame parameter lambda of elastic wave of model.
mu = rho.*vs.^2;                                                            % the lame parameter mu of elastic wave of model.

f0 = min(vs(:))/(max(dx,dz)*nppw);                                          % the dominant frequency of source Ricker wavelet; Unit: Hz.
t0 = 1/f0;                                                                  % the time shift of source Ricker wavelet; Unit: s; Suggest: 0.02 if fm = 50, or 0.05 if fm = 20.
t = dt*(1:1:nt);
src = (1 - 2*(pi*f0.*(t - t0)).^2).*exp( - (pi*f0*(t - t0)).^2);           	% the time series of source wavelet.
% The source wavelet formula refers to the equations (18) of Collino and Tsogka, 2001.

%% Perfectly matched layer absorbing factor

% R = 1e-6;                                                                   % Recommend: $R = 1e-2, if $npmlr = 5; $R = 1e-3, if $npmlr = 10; $R = 1e-4, if $npmlr = 20.
% dpml0z = log(1/R)*3*max(vs(:))/(2*npmlz);
dpml0z = 3*max(vs(:))/dz*(8/15 - 3/100*npmlz + 1/1500*npmlz^2);
dpmlz = zeros(Nz,Nx);
dpmlz(1:npmlz,:) = (dpml0z*((npmlz: - 1:1)./npmlz).^2)'*ones(1,Nx);
dpmlz(npmlz + nz + 1:Nz,:) = dpmlz(npmlz: - 1:1,:);
dpml0x = 3*max(vs(:))/dx*(8/15 - 3/100*npmlx + 1/1500*npmlx^2);
dpmlx = zeros(Nz,Nx);
dpmlx(:,1:npmlx) = ones(Nz,1)*(dpml0x*((npmlx: - 1:1)./npmlx).^2);
dpmlx(:,npmlx + nx + 1:Nx) = dpmlx(:,npmlx: - 1:1);
% The PML formula refers to the equations (2) and (3) of Marcinkovich and Olsen, 2003.

%% Wavefield calculating

Coeffi1 = (2 - dt.*dpmlx)./(2 + dt.*dpmlx);
Coeffi2 = (2 - dt.*dpmlz)./(2 + dt.*dpmlz);
Coeffi3 = 2*dt./(2 + dt.*dpmlx)./rho./dx;
Coeffi4 = 2*dt./(2 + dt.*dpmlz)./rho./dz;
Coeffi5 = 2*dt./(2 + dt.*dpmlx).*(lmd + 2.*mu)./dx;
Coeffi6 = 2*dt./(2 + dt.*dpmlz).*lmd./dz;
Coeffi7 = 2*dt./(2 + dt.*dpmlx).*lmd./dx;
Coeffi8 = 2*dt./(2 + dt.*dpmlz).*(lmd + 2.*mu)./dz;
Coeffi9 = 2*dt./(2 + dt.*dpmlx).*mu./dx;
Coeffi0 = 2*dt./(2 + dt.*dpmlz).*mu./dz;

% +++++++++++++++++++++++++++++++++++++ approximate coeffient ++++++++++++++++++++++++++++++++++++++
% Coeffi1 = 1 - dt.*dpmlx;
% Coeffi2 = 1 - dt.*dpmlz;
% Coeffi3 = dt./rho./dx;
% Coeffi4 = dt./rho./dz;
% Coeffi5 = (lmd + 2.*mu).*dt./dx;
% Coeffi6 = lmd.*dt./dz;
% Coeffi7 = lmd.*dt./dx;
% Coeffi8 = (lmd + 2.*mu).*dt./dz;
% Coeffi9 = mu.*dt./dx;
% Coeffi0 = mu.*dt./dz;
% --------------------------------------------------------------------------------------------------

NZ = Nz + 2*nodr;
NX = Nx + 2*nodr;

Znodes = nodr + 1:NZ - nodr;
Xnodes = nodr + 1:NX - nodr;
znodes = nodr + npmlz + 1:nodr + npmlz + nz;
xnodes = nodr + npmlx + 1:nodr + npmlx + nx;
nsrcz = nodr + npmlz + sz;
nsrcx = nodr + npmlx + sx;

vxt = zeros(NZ,NX);
vxx = zeros(NZ,NX);
vxz = zeros(NZ,NX);
vzt = zeros(NZ,NX);
vzx = zeros(NZ,NX);
vzz = zeros(NZ,NX);
txxt = zeros(NZ,NX);
txxx = zeros(NZ,NX);
txxz = zeros(NZ,NX);
tzzt = zeros(NZ,NX);
tzzx = zeros(NZ,NX);
tzzz = zeros(NZ,NX);
txzt = zeros(NZ,NX);
txzx = zeros(NZ,NX);
txzz = zeros(NZ,NX);
Psum = NaN*ones(Nz,Nx);

P = NaN*ones(nz,nx,nt);

tic;
for it = 1:1:nt
    fprintf('The calculating time node is: it = %d\n',it);
    %% load source
    txxx(nsrcz,nsrcx) = txxx(nsrcz,nsrcx) + ampl*src(it)./4;
    txxz(nsrcz,nsrcx) = txxz(nsrcz,nsrcx) + ampl*src(it)./4;
    tzzx(nsrcz,nsrcx) = tzzx(nsrcz,nsrcx) + ampl*src(it)./4;
    tzzz(nsrcz,nsrcx) = tzzz(nsrcz,nsrcx) + ampl*src(it)./4;
    txxt(:,:) = txxx(:,:) + txxz(:,:);
    tzzt(:,:) = tzzx(:,:) + tzzz(:,:);
    P(:,:,it) = txxt(znodes,xnodes);
%     P(:,:,it) = tzzt(znodes,xnodes);
%     P(:,:,it) = txzt(znodes,xnodes);
    %% calculate $v_x$
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(txxt(Znodes,Xnodes + i - 1) - txxt(Znodes,Xnodes - i));
    end
    vxx(Znodes,Xnodes) = Coeffi1.*vxx(Znodes,Xnodes) + Coeffi3.*Psum;
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(txzt(Znodes + i - 1,Xnodes) - txzt(Znodes - i,Xnodes));
    end
    vxz(Znodes,Xnodes) = Coeffi2.*vxz(Znodes,Xnodes) + Coeffi4.*Psum;
    vxt(:,:) = vxx(:,:) + vxz(:,:);
%     P(:,:,it) = vxt(znodes,xnodes);
    %% calculate $v_z$
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(txzt(Znodes,Xnodes + i) - txzt(Znodes,Xnodes - i + 1));
    end
    vzx(Znodes,Xnodes) = Coeffi1.*vzx(Znodes,Xnodes) + Coeffi3.*Psum;
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(tzzt(Znodes + i,Xnodes) - tzzt(Znodes - i + 1,Xnodes));
    end
    vzz(Znodes,Xnodes) = Coeffi2.*vzz(Znodes,Xnodes) + Coeffi4.*Psum;
    vzt(:,:) = vzx(:,:) + vzz(:,:);
%     P(:,:,it) = vzt(znodes,xnodes);
    %% calculate $\tau_{xx}$ and $\tau_{zz}$
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(vxt(Znodes,Xnodes + i) - vxt(Znodes,Xnodes - i + 1));
    end
    txxx(Znodes,Xnodes) = Coeffi1.*txxx(Znodes,Xnodes) + Coeffi5.*Psum;
    tzzx(Znodes,Xnodes) = Coeffi1.*tzzx(Znodes,Xnodes) + Coeffi7.*Psum;
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(vzt(Znodes + i - 1,Xnodes) - vzt(Znodes - i,Xnodes));
    end
    txxz(Znodes,Xnodes) = Coeffi2.*txxz(Znodes,Xnodes) + Coeffi6.*Psum;
    tzzz(Znodes,Xnodes) = Coeffi2.*tzzz(Znodes,Xnodes) + Coeffi8.*Psum;
    txxt(:,:) = txxx(:,:) + txxz(:,:);
    tzzt(:,:) = tzzx(:,:) + tzzz(:,:);
    %% calculate $\tau_{xz}$
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(vzt(Znodes,Xnodes + i - 1) - vzt(Znodes,Xnodes - i));
    end
    txzx(Znodes,Xnodes) = Coeffi1.*txzx(Znodes,Xnodes) + Coeffi9.*Psum;
    Psum(:,:) = 0;
    for i = 1:1:nodr
        Psum = Psum + C(i).*(vxt(Znodes + i,Xnodes) - vxt(Znodes - i + 1,Xnodes));
    end
    txzz(Znodes,Xnodes) = Coeffi2.*txzz(Znodes,Xnodes) + Coeffi0.*Psum;
    txzt(:,:) = txzx(:,:) + txzz(:,:);
end
toc;

%% Plotting

% Wavefield Snapshot
figure;% colormap gray;
clims = [min(P(:)) max(P(:))]./5;
for it = 1:5:nt
    imagesc((0:nx - 1).*dx,(0:nz - 1).*dz,P(:,:,it));%,clims);
    set(gca,'xaxislocation','top'); axis equal; axis([0 (nx - 1)*dx 0 (nz - 1)*dz]);
    colorbar; xlabel('x distance (m)'); ylabel('z depth (m)');
    title(sprintf('the snapshot of %.1f ms',it*dt*1e3),'position',[(nx - 1)*dx/2,(nz - 1)*dz*(1 + 0.07)]);
    pause(0.01);
end

% Synthetic Seismogram
syngram(:,:) = P(1,xrcvr,:);
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
