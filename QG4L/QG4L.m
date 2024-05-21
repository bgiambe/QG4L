%-------------------------------------------------------------------------
% Code used for the simulations presented in the paper "Semi-analytical 
% approach to study the role of abyssal stratification in the propagation 
% of potential vorticity in a four-layer ocean basin" submitted to OS 
% (Ocean Science). The code describes vorticity advection using Quasi 
% Geostrophic (QG) formalism, with a 4-layer scheme based on in-situ 
% observation (data are available at http://doi.org/10.5281/zenodo.7871735).
% The model equations are derived with a z-coordinate scheme for the 
% derivative discreatization instead of the density coordinate as is done 
% tipically in QG layered models, to account for bigger density jumps.
% 
% Beatrice Giambenedetti, Last modified May 2024
%-------------------------------------------------------------------------

clear variables, close all;

%----------------------------------------------Settings and Scaling numbers

% location ER-0121 (Western Ionian Sea)
lat= 36.3; 
long=16.1;

T =100000; %(s)
L=8000;  %(m)
H=3000; %(m)
U0 = L/T;
W = H/T;
St = T*U0/L;
Sw = T*W/H;
omegascale=1/T;
psiscale=L*L*omegascale;
rho0= 1.0283e+03; % kg/m^(-3)
f = gsw_f(lat); % TEOS-10 function http://www.teos-10.org/
g = gsw_grav(lat); % TEOS-10 function http://www.teos-10.org/
gs = g*H/U0^2; %scaled gravity
Tf = T*f; % time scaling factor


% Grid
M=125*2;
N=125*2;

% Integration Time
NT=50000;

% Resolution 
dx=M/.125; 
dy=dx;

% Time step based on estimate of velocity scale and CFL value
dt=0.1*dx*dy/psiscale;

% Diffusion and dissipation parameters 
AH=0.1;
BH=10*dx*dy;
dissrate=0.0;

% Boundary conditions 
% bcflag=1: perdiodic
% bcflag=2: dirichlet with bcval value everywhere
% bcflag=3: periodic in x
bcflag=2;
bcval=0; 

% Perturbation mode
imode=2; 
epsilon=0;
 
% layers thicknesses
h = [200 500 1000 1300];  

% layers densities
rho = [1.0292e+03, 1.0322e+03, 1.0372e+03,  1.0416e+03]; 


% Vortex in layer 1
rotation = 'c'; % 'c' for cyclones 
% normalized layer-wise velocities
if strcmp(rotation,'c')
    U=[1 0 0 0];
else
    U=[-1 0 0 0];
end


% coefficients for coupling matrix
F = @(h,r1,r2)((H*rho0)/(h*(r2-r1)));

% coupling matrix

A = [-(rho(1)/rho0)*F(h(1),rho0,rho(1)), F(h(1),rho0,rho(1)),0,0;
    F(h(2),rho(1),rho(2)), -2*F(h(2),rho(1),rho(2)),...
    F(h(2),rho(1),rho(2)),0;
    0,F(h(3),rho(2),rho(3)),-2*F(h(3),rho(2),rho(3)),...
    F(h(3),rho(2),rho(3));
    0,0,F(h(4),rho(3),rho(4)),-F(h(4),rho(3),rho(4))];

% adimentionalization for computation
A=-A*(Tf^2/St^2/gs); 


%-------------------------------------------------------------Core function
[qphys, psiphys] = QG4L_evolution(A,U, M,N,NT,epsilon,imode,...
    bcflag,bcval,dissrate,BH,AH,dt,dx,dy,L, psiscale,omegascale);


%--------------------------------------------------------------Example Plot
% Settings 
load('cork.mat'); %  https://www.fabiocrameri.ch/colourmaps/
names = {'Layer 1', 'Layer 2','Layer 3', 'Layer 4'};
lim_val = [-10, 10];
fontsize_val = 15;
% grid definition
[ic, jc] = meshgrid(1:M, 1:N);
xx = (ic-M/2)*dx/L;
yy = (jc-N/2)*dy/L;

%iteration step in [1,10]
kk=1;

figure
tiledlayout(2,2)
for i = 1:4
    nexttile
    pcolor(xx,yy,qphys{kk}.(i)./omegascale);
    set(gca,'xlim',lim_val,'ylim',lim_val,'FontSize',fontsize_val)
    shading('interp')
    colormap(cork)
    if i==1
        caxis([-0.1,0.1])
    elseif i==2
        caxis([-0.01,0.01])
    elseif i==3
        caxis([-0.001,0.001])
    else
        caxis([-0.0005,0.0005])
    end
    c = colorbar; 
    c.Label.FontSize = fontsize_val;
    title(join(['Vorticity ', names{i}]))
end


