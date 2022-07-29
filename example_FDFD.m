clc;clear

% This example shows how to use FDFD_ver02 package to solve 2D
% visco-acoustic wave equation. First, we need to define model properties
% such as Vp, rho and Q factor models, then we need to set numerical modeling
% parameters. 
% I tried to explain this example with appropriate comments. All input and output
% variables of the functions are also explained briefly inside each function. 
%
% For more information about theory of this work please see the following paper:
%
% Amini, N. and Javaherian, A., “A MATLAB-based frequency-domain 
% finite-difference package for solving 2D visco-acoustic  wave equation”,
% Waves in Random and Complex Media, vol. 21, no. 1, pp. 161–183, 2011.
% doi:10.1080/17455030.2010.537708.
%
% Please cite the above paper when reporting, reproducing or extending the results.

%% model properties 
% Homogeneous model
nx_org = 200; % horizontal dim of model
nz_org = 200; % vertical dim of model
vp =  2000*ones(nz_org,nx_org); % velocity model
rho = 2000*ones(nz_org,nx_org); % density model
q =    100*ones(nz_org,nx_org); % Q-factor model. For visco-acoustic modeling "atten_opt" should be set to 'KF'.

% Marmousi model
% load('data/Marmousi_vp.mat')
% vp =   MARM_vp(1:2:end,1:2:end);
% rho =  310*vp.^0.25; % density model (using Gardner's equation)
% q =    vp/15; % Q-factor model (just an assumption!). For visco-acoustic modeling "atten_opt" should be set to 'KF'.

%% modeling parameters

%%% boundary conditions
% hybrid PML+ABC absorbing boundary is used for attenuation of waves reflected from
% boundaries. Absorbing boundary is applied to the left, right and bottom boundaries automatically. 
L = 30; % width of PML
alpha = 180; % amplitude of PML damping cosine function
% For top boundary we have 3 options: 
% 'PML' : hybrid PML+ABC absorbing boundary
% 'Dirichlet' : fixed boundary  
% 'Neumann' : free boundary
top_bc = 'Neumann'; % boundary condition at top of model 

%%% source parameters
sname = 'ricker'; % type of source wavelet ('ricker' or 'gaussian')
f0 = 20; % dominant frequency of the source wavelet

%%% attenuation parameters
atten_opt = 'KF'; % attenuation option ('no_atten' for non-attenuation or 'KF' for Kolsky-Futterman)
wref= 5; % reference frequency for Kolsky-Futterman attenuation mechanism

%%% other parameters
fmax = 3*f0; % maximum frequency of the simulation (usually fmax = 3*fdom for 'ricker' and 5*fdom for 'gaussian')
freq_zpad = 0; % number of zeros to pad in the frequency domain to have a smooth waveform in time domain
twrap = 5;  % damping to suppress time aliasing 

%%% discretization parameters
G = 8; % number of samples per minimum wavelength
tmax = 1; % time of simulation
lambda = vp/fmax; 
lambda_min = min(lambda(:));
dx = lambda_min/G;

%%% parallelize over frequencies
use_parfor  = 0; % if 0 sequential "for" loop over frequencies and if 1 parallelize over frequencies using "parfor" loop

%% extend grids for boundary condition
vp = ext_pml(vp,L,top_bc);
rho = ext_pml(rho,L,top_bc);
q   = ext_pml(q,  L,top_bc);
[nz,nx] = size(vp);  % extended model sizes 

%% Sources and Receivers positions

%%% Sources positions
%Sx = 10:10:nx;           % muti sources
%Sz = 35*ones(size(Sx)); % muti sources

Sx = fix(nx/2)+1; % single source
Sz = 70; % single source

%%% Receivers positions 
Rx = 1:2:nx;                    % receivers on some of grids
Rz = 25*ones(size(Rx));     % receivers on some of grids

[Rx,Rz] = meshgrid(1:nx,1:nz); % receivers on all of grids. This case is 
% good for generating wave propagation snapshots and visualizing as 
% animation. this case requires a large amount of memory for many sources

%% FDFD modeling
[pf_r,w] = fdfd(vp,rho,q,atten_opt,wref,dx,tmax,twrap,...
               sname,f0,fmax,L,alpha,top_bc,Rx(:),Rz(:),Sx(:),Sz(:),use_parfor);
                                             
%% convert frequency-domain to time-domain
[pt,t] = four2time(pf_r,tmax,twrap,freq_zpad);
dt = t(2) - t(1);

%% display options

%%% display wave propagation animation

shot_id = 1; % index number of source

if size(pt,1) == nz*nx % check if receivers are on all of grids
    figure
    for i = 1:length(t)
        tmp = reshape(pt(:,i,shot_id),nz,nx);
        imagesc(0:dx:(nx-1)*dx,0:dx:(nz-1)*dx,tmp)
        colorbar
        axis equal
        axis tight
        %caxis([-5 5]*2000) % for 'gaussian' source
        caxis([-5 5]) % for 'ricker' source
        title(['shot number = ' , num2str(shot_id) ...
               '   , time (s) = ' , num2str((i-1)*dt)])
        xlabel('x (m)'), ylabel('z (m)')   
        pause(0.1)
    end
else
    warning('you need put to receivers at all of grids for animation!')
end

%%% display shot records

shot_id = 1; % index number of source (can be 1 to length(sx))

figure
if size(pt,1) == nz*nx % check if receivers are on all of grids
    % in this case to display shot record we need to define locations
    % of the receivers
    Rx = 1:1:nx;              % horizontal index of receivers
    Rz = (5)*ones(size(Rx)); % vertical index of receivers
    
    rec_shot_ind = sub2ind([nz,nx],Rz,Rx);
    shot = squeeze(pt(rec_shot_ind,:,shot_id))';
    imagesc((Rx-1)*dx,t,-shot(:,:))
    xlabel('x (m)'), ylabel('t (s)')
    colorbar
    %caxis([-5 5]*2000) % for 'gaussian' source
    caxis([-5 5]) % for 'ricker' source
    grid on
    title(['shot number  = ' , num2str(shot_id)])
    
else % we already have defined locations of the receivers
    
    shot = squeeze(pt(:,:,shot_id))';
    imagesc((Rx-1)*dx,t,-shot(:,:))
    xlabel('x (m)'), ylabel('t (s)')
    colorbar
    %caxis([-5 5]*2000) % for 'gaussian' source
    caxis([-5 5]) % for 'ricker' source
    grid on
    title(['shot number  = ' , num2str(shot_id)])
end
