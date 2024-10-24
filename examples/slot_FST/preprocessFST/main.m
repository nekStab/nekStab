%% ======================================================================= %
% This script is used to generate the Free Stream Turbulence using 200 
% modes of the continuous spectrum of the Orr-Sommerfeld-Squire 
% equation (Brandt et al. 2004).
%
% The script is divided in two macro operation:
%    1) generation of homogeneous distribution of 10 (\omega, \gamma, 
%       \beta) values on 20 shell at 20 differents wavenumbers k
%    2) Knowing the analytical value of the eigenvalue (Jacobs & 
%       Durbin 1998) the OSS system is solved to find the eigenmode.
%       The eigenmode is normalized to unit enegy.
%
%    INPUT  for (1): kkini = smallest wavenumber
%                    kkfin = largest wavenumber
%                    numk  = number of shell between kkini and kkfin
%    OUTPUT for (1): ./RESULTS/wavenumber(nfile).dat = file with \omega,
%                    \gamma, \beta
%
%    INPUT  for (2): Re = Raynolds number based on delta*
%                    lx1
%                    Ly = height of the domain.
%
%    OUTPUT for (2): ./FST_data/velocity(nfile).dat = file with y 
%                    discretization and (u,v,w) velocity profiles
%% ======================================================================= %
close all; clear all; plotres = 'false';  % for plotting pourpose
tic; fprintf('Start to find (omega, gamma, beta) basis\n')
%% == FIRTS  STEP == %

Re = 495 % Reynolds delta* at inflow
lx1 = 6 
Ney = 35 % number of elements in the vertical direction
Ny = Ney*lx1; 

delta_star = 0.28 % at inflow plane
Ly = 20           % at inflow plane
dy_min = 0.001    % at inflow plane

kkini = (2*pi)/(Ly*delta_star); %% upper bound 
kkfin = (2*pi)/(dy_min*delta_star); %%% lower bound
deltak = 0.1 % following 0.14 from Brandt 2004
numk = int8((kkfin - kkini)/deltak); % energy shells distribution

%% == SECOND  STEP == %
[nmode]=wavenumber(numk,kkini,kkfin,plotres); time = toc;
fprintf('nmode=',nmode);
fprintf('done in %f seconds \n',time); tic; fprintf('Start to find OSS modes\n')
FST_modes(Re, Ny, Ly,plotres,numk,nmode); time = toc; fprintf('done in %f minutes \n',time/60)
%% % ====== CHANGE THE DIMENSIONLESS IF YOU NEED from delta to H
% delta = 0.003679225160663; %dimensional delta
% H = 0.01;                  %dimensional h of the cylinder
% dimensionless(numk,10,delta,H)
% fprintf(' kmax = %4f and kmin = %4f \n',kkfin/delta*H,kkini/delta*H)
