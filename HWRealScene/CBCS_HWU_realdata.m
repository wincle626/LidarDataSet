%Andreas Assmann
%Heriot-Watt University and ST
%EngD in Applied Photonics, 2020
%% Data import: CBCS for real LiDAR data 
% Outline:
%           1. Select and load appropriate LiDAR scene (5 available)
%                   2 options: 25% pattern density and 50% pattern density at 64 patterns total                
%           2. Reformat data into A, YQ, YI 
% To adjust sample numbers for processing
% e.g. m<n with m=8
%           A(:,1:8,:), YQ(:,1:8), YI(:,1:8)
%           samples can be chosen randomly and/or from different start points
%           e.g XI = admm_lasso(A, YI, lambda, rho, alpha); for m=64 > n=16
%               -> XI = admm_lasso(A(:,4:12,:), YI(:,4:12), lambda, rho, alpha); for m=8 < n=16
%%
close all
clear all
clc
%%
ratio = '50perc'; % options: '25perc' and '50perc' - ratio of non-zeros/zeros in pattern
scene_select = 3;
%
switch scene_select
    case 1
        load(['neopropeneBCS_64_' ratio '.mat']);
    case 2
        load(['concblocksBCS_64_' ratio '.mat']);
    case 3
        load(['floatingBCS_64_' ratio '.mat']);
    case 4
        load(['sandblocksBCS_64_' ratio '.mat']);
    case 5
        load(['pipesBCS_64_' ratio '.mat']);
end
%% Reformat data
A = cbcs_data.A; YQ = cbcs_data.YQ; YI = cbcs_data.YI;
dim = cbcs_data.dim; dim_cb = cbcs_data.dim_cb;
XDref = cbcs_data.dsparseXD; % reference dSparse solution for full sequence
clim = cbcs_data.clim; % colour bar limits for scene
% Show pre-computed 
figure
imagesc(XDref);
title('Depth (dSparse) - pre-computed')
title(colorbar,'Distance, cm')
caxis manual
caxis(clim)
pbaspect([dim(1)/dim(2) 1 1])
%% END OF CODE