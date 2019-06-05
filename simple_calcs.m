%% simple_calcs.m
% Do:
% a) extreme SSG warming calculation
% b) optimal SSG warming calculation
% See Methods section of Wordsworth, Cockell & Kerber 2019, Nature Astronomy.
%
% Robin Wordsworth 16/10/18

close all
clear all

%% extreme SSG warming calculation

% The IR radiation cannot escape at all.
% Therefore we need to radiate what we received back in the near-IR/visible
% (i.e. at wavelengths of 2 um and shorter)

Fa     = 150; % absorbed visible flux [W/m2]
lam_c  = 2; % cutoff wavelength [um]
lam_VI = linspace(1e-8,lam_c,1e4); % wavelength array [um]
B_VI   = @(T) pi*blackbody_gen_fn(T,lam_VI,1); % spectral flux [W/m2/um]
bal    = @(T) trapz(lam_VI,B_VI(T)) - Fa; % balance equation [W/m2]

T0     = 100; % starting guess for surface temperature [K]
Ts     = fzero(bal,T0) % resultant surface temperature [K]

%% optimal SSG warming calculation

T     = [0.8 0.6] % transmission of 1-cm aerogel layer (tiles, particles)
tau   = -log(T) % vertical path optical depth of 1-cm aerogel layer 
h     = 0.01 % aerogel layer height [m]
alph  = tau/h % layer extinction coefficient [1/m]
muav  = 1.0; % cosine of zenith angle [rad]
hm    = muav./alph % optimal thickness for aerogel layer
kap   = 0.02; % aerogel thermal conductivity @ 1 bar pressure [W/m/K]
dT    = muav*Fa*exp(-1)./(alph*kap) % temperature difference [K]

