%% Diffusion_1D_Mars.m
%
% Simulates the 1D sub-surface heat diffusion equation with varying
% insolation and a surface solid-state greenhouse effect under martian
% conditions. Simple forward Euler is currently used for the solver. The
% code assumes a mix of basalt and H2O ice in the regolith,
% corresponding to an ice-rich region such as Deuteronilus Mensae. It takes
% into account H2O ice <-> water phase changes. The code conserves energy
% and H2O mass to machine precision. Atmospheric radiative effects (visible
% attenuation and infrared warming) are taken into account using Mars
% Climate Database data. Sensible heat exchange with the atmosphere is neglected.
%
% Robin Wordsworth 19/10/18

close all
clear all

%% set up the system

% general options
test_gaussian  = 0; % test propagation of internal gaussian pulse?
include_latent = 1; % include H2O latent heat of fusion effects?
high_lat       = 1; % choose high-latitude location?
use_MCD        = 1; % use Mars Climate Database for atmospheric forcing?
SSG_effect     = 1; % include solid-state greenhouse effect?
low_res        = 0; % use low spatial resolution?

% orbital parameters and functions
ecc   = 0.0934;         % eccentricity []
gam   = deg2rad(25.19); % obliquity []
Ls_p  = deg2rad(251);   % perihelion Ls (near summer solstice)
F0    = 1366/1.524^2;   % solar flux at Mars orbit [W/m2]

if high_lat
    phi0  = +40*pi/180; % Deuteronilus latitude [rad]
    if(use_MCD)
        load deuteronilus_MCD_data_40.mat
    end
else
    phi0  = +25*pi/180; % Arabia latitude [rad]
    if(use_MCD)
        load arabia_MCD_data.mat
    end
end
if(use_MCD)
    Ls_MCD  = [-0.5 Ls 360.5];
    GSR_MCD = [GSR(end) GSR' GSR(1)]';
    GLR_MCD = [GLR(end) GLR' GLR(1)]';
    clear GLR GSR Ls psurf Tatm taud
end

% create Insolation object
% for simplicity, we always start the run at perihelion
s2day = 86400;                % 1 Earth day [s]
Torb  = 687*s2day;            % Mars orbital period [s]
n     = 2*pi/Torb;            % mean motion [rad/s]
in    = Insolation(ecc,gam,Ls_p,Torb,F0);

% set up time/space quantities
if(low_res)
    nt    = 2e3;              % timesteps
    nz    = 100;              % layers in z (there are nz+1 levels)
else
    nt    = 8e4;              % timesteps
    nz    = 400;              % layers in z (there are nz+1 levels)
end
Lz    = 80;                   % size of domain [m]
dz    = Lz/nz;                % vertical discretization depth [m]
zm    = (0:dz:Lz)';           % level depth array [m]
z     = zm(1:nz) + dz/2;      % layer depth array [m]
if(test_gaussian)
    t_f = Torb;               % total time [s]
else
    t_f = 20*Torb;            % total time [s]
end
dt    = t_f/nt;               % time discretization interval [s]
t_a   = (1:nt)*dt;            % time array [s]

% main parameters
if(include_latent)
    Phi     = 0.5;                % regolith ice volume fraction [m3/m3]
else
    Phi     = 0.0;                % regolith ice volume fraction [m3/m3]
end
rhor        = 3e3;                % regolith rock density [kg/m3]
rhoi        = (917 + 1e3)/2;      % ice/water density (just take the approx. mean here) [kg/m3]
ch_r        = 840;                % basalt heat capacity [J/kg/K]
ch_i        = 2100;               % ice/water heat capacity [J/kg/K]
% we ignore the doubling of ch_i that occurs when ice melts to water
% this is a conservative assumption.
s.rho       = rhor*(1-Phi) + rhoi*Phi; % mean regolith density [kg/m3]
s.ch        = (ch_r*rhor*(1-Phi) + ch_i*rhoi*Phi)/s.rho; % mean regolith heat capacity [J/kg/K]
s.dz_top    = 0.025;              % aerogel layer thickness [m]
s.tau_aero  = 0.2*s.dz_top*1e2;   % aerogel visible optical depth []
s.heatcap   = s.rho*s.ch;         % kg/m3 J/kg/K = J/m3/K
s.K0        = 2.0;                % approx. thermal conductivity of frozen soil/regolith (Clifford 1993) [W/m/K]
s.Ka        = 0.01;               % approx. thermal conductivity of aerogel (Dorcheh et al. 2007) [W/m/K]
s.Fgeo      = 30e-3;              % geothermal heat flux (Clifford 1993) [W/m2]
s.sigma     = 5.67e-8;            % Stefan-Boltzmann constant [W/m2/K4]
s.L_lH2O    = 333.55*1e3;         % latent heat of fusion (H2O) [J/kg]
Twater      = 273.15;             % freezing temperature for water [K]
Ti          = zeros(nz,1);        % initial temperature profile [K]
Teq_r       = 214;                % est. regolith equilibrium temperature [K]
K           = zeros(nz,1) + s.K0; % layer thermal conductivity [W/m/K]
dTdt        = zeros(nz,1);        % rate of change of temperature [K/s]
Km          = interp1(z,K,zm);    % level thermal conductivity [W/m/K]
Km(1)       = Km(2);              % make Km const. at boundaries to avoid problems
Km(nz+1)    = Km(nz);
Am          = Km/s.heatcap;       % thermal diffusivity [m2/s]
dTdz_geo    = s.Fgeo/(s.K0);      % geothermal temperature gradient [K/m]

% set up arrays
SZA_a   = zeros(1,nt);            % solar zenith angle [degrees]
Fin_a   = zeros(1,nt);            % visible flux incident on top of SSG layer [W/m2]
Fabs_a  = zeros(1,nt);            % visible flux incident on base of SSG layer [W/m2]
GLR_a   = zeros(1,nt);            % infrared flux from atmosphere incident on top of SSG layer [W/m2]
Fout_a  = zeros(1,nt);            % net infrared flux radiated from top of SSG layer [W/m2]
Etot_a  = zeros(1,nt);            % total column energy [J/m2]
T_a     = zeros(nt,nz);           % temperature [K]
u_ice_a = zeros(nt,nz);           % local ice column density [kg/m2]
u_wat_a = zeros(nt,nz);           % local water column density [kg/m2] 


%% solve the system

if(test_gaussian)
    % define analytic function for evolution of Gaussian pulse [K]
    Tana  = @(ti,tf) 270 + 20*sqrt(ti/tf)*exp(-(z-Lz/2).^2/(4*tf*s.K0/s.heatcap));
    Ti(:) = Tana(0.1*Torb,0.1*Torb);
else
    Ti(:) = Teq_r + dTdz_geo*(z(nz)-z);
end
T = Ti;                                % initial temperature [K]

qice_i = zeros(nz,1) + Phi*rhoi/s.rho; % initial specific ice fraction [kg/kg]
qwat_i = zeros(nz,1) + 0.0;            % initial specific water fraction [kg/kg]
u_ice  = qice_i*s.rho*dz;              % initial ice column density, by layer [kg/m2]
u_wat  = qwat_i*s.rho*dz;              % initial water column density, by layer [kg/m2]

tic
for it=1:nt
    
    % calculate solar zenith angle
    M         = n*t_a(it);           % mean anomaly [rad]
    Ls        = in.kepler_eqn_fn(M); % solar longitude [rad]
    [ff, dd]  = in.insolation_fn(phi0,Ls);
    SZA       = in.meanSZA_fn(phi0,Ls);
    
    % calculate incoming solar OR interpolate MCD GSR and GLR data
    if(use_MCD)
        F_sw_dn = interp1(Ls_MCD,GSR_MCD,rad2deg(Ls));
        GLR     = interp1(Ls_MCD,GLR_MCD,rad2deg(Ls));
    else
        F_sw_dn = ff*in.F0;
        GLR     = 0;
    end
    
    % calculate infrared flux to atmosphere at top of SSG layer
    % and attenuation of visible radiation by the layer
    % we assume that the layer is in thermal equilibrium at all times,
    % and that its internal thermal balance is dominated by conduction.
    %
    % .................
    % atmosphere
    % ---------------Ta
    % SSG layer
    % ---------------Tb=T(nz)
    % subsurface
    %
    if(SSG_effect)
        f         = @(Ta) s.sigma*Ta^4 - GLR - s.Ka*(T(nz) - Ta)/s.dz_top;
        Ta_e      = fzero(f,T(nz)); % get equilibrium temperature at top of layer [K]
        % T(nz) is inputted as the starting value for root-finder function fzero
        F_ir_up   = s.sigma*Ta_e^4 - GLR;
        F_sw_abs  = F_sw_dn*exp(-s.tau_aero/cos(pi*SZA/180));
    else
        F_ir_up   = s.sigma*T(nz)^4 - GLR;
        F_sw_abs  = F_sw_dn;
    end
    
    % calculate sub-surface thermal diffusion
    % interior points
    for i=2:nz-1
        dTdt(i) = (Am(i+1)*(T(i+1) - T(i)) - Am(i)*(T(i) - T(i-1)))/dz^2;
    end
    % boundary points
    if(test_gaussian)
        % bottom BC: no flux
        dTdt(1)  = Am(2)*(T(2) - T(1))/dz^2;
        % top BC: no flux
        dTdt(nz) = -Am(nz)*(T(nz) - T(nz-1))/dz^2;
    else
        % bottom BC: geothermal heat flux
        dTdt(1)  = Am(2)*(T(2) - T(1))/dz^2 + s.Fgeo/(dz*s.heatcap);
        % top BC: SSG / insolation
        dTdt(nz) = -Am(nz)*(T(nz) - T(nz-1))/dz^2 + (F_sw_abs - F_ir_up)/(dz*s.heatcap);
    end
    
    % update temperature
    T = T + dTdt*dt;
    
    % calculate latent heat consumed for ice melting
    if(include_latent)
        
        % figure out if melting / freezing is occurring
        Tnew  = T;
        melt  = (Tnew>Twater & u_ice>0);
        frez  = (Tnew<Twater & u_wat>0);
        
        % melt the ice
        dT_try  = melt.*(Tnew-Twater);
        du_ice  = dT_try*dz*s.heatcap/s.L_lH2O;     % ice loss increment [K m kg/m3 J/kg/K = kg/m2]
        du_ice(du_ice>u_ice) = u_ice(du_ice>u_ice); % don't remove more ice than exists
        
        % freeze the water
        dT_try  = frez.*(Twater-Tnew);
        du_wat  = dT_try*dz*s.heatcap/s.L_lH2O;     % water loss increment [K m kg/m3 J/kg/K = kg/m2]
        du_wat(du_wat>u_wat) = u_wat(du_wat>u_wat); % don't remove more liquid water than exists
        
        % update ice and water amounts
        u_ice = u_ice - du_ice + du_wat;
        u_wat = u_wat - du_wat + du_ice;
        
        % now calculate actual temperature tendency per layer [K/s]
        dTdtl = s.L_lH2O*(du_wat-du_ice)/(dz*s.heatcap*dt);
        
        % final update to temperature
        T     = T + dTdtl*dt;
        
    end
    
    % calculate total energy column density [J/m2] of system
    Etot  = s.heatcap*sum(T)*dz + s.L_lH2O*sum(u_wat); 
    
    % update arrays
    SZA_a(it)     = SZA;      % solar zenith angle [degrees]
    Fin_a(it)     = F_sw_dn;  % visible flux incident on top of SSG layer [W/m2]
    Fabs_a(it)    = F_sw_abs; % visible flux incident on base of SSG layer [W/m2]
    GLR_a(it)     = GLR;      % infrared flux from atmosphere incident on top of SSG layer [W/m2]
    Fout_a(it)    = F_ir_up;  % net infrared flux radiated from top of SSG layer [W/m2]
    Etot_a(it)    = Etot;     % total column energy [J/m2]
    T_a(it,:)     = T;        % temperature [K]
    u_ice_a(it,:) = u_ice;    % local ice column density [kg/m2]
    u_wat_a(it,:) = u_wat;    % local water column density [kg/m2]
    
end
toc

%% display results
d       = z(nz) - z;                % depth (0=surface) [m]
T_f     = T;                        % final temperature [K]
F_b_net = Fabs_a + s.Fgeo - Fout_a; % net flux vs. time [W/m2]

if(test_gaussian)
    F_b_net = 0.0;
end

% temperature and water
figure(1)

subplot(2,1,1)
plot(Ti,d,'k:',T_f,d,'k',Twater+d*0,d,'k--'); axis ij; hold on
if(test_gaussian)
    plot(Tana(0.1*Torb,1.1*Torb),d,'g'); 
    legend('initial','final numerical','melting point','final analytical')
else
    legend('initial','final','melting point')
end
xlabel('temperature [K]')
ylabel('z [m]')

subplot(2,1,2)
plot(u_ice,d,'b',u_wat,d,'r'); axis ij
xlabel('[kg/m^2]')
legend('ice','water')
ylabel('z [m]')

c = colormap(parula);

% energy + H2O budget of model
figure(2)

subplot(2,2,1)
plot(t_a/Torb,Fabs_a,'b',t_a/Torb,GLR_a,'r',t_a/Torb,s.Fgeo + t_a*0,'k'); hold on
xlabel('time [Mars year]')
ylabel('flux [W/m^2]')
legend('F_{abs}','GLR','F_{geo}')

subplot(2,2,2)
plot(t_a/Torb,F_b_net,'r',t_a(1:nt-1)/Torb,diff(Etot_a)/dt,'k'); hold on
xlabel('time [Mars year]')
ylabel('flux [W/m^2]')
legend('net boundary flux','rate of internal energy change')

subplot(4,2,5)
plot(t_a/Torb,(cumsum(F_b_net)*dt-F_b_net(1)*dt - (Etot_a-Etot_a(1)))./abs(Etot_a),'k');
xlabel('time [Mars year]')
ylabel('energy error')

subplot(4,2,7)
plot(t_a/Torb,(sum(u_ice_a+u_wat_a,2)*dz - mean(sum(u_ice_a+u_wat_a,2)*dz) )/mean(sum(u_ice_a+u_wat_a,2)*dz),'b');
xlabel('time [Mars year]')
ylabel('H_2O error')

subplot(2,2,4)
plot(t_a/Torb,sum(u_ice_a+u_wat_a,2)*dz,'k'); hold on
plot(t_a/Torb,sum(u_ice_a,2)*dz,'r')
plot(t_a/Torb,sum(u_wat_a,2)*dz,'b')
xlabel('time [Mars year]')
ylabel('column mass [kg/m^2]')
legend('total','ice','water')

% main plot for the paper (Fig. 4)
figure(3)
clev = [210 230 250 Twater 280 290 300];
contourf(t_a/Torb,d,T_a',clev); axis ij; colorbar vert; hold on; colormap(c)
caxis([210 300])
contour(t_a/Torb,d,T_a',[1 1]*Twater,'w','LineWidth',2);
xlabel('time [Mars year]')
ylabel('depth [m]')
title('temperature [K]')
axis([0 t_f/Torb 0 20])
axis([0 t_f/Torb 0 30])
