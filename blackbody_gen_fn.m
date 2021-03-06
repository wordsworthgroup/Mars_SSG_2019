%% blackbody_gen_fn.m
% calculate the Planck spectral irradiance
% type 1: B_lam [W/m2/um/sr]
% type 2: B_nu [W/m2/Hz/sr]
% type 3: B_nucm [W/m2/cm^-1/sr]
%
% Robin Wordsworth 9/5/18

function B = blackbody_gen_fn(T,x,type)

c     = 299792458;      % speed of light in a vacuum [m/s]
h     = 6.626070e-34;   % Planck constant [J.s]
kB    = 1.38064852e-23; % Boltzman constant [J/K]

if(type==1)
    lam = x*1e-6; % x is wavelength [um]. convert to SI units [m].
    B   = 1e-6*(2*h*c^2./lam.^5)./(exp(h*c./(lam*kB*T)) - 1); % [W/m2/um/sr]
elseif(type==2)
    nu  = x; % x is frequency [Hz].
    B   = (2*h*nu.^3/c^2)./(exp(h*nu./(kB*T)) - 1); % [W/m2/Hz/sr]
elseif(type==3)
    nu_cm = x; % x is wavenumber [cm^-1]
    B     = 1e8*(2*h*c^2*nu_cm.^3)./(exp(h*c*100*nu_cm./(kB*T)) - 1); % [W/m2/cm^-1/sr]
    % factor 1e8 to convert h*c^2 to cm units
    % remember B_nu*dnu = B_nucm*dnu_cm 
else
   fprintf('Unknown type in blackbody_gen_fn.m.') 
end

return

