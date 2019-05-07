%% Insolation.m
% contains methods for calculating Mars solar flux vs. orbital parameters
% Based on calculations in (and validated against): 
% Levine et al., Icarus (1977), 
% Pierrehumbert, Principles of Planetary Climate (2010),
% and the Mars Climate Database v5.3.
%
% Robin Wordsworth 25/10/18

classdef Insolation
    
    properties (SetAccess = immutable)
        ecc;  % eccentricity []
        gam;  % obliquity [rad]
        Ls_p; % Ls of perihelion [rad]
        Torb; % orbital period [s]
        F0;   % solar constant at a given semi-major axis [s]
    end
    
    methods
        
        function in = Insolation(ecc,gam,Ls_p,Torb,F0)
            % construct Insolation object
            %
            % INPUTS
            %   ecc:  eccentricity []
            %   gam:  obliquity [rad]
            %   Ls_p: Ls of perihelion [rad]
            %   Torb: orbital period [s]
            %   F0:   solar constant at a given semi-major axis [s]
            %   Ls0:  starting Ls [rad]
            
            in.ecc  = ecc;
            in.gam  = gam;
            in.Ls_p = Ls_p;
            in.Torb = Torb;
            in.F0   = F0;
            
        end
        
        function [f, del] = insolation_fn(in,phi,Ls)
            % calculate insolation fraction vs. latitude and season
            %
            % INPUTS
            % phi: latitude [rad]
            % Ls: solar longitude [rad]
            % 
            % OUTPUTS
            % f:   diurnal mean fraction of stellar flux to surface []
            % del: substellar latitude [rad]
            
            kappa = Ls - in.Ls_p; % true anomaly [rad]

            % declination / subsolar latitude (Levine et al. 1977 eqn. 5) [rad]
            % corresponds to northern summer solstice at Ls = 90
            % corresponds to southern summer solstice at Ls = 270
            del = asin(sin(Ls)*sin(in.gam));
            
            % hour angle at terminator (PPC eqn. 7.8 p. 437) [rad]
            HT  = real(acos(-tan(phi).*tan(del)));
       
            % orbital distance / semi-major axis (c.f. Levine et al. 1977 eqn. 2) []
            % minimum when Ls = Ls_p. 
            r_div_a = (1 - in.ecc^2)/(1 + in.ecc*cos(kappa)); 
            
            % insolation fraction (c.f. Levine et al. 1977 eqn. 4) [rad]
            f   = (cos(phi).*cos(del).*sin(HT) + sin(phi).*sin(del).*HT)/pi;
            f   = real(f)/r_div_a^2;
            
            return
            
        end
        
        function SZA_mean = meanSZA_fn(in,phi,Ls)
            % calculate diurnal mean solar zenith angle vs. latitude and season
            % 
            % INPUTS
            % phi: latitude [rad]
            % Ls: solar longitude [rad]
            % 
            % OUTPUTS
            % SZA_mean: diurnal mean solar zenith angle [degrees]

            nint        = 1e2;                                             % number of integration steps in day
            del         = asin(sin(Ls)*sin(in.gam));                       % declination / subsolar latitude (Levine et al. 1977 eqn. 5) [rad]
            HT          = real(acos(-tan(phi).*tan(del)));                 % hour angle at terminator (PPC eqn. 7.8) [rad]
            h           = linspace(-HT,HT,nint);                           % hour angle during daytime [rad]
            cosSZA      = sin(phi).*sin(del) + cos(phi).*cos(del).*cos(h); % solar zenith angle cosine (PPC eqn. 7.7) [rad]
            cosSZA_mean = trapz(h,cosSZA.^2)./trapz(h,cosSZA);             % daily mean solar zenith angle cosine (Hartman eqn. 2.18) [rad]
            SZA_mean    = acos(cosSZA_mean).*180/pi;                       % daily mean solar zenith angle [deg]
            
            return
        end
        
        function kappa = kepler_eqn_fn(in,M0)
            % get true anomaly kappa from mean anomaly M
            % kappa relates to the solar longitude as kappa = Ls - Ls_p
            %
            % INPUTS
            % M0: mean anomaly, or normalized time from when Ls = Ls_p [rad]
            % 
            % OUTPUTS
            % kappa: true anomaly [degrees]

            e      = in.ecc;                             % eccentricity []
            M      = @(E) E - e*sin(E) - M0;             % mean anomaly equation [rad]   
            Mguess = M0;                                 % starting guess for E in root-finding algorithm [rad]
            E      = fzero(M,Mguess);                    % eccentric anomaly [rad]
            kappa  = 2*atan(sqrt((1+e)/(1-e))*tan(E/2)); % true anomaly [rad]
            if(kappa<0 || kappa > 2*pi)
                kappa = mod(kappa,2*pi);
            end

            return
            
        end
        
    end
    
end