function params = params_JE_hillr3bp()

%==========================================================================
%
% Builds parameter struct for the Jupiter–Europa Hill Restricted Three-Body
% Problem (HillR3BP) in the synodic frame, using nondimensional units.
%
% MODEL ID: HILLR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 5, 2026
%
% OUTPUT:              Description                                   Units
%
%  params     -   parameter struct                                  [-]
%
%==========================================================================

    %-------------------------------
    % Physical constants
    %-------------------------------
    mu_E = 3.201e3;        % Europa gravitational parameter  [km^3/s^2]
    mu_J = 1.267e8;        % Jupiter gravitational parameter [km^3/s^2]

    R_E = 1565;            % Europa radius                  [km]
    T_E = 3.552;           % Europa orbital period           [days]
    N_E = 2.05e-5;         % Europa orbital rate             [rad/s]

    % Data from HillR3BP document:
    %   mu_E = 3.201e3 km^3/s^2
    %   mu_J = 1.267e8 km^3/s^2
    %   R_E  = 1565 km
    %   T_E  = 3.552 days
    %   N_E  = 2.05e-5 rad/s
    %   Hill length scale approximately 19.600e3 km
    %   normalized Europa radius approximately 0.0796

    %-------------------------------
    % HILLR3BP scaling
    %-------------------------------
    l_hill = (mu_E/N_E^2)^(1/3);     % Hill length scale [km]

    %-------------------------------
    % Model definition
    %-------------------------------
    params.model.name = 'HILLR3BP';
    params.model.formulation = 'lagrangian';
    params.model.frame = 'synodic';
    params.model.units = 'nd';
    
    %-------------------------------
    % Physical system data
    %-------------------------------
    params.phys.primary.name = 'Jupiter';
    params.phys.secondary.name = 'Europa';

    params.phys.mu_primary = mu_J;       % [km^3/s^2]
    params.phys.mu_secondary = mu_E;     % [km^3/s^2]
    params.phys.radius_secondary = R_E;  % [km]
    params.phys.period_secondary = T_E;  % [days]
    params.phys.n_secondary = N_E;       % [rad/s]

    %-------------------------------
    % Functions
    %-------------------------------
    params.fun.eom = @eom_hillr3bp;
    params.fun.integrate = @integrate_hillr3bp;

    %-------------------------------
    % Scaling
    %-------------------------------
    params.scale.length = l_hill;        % [km]
    params.scale.time = 1/N_E;           % [s]
    params.scale.velocity = l_hill*N_E;  % [km/s]
    params.scale.n = N_E;                % [rad/s]

    % Useful nondimensional quantities
    params.scale.radius_secondary_nd = R_E/l_hill;

    %-------------------------------
    % ODE options
    %-------------------------------
    params.ode.options = odeset( ...
        'RelTol', 1e-12, ...
        'AbsTol', 1e-12);

end