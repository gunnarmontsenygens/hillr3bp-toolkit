function dx_dt_vec =  eom_hillr3bp(t, x_vec, params)
%==========================================================================
%
% Computes the dimensionless Hill Restricted Three-Body Problem (HillR3BP)
% equations of motion in the uniformly rotating synodic frame, using the
% Lagrangian (position–velocity) formulation.
%
% MODEL DESCRIPTION:
% The HillR3BP describes the motion of a massless spacecraft in the vicinity
% of a planetary satellite that orbits a much larger primary body. The model
% is obtained as a local approximation of the restricted three-body problem
% near the smaller primary, assuming that the satellite follows a near-circular
% orbit about the planet and that the satellite-to-planet mass ratio is small.
%
% The equations are written in a synodic frame centered at the planetary
% satellite. In this frame:
%   - x is measured from the satellite center along the planet-satellite line
%   - y completes the rotating orbital-plane triad
%   - z is measured normal to the satellite orbital plane
%
% NORMALIZATION:
% The system is nondimensionalized such that:
%   - The satellite orbital angular rate about the planet is N_S = 1
%   - Time is scaled by 1/N_S
%   - Distance is scaled by the Hill length scale:
%
%       l = (mu_S/N_S^2)^(1/3)
%
% where mu_S is the gravitational parameter of the satellite.
%
% In normalized form, the HillR3BP contains no system parameter. Dynamical
% results may therefore be scaled to different planetary-satellite systems
% satisfying the Hill approximation.
%
% STATE DEFINITION (LAGRANGIAN FORM):
%   x_vec = [x; y; z; v_x; v_y; v_z]
%
% where (x,y,z) is position and (v_x,v_y,v_z) is velocity in the synodic frame.
%
% EQUATIONS OF MOTION:
%
%   x_ddot - 2*y_dot = -x/r^3 + 3*x
%   y_ddot + 2*x_dot = -y/r^3
%   z_ddot           = -z/r^3 - z
%
% where:
%
%   r = sqrt(x^2 + y^2 + z^2)
%
% Equivalently:
%
%   x_ddot =  2*y_dot + 3*x - x/r^3
%   y_ddot = -2*x_dot       - y/r^3
%   z_ddot =                 - z/r^3 - z
%
% NOTES:
% - This is the LAGRANGIAN (position–velocity) formulation of the HillR3BP.
% - The system is autonomous in the synodic frame.
% - Coriolis and tidal terms arise from the rotating Hill approximation.
% - The normalized equations contain no mass parameter.
% - The singularity at r = 0 corresponds to collision with the central
%   satellite.
%
% MODEL ID: HILLR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 5, 2026
%
% INPUT:               Description                                   Units
%
%  t         -   time (unused, included for ODE solver compatibility) [-]
%  x_vec     -   state vector (6x1)                                   [-]
%  params    -   struct containing model/integration parameters       [-]
%
% OUTPUT:              Description                                   Units
%
%  dx_dt_vec -   time derivative of state (6x1)                       [-]
%
%==========================================================================

    % Initialization
    x_vec = x_vec(:);
    dx_dt_vec = zeros(size(x_vec));

    % Extract values from x_vec
    x = x_vec(1);
    y = x_vec(2);
    z = x_vec(3);
    v_x = x_vec(4);
    v_y = x_vec(5);
    v_z = x_vec(6);
    
    % Calculate important quantities
    r_vec = [x; y; z]; r = norm(r_vec);

    % Lagrangian EoM
    dx_dt_vec(1) = v_x;
    dx_dt_vec(2) = v_y;
    dx_dt_vec(3) = v_z;
    dx_dt_vec(4) =   2*v_y - x/r^3 +3*x;
    dx_dt_vec(5) = - 2*v_x - y/r^3;
    dx_dt_vec(6) = - z - z/r^3;

end