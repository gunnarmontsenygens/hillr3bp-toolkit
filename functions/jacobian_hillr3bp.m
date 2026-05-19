function A = jacobian_hillr3bp(t, x_vec, params)
%==========================================================================
%
% Computes the Jacobian matrix for the dimensionless Hill Restricted
% Three-Body Problem (HillR3BP) equations of motion in the uniformly rotating
% synodic frame, using the Lagrangian (position-velocity) formulation.
%
% MODEL DESCRIPTION:
% The HillR3BP describes the motion of a massless particle in the vicinity of
% a planetary satellite that orbits a much larger primary body. The model is
% obtained as a local approximation of the restricted three-body problem near
% the smaller primary, assuming that the satellite moves on a near-circular
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
% where (x,y,z) is position and (v_x,v_y,v_z) is velocity in the synodic
% frame.
%
% JACOBIAN DEFINITION:
% The Jacobian matrix is defined as:
%
%   A = df/dx
%
% where f is the first-order HillR3BP vector field:
%
%   d/dt(x_vec) = f(x_vec)
%
% This matrix is used to propagate the State Transition Matrix (STM):
%
%   Phi_dot = A*Phi
%
% The lower-left 3x3 block of A contains the position partial derivatives of
% the HillR3BP acceleration field:
%
%   x_ddot =  2*y_dot + 3*x - x/r^3
%   y_ddot = -2*x_dot       - y/r^3
%   z_ddot =                 - z/r^3 - z
%
% where:
%
%   r = sqrt(x^2 + y^2 + z^2)
%
% NOTES:
% - This is the Lagrangian, not Hamiltonian, HillR3BP Jacobian.
% - The system is autonomous in the synodic frame, so A depends on the
%   current state but not explicitly on time.
% - Coriolis terms appear as constant velocity partial derivatives.
% - Tidal terms appear in the position partial derivatives.
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
%  x_vec     -   state vector [x;y;z;v_x;v_y;v_z] (6x1)               [-]
%  params    -   struct containing model/integration parameters       [-]
%
% OUTPUT:              Description                                   Units
%
%  A         -   Jacobian matrix df/dx (6x6)                          [-]
%
%==========================================================================

    % Initialization
    x_vec = x_vec(:);

    % Extract values from x_vec
    x = x_vec(1);
    y = x_vec(2);
    z = x_vec(3);

    % Calculate important quantities
    r_vec = [x; y; z];
    r = norm(r_vec);

    if r == 0
        error('HillR3BP singularity: r = 0.');
    end

    % Lower-left block partial derivatives
    A41 =  3 - 1/r^3 + 3*x^2/r^5;
    A42 =  3*x*y/r^5;
    A43 =  3*x*z/r^5;

    A51 =  3*x*y/r^5;
    A52 = -1/r^3 + 3*y^2/r^5;
    A53 =  3*y*z/r^5;

    A61 =  3*x*z/r^5;
    A62 =  3*y*z/r^5;
    A63 = -1 - 1/r^3 + 3*z^2/r^5;

    % Dynamics matrix A
    A = [0,   0,   0,   1,  0, 0;
         0,   0,   0,   0,  1, 0;
         0,   0,   0,   0,  0, 1;
         A41, A42, A43, 0,  2, 0;
         A51, A52, A53, -2, 0, 0;
         A61, A62, A63, 0,  0, 0];

end