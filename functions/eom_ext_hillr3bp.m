function dX_dt_vec = eom_ext_hillr3bp(t, X_vec, params)
%==========================================================================
%
% Computes the dimensionless Hill Restricted Three-Body Problem (HillR3BP)
% equations of motion in the uniformly rotating synodic frame together with
% the State Transition Matrix (STM) dynamics, using the Lagrangian
% (position–velocity) formulation.
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
% STATE DEFINITION (AUGMENTED SYSTEM):
% The augmented state combines the physical state and the STM:
%
%   X_vec = [x; y; z; v_x; v_y; v_z; vec(Phi)]
%
% where:
%   - x_vec = [x; y; z; v_x; v_y; v_z] is the 6x1 Lagrangian state
%   - Phi is the 6x6 State Transition Matrix
%   - vec(Phi) is the column-wise vectorization of Phi (36x1)
%
% Total state dimension: 6 + 36 = 42
%
% DYNAMICS:
% The physical state evolves according to:
%
%   x_dot = f(x)
%
% and the STM evolves according to:
%
%   Phi_dot = A(x) * Phi
%
% where:
%   - f(x) is the HillR3BP vector field
%   - A(x) = df/dx is the Jacobian matrix
%
% The HillR3BP equations of motion are:
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
% - This is the LAGRANGIAN formulation (not Hamiltonian).
% - The system is autonomous; A depends on state but not explicitly on time.
% - The STM captures the linearized flow and is used for:
%     • stability analysis
%     • monodromy matrix computation
%     • invariant manifolds
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
%  X_vec     -   augmented state (42x1)                               [-]
%  params    -   struct containing model/integration parameters       [-]
%
% OUTPUT:              Description                                   Units
%
%  dX_dt_vec -   time derivative of augmented state (42x1)            [-]
%
%==========================================================================

    % Initialization
    X_vec = X_vec(:);
    
    % Extract state and STM vector
    x_vec = X_vec(1:6);
    Phi_vec = X_vec(7:42);

    % EoM
    dx_dt_vec = eom_hillr3bp(t, x_vec, params);

    % STM
    Phi_mtx = reshape(Phi_vec, 6, 6);
    A_t = jacobian_hillr3bp(t, x_vec, params);
    dPhi_dt_mtx = A_t*Phi_mtx;

    % Put vectors back into Y
    dx_dt_vec = dx_dt_vec(:);
    dPhi_dt_vec = dPhi_dt_mtx(:);
    dX_dt_vec = [dx_dt_vec; dPhi_dt_vec];
end