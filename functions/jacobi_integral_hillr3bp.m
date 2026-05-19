function J_hist = jacobi_integral_hillr3bp(t_hist, x_vec_hist, params)
%==========================================================================
%
% Computes the Jacobi integral for the dimensionless Hill Restricted
% Three-Body Problem (HillR3BP) along a trajectory in the uniformly rotating
% synodic frame.
%
% MODEL DESCRIPTION:
% The HillR3BP describes the motion of a massless particle in the vicinity of
% a planetary satellite that orbits a much larger primary body. The model is
% obtained as a local approximation of the restricted three-body problem near
% the smaller primary, assuming that the satellite moves on a near-circular
% orbit about the planet and that the satellite-to-planet mass ratio is small.
%
% The equations are written in a synodic frame centered at the planetary
% satellite, resulting in an autonomous dynamical system.
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
% In normalized form, the HillR3BP contains no system parameter.
%
% STATE DEFINITION:
%   x_vec = [x; y; z; v_x; v_y; v_z]
%
% where (x,y,z) is position and (v_x,v_y,v_z) is velocity in the synodic frame.
%
% JACOBI INTEGRAL:
% The Jacobi integral for the HillR3BP is defined as:
%
%   J = 1/2*(v_x^2 + v_y^2 + v_z^2) - V(x,y,z)
%
% where the Hill effective potential is:
%
%   V(x,y,z) = 1/r + 1/2*(3*x^2 - z^2)
%
% and:
%
%   r = sqrt(x^2 + y^2 + z^2)
%
% The Jacobi integral is conserved along trajectories of the HillR3BP and is
% used to characterize the accessible regions of motion (zero-velocity surfaces)
% and the energy level of the trajectory.
%
% NOTES:
% - This quantity is constant along a trajectory in the absence of numerical error.
% - The formulation differs from the CR3BP Jacobi constant due to the Hill
%   approximation and absence of a mass parameter.
% - Deviations from constancy can be used to assess numerical integration accuracy.
%
% MODEL ID: HILLR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 5, 2026
%
% INPUT:               Description                                   Units
%
%  t_hist     -   time vector (Nx1)                                  [-]
%  x_vec_hist -   state history (Nx6)                                [-]
%  params     -   struct containing model parameters                 [-]
%
% OUTPUT:              Description                                   Units
%
%  J_hist     -   Jacobi integral history (Nx1)                      [-]
%
%==========================================================================

    % Initialization
    N = length(t_hist);
    J_hist = zeros(N,1);

    % Loop
    for i = 1 : N

        % Extract variables
        x_vec = x_vec_hist(i,:);
        x = x_vec(1); y = x_vec(2);  z = x_vec(3); r = sqrt(x^2+y^2+z^2);
        v_x = x_vec(4); v_y = x_vec(5); v_z = x_vec(6); 

        % Calculate J and append to list
        J = 0.5*(v_x^2+v_y^2+v_z^2) - (1/r + 0.5*(3*x^2 - z^2));
        J_hist(i) = J;

    end

end