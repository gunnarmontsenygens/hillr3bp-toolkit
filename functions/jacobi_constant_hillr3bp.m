function C_hist = jacobi_constant_hillr3bp(t_hist, x_vec_hist, params)
%==========================================================================
%
% Computes the Jacobi constant for the dimensionless Hill Restricted
% Three-Body Problem (HillR3BP) in the uniformly rotating synodic frame,
% using the Lagrangian (position–velocity) formulation.
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
% STATE DEFINITION (LAGRANGIAN FORM):
%   x_vec = [x; y; z; v_x; v_y; v_z]
%
% where (x,y,z) is position and (v_x,v_y,v_z) is velocity in the synodic frame.
%
% JACOBI CONSTANT:
% The Jacobi constant is defined as:
%
%   C = -2J
%
% where the energy-like quantity J is:
%
%   J = 1/2*(v_x^2 + v_y^2 + v_z^2) - V(x,y,z)
%
% with the Hill effective potential:
%
%   V(x,y,z) = 1/r + 1/2*(3*x^2 - z^2)
%
% and:
%
%   r = sqrt(x^2 + y^2 + z^2)
%
% Equivalently:
%
%   C = 2*V(x,y,z) - (v_x^2 + v_y^2 + v_z^2)
%
% NOTES:
% - The Jacobi constant is conserved along trajectories of the HillR3BP.
% - It defines zero-velocity surfaces and allowed regions of motion.
% - The formulation differs from the CR3BP due to the Hill approximation
%   and absence of a mass parameter.
% - Deviations from constancy can be used to assess numerical integration accuracy.
%
% MODEL ID: HILLR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 5, 2026
%
% INPUT:               Description                                   Units
%
%  t_hist      -   time vector (Nx1)                                 [-]
%  x_vec_hist  -   state history (Nx6)                               [-]
%  params      -   struct containing model parameters                [-]
%
% OUTPUT:              Description                                   Units
%
%  C_hist      -   Jacobi constant history (Nx1)                     [-]
%
%==========================================================================

C_hist = -2*jacobi_integral_hillr3bp(t_hist, x_vec_hist, params);