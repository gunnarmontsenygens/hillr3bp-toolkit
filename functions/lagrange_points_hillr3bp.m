function L_pts = lagrange_points_hillr3bp()
%==========================================================================
%
% Computes the equilibrium points of the dimensionless Hill Restricted
% Three-Body Problem (HillR3BP) in the uniformly rotating synodic frame.
%
% MODEL DESCRIPTION:
% The Hill Restricted Three-Body Problem (HillR3BP) is obtained as a local
% approximation of the Circular Restricted Three-Body Problem near the
% secondary primary. The equations describe the motion of a massless
% particle in a rotating frame centered on the secondary body.
%
% The Hill approximation assumes:
%   - The secondary mass is much smaller than the primary mass
%   - Motion is localized near the secondary body
%   - The gravitational influence of the primary is approximated through
%     tidal terms
%
% EQUILIBRIUM POINTS:
% The equilibrium points are defined as the solutions to:
%
%   dV/dx = 0,  dV/dy = 0,  dV/dz = 0
%
% where V(x,y,z) is the Hill effective potential:
%
%   V(x,y,z) = 1/r + 1/2*(3*x^2 - z^2)
%
% with:
%
%   r = sqrt(x^2 + y^2 + z^2)
%
% In the planar Hill problem, the equilibrium points lie along the x-axis
% and correspond to the collinear libration points:
%
%       L1 = [ (1/3)^(1/3), 0 ]
%       L2 = [-(1/3)^(1/3), 0 ]
%
% NOTES:
% - The HillR3BP admits only two equilibrium points.
% - Both equilibrium points are saddle-center type equilibria.
% - These points are commonly used as reference locations for periodic
%   orbit and invariant manifold computations in Hill dynamics.
%
% MODEL ID: HILLR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 7, 2026
%
% OUTPUT:              Description                                   Units
%
%  L_pts     -   struct containing equilibrium points                [-]
%                L1, L2 as 6D equilibrium states
%
%==========================================================================

    L_pts.L1 = [(1/3)^(1/3); 0; 0; 0; 0; 0]';
    L_pts.L2 = [-(1/3)^(1/3); 0; 0; 0; 0; 0]';

end