function dCjdx_vec = jacobi_gradient_hillr3bp(x_vec, params)
%==========================================================================
%
% Computes the gradient of the Jacobi constant for the dimensionless
% Hill Restricted Three-Body Problem (HillR3BP) in the uniformly rotating
% synodic frame, using the Lagrangian (position–velocity) formulation.
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
% STATE DEFINITION (LAGRANGIAN FORM):
%   x_vec = [x; y; z; v_x; v_y; v_z]
%
% where (x,y,z) is position and (v_x,v_y,v_z) is velocity in the rotating
% Hill frame.
%
% JACOBI CONSTANT:
% The Jacobi constant is defined as:
%
%   C_J = 2*V(x,y,z) - (v_x^2 + v_y^2 + v_z^2)
%
% where:
%
%   V(x,y,z) = 1/r + 1/2*(3*x^2 - z^2)
%
% and:
%
%   r = sqrt(x^2 + y^2 + z^2)
%
% The gradient of the Jacobi constant is:
%
%   grad(C_J) = [2*V_x;
%                2*V_y;
%                2*V_z;
%               -2*v_x;
%               -2*v_y;
%               -2*v_z]
%
% NOTES:
% - The Jacobi gradient is commonly used in differential correction and
%   continuation methods for periodic orbit computation.
% - This implementation assumes the Lagrangian HillR3BP state convention.
%
% MODEL ID: HILLR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 7, 2026
%
% INPUT:               Description                                   Units
%
%  x_vec      -   state vector [x;y;z;v_x;v_y;v_z]                  [-]
%  params     -   parameter struct                                   [-]
%
% OUTPUT:              Description                                   Units
%
%  dCjdx_vec  -   gradient of Jacobi constant wrt x_vec             [-]
%
%==========================================================================
    
    % Extract state variables
    x   = x_vec(1);
    y   = x_vec(2);
    z   = x_vec(3);
    v_x = x_vec(4);
    v_y = x_vec(5);
    v_z = x_vec(6);

    % Distances to primary 2
    r = sqrt(x^2+y^2+z^2);


    % Derivatives of potential V
    V_x = -x/r^3+3*x;
    V_y = -y/r^3;
    V_z = -z/r^3 - z;

    % Gradient of Jacobi constant
    dCjdx_vec = 2* [ V_x;
                  V_y;
                  V_z;
                 -v_x;
                 -v_y;
                 -v_z ];

end