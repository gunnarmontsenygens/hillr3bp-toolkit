function [t_hist, x_vec_hist, Phi_mtx_hist, i_e, t_e, x_e_vec, Phi_mtx_e] = integrate_hillr3bp(t_span, x_0_vec, params, event_fun)
%==========================================================================
%
% Integrates the dimensionless Hill Restricted Three-Body Problem (HillR3BP)
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
% The system is integrated by solving:
%
%   x_dot   = f(x)
%   Phi_dot = A(x) * Phi
%
% where:
%   - f(x) is the HillR3BP vector field
%   - A(x) = df/dx is the Jacobian matrix
%
% The STM is initialized as:
%
%   Phi(0) = I
%
% EVENT HANDLING:
% An optional event function can be provided to detect specific conditions
% during integration (e.g., plane crossings, impacts, section conditions).
% If provided, the integrator returns:
%   - event times
%   - event states
%   - corresponding STM values
%
% NOTES:
% - This is the LAGRANGIAN formulation (not Hamiltonian).
% - The system is autonomous; dynamics depend on state but not explicitly on time.
% - The STM captures the linearized flow and is used for:
%     • stability analysis (eigenvalues of monodromy matrix)
%     • invariant manifold computation
%     • sensitivity analysis and targeting
% - The normalized equations contain no mass parameter.
%
% MODEL ID: HILLR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 5, 2026
%
% INPUT:                 Description                                   Units
%
%  t_span      -   time span for integration [t0 tf]                   [-]
%  x_0_vec     -   initial state vector (6x1)                          [-]
%  params      -   parameter struct                                    [-]
%                  params.ode.options - ODE solver options             [-]
%  event_fun   -   optional event function handle                      [-]
%
% OUTPUT:                Description                                   Units
%
%  t_hist       -   time vector (Nx1)                                  [-]
%  x_vec_hist   -   state history (Nx6)                                [-]
%  Phi_mtx_hist -   STM history (Nx6x6)                                [-]
%  i_e          -   event indices                                      [-]
%  t_e          -   event times                                        [-]
%  x_e_vec      -   state at events (Ne x 6)                           [-]
%  Phi_mtx_e    -   STM at events (Ne x 6 x 6)                         [-]
%
%==========================================================================

    % Check if there's any events
    if nargin < 4
        event_fun = [];
    end
    
    % Extract base ODE options from params 
    ode_options = params.ode.options;
    
    % Implement event in case it's not empy
    if ~isempty(event_fun)
        ode_options = odeset(ode_options, 'Events', @(t,X) event_fun(t,X,params));
    end
    
    % Initialization
    x_0_vec = x_0_vec(:);
    n = length(x_0_vec);
    Phi_0_mtx = eye(n);
    X_0_vec = [x_0_vec; Phi_0_mtx(:)];
    
    % Integration
    if isempty(event_fun)

        [t_hist, X_vec_hist] = ode113( @(t,X) eom_ext_hillr3bp(t, X, params), t_span, X_0_vec, ode_options);

        % Empty event outputs
        t_e = [];
        x_e_vec = [];
        Phi_mtx_e = [];
        i_e = [];

    else
        [t_hist, X_vec_hist, t_e, X_e_vec, i_e] = ode113(@(t,X) eom_ext_hillr3bp(t, X, params), t_span, X_0_vec, ode_options);
    
        x_e_vec = X_e_vec(:,1:n);
        N_e = size(X_e_vec,1);
        Phi_mtx_e = zeros(N_e,n,n);
    
        for k = 1:N_e
            Phi_mtx_e(k,:,:) = reshape(X_e_vec(k,n+1:n+n^2), n, n);
        end
    end

    % Recover state history
    x_vec_hist = X_vec_hist(:,1:n);

    % Recover STM history as (N x n x n)
    N = size(X_vec_hist,1);
    Phi_mtx_hist = zeros(N,n,n);
    
    for k = 1:N
        Phi_mtx_hist(k,:,:) = reshape(X_vec_hist(k,n+1:n + n^2),n,n);
    end
end