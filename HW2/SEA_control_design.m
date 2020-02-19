clear

%% Define constants
% We're defining numerical values for our constants because symbolically
% computing eigenvalues, defining controllers, etc is hard for the symbolic
% toolbox to handle with matrices of this size
L1 = 1; L2 = 1; m1 = 1; m2 = 1; I1 = 1; I2 = 1;
J1 = 0.5; J2 = 0.5; K1 = 1; K2 = 1; B1 = 1; B2 = 1;
g = 9.81;

%% Set up Symbolic Variables

syms q1 q2 q3 q4 real
syms dq1 dq2 dq3 dq4 real
syms ddq1 ddq2 ddq3 ddq4 real
syms tau1 tau2 real

q = [q1; q2; q3; q4];
dq = [dq1; dq2; dq3; dq4];
ddq = [ddq1; ddq2; ddq3; ddq4];
constants = [L1, L2, m1, m2, I1, I2, J1, J2, K1, K2, B1, B2, g];
tau = [tau1; tau2];

%% Extract the dynamics from our previous solutions
% Here I'm using Roy Featherstone's notation:
%
%                            B*tau = D*ddq + H
%
% Where H = C*dq + N     (and N = G + damping)
% And B is the input matrix, which is different from the B used in damping
% The input matrix allows us to specify which degrees of freedom are
% actuated, and which are not

LHS = SEA_dynamics_dampers_gen(q, dq, ddq, constants);
D = jacobian(LHS, ddq); % This is the inertia matrix
H = LHS - D*ddq; % This is the drift term (everything else)
B = [1, 0; 0, 0; 0, 1; 0, 0]; % This is the input matrix

%% Forward Dynamics
% Here we define the forward dynamics of the system in terms of the full 
% state x = [q; dq]
% Here you should define the forward dynamics of your system dx = F(x, u)

x = [q; dq];
dx = [dq; ddq];

%%%%%%% YOUR WORK HERE %%%%%%%

ForwardDynamics = [dq;pinv(D) * (B*tau - H)];


%% Linearize Dynamics and Check Eigenvalues
% Now linearize about the equilibrium point (x, tau) = (0, 0)
% Lecture 4 and 5 slides both contain the formula
% The functions "jacobian" and "subs" will be very useful

% After this you should check the eigenvalues of your A matrix (the "eig" 
% function might be helpful). How many of them have positive real 
% components?


%%%%%%% YOUR WORK HERE %%%%%%%

Alin = jacobian(ForwardDynamics,x);
Blin = jacobian(ForwardDynamics,tau);
Alin = subs(Alin,x,[0,0,0,0,0,0,0,0]');
Blin = subs(Blin,tau,[0,0]');
eigenvalues = eig(Alin);

%% Cast variables
% We need to cast the variables to doubles rather than symbolic expressions
% in order for the future parts to work
Alin = double(Alin);
Blin = double(Blin);
eigenvalues = double(eigenvalues);

%% Controller design
% Now that you've proved that the linearization is unstable (and therefore
% that the nonlinear dynamics are unstable), you can design a controller to
% stabilize these linear dynamics. Design an LQR controller using the
% following gains:
Q = eye(numel(x));
R = eye(numel(tau));
% What are the closed loop eigenvalues of your system? Is the controlled
% system stable?

%%%%%%% YOUR WORK HERE %%%%%%%

[K,S,P] = lqr(Alin,Blin,Q,R);
closed_loop_eigenvalues = P;

%% Stability Bound
% What is the lower (slowest) bound of your controller's rate of
% convergence close to the origin?

%%%%%%% YOUR WORK HERE %%%%%%%

lower_bound = min(real(-closed_loop_eigenvalues));

%% Generate Matlab Function and Save Variables

matlabFunction(ForwardDynamics, 'File', 'ForwardDynamics_gen', 'Vars', {x, tau});

save('AG_vars.mat', 'Alin', 'Blin', 'eigenvalues', ...
    'K', 'closed_loop_eigenvalues', 'lower_bound');

