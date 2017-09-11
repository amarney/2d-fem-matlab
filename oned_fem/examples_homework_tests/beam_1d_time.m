%space
domain = [0,1];
nElement = 20; 
Mesh = mesh_generator_1d(domain, nElement);

iDegree = 2;
Fem = fem_generator_1d_lagrange(Mesh, iDegree);

iBoundaryCondition = [0,0];
[dof_u, nt] = bc_array_generator_1d(Fem, domain, iBoundaryCondition);
Fem.dof_u = dof_u; Fem.nt = nt;

% time 
tspan = [0,2]; 
tau = 1/40; % time step size, needed for fin-diff
Ntime = (tspan(2) - tspan(1))/tau; 
u0Name = @(x) zeros(size(x));

% prepare the matrices
mName = @(x) ones(size(x));
iDerivative = 2;
nGaussPoint = 4;
Mh_22 = stiffness_matrix_assembler_1d_lagrange_global(mName,...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
M_22 = Mh_22(dof_u, dof_u);

iDerivative = 0; 
nName = @(x) ones(size(x));
Nh_00 = stiffness_matrix_assembler_1d_lagrange_global(nName,...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
N_00 = Nh_00(dof_u, dof_u);


% normally do this in time loop if B/C is time dep.
essName = @(x) zeros(size(x)); % zero on l/r boundary
ess_ttName = @(x) zeros(size(x));
I = find(nt == 0);
u_e = zeros(length(Fem.point), 1);
u_e(I) = essName(Fem.point(I));
u_ett = zeros(length(Fem.point), 1); 
u_ett(I) = ess_ttName(Fem.point(I));

tmp = -Mh_22*u_e; bc_eM22 = tmp(dof_u);
tmp = -Nh_00*u_ett; bc_eN00 = tmp(dof_u);
tmp = u0Name(Fem.point);
u0 = tmp(dof_u)';
u0_double = zeros(length(u0)*2, 1);

qName = @(x,t) pi^4/100*t.^2.*sin(pi*x) + sin(pi*x)/50;








% solve semi-discrete FE by ODE solver:
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[t, u] = ode45(@beam_semi_FE_1d, tspan, u0_double, options, ...
    Mesh, Fem, qName, bc_eM22, bc_eN00, N_00, M_22);

% % form the FE solution:
% u_fe = zeros(length(t), length(Fem.point));
% I = find(Fem.nt == 0);
% for k = 1:length(t)
%     u_e = zeros(1, 2*(length(Fem.point) - 2));
%     tmp = u_e; tmp(dof_u) = u(k, :);
%     u_fe(k, :) = tmp;
% end

u_fe = u;

% post processing: evaluate at 


