addpath('../oned_fem')

%% Stokes problem:

%% Make the fluid mesh:
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 5; ny = 5;
MeshFluid = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,ny);
[edge, M_e, M_ne] = mesh_edge_generator(MeshFluid,1,1);
MeshFluid.edge = edge; MeshFluid.M_e = M_e; MeshFluid.M_ne = M_ne;

%% Make the Fem data structures using Taylor-hood P2-P1:
uDegree = 2; pDegree = 1;
uFem = fem_generator_2d_lagrange(MeshFluid,uDegree);
pFem = fem_generator_2d_lagrange(MeshFluid,pDegree);

%% Locate the boundary nodes to make dof_u, dof_p, nt_u, and nt_p arrays:
fluidBoundaryName = @(x,y) (x-xmin).*(x-xmax).*(y-ymin).*(y-ymax);
fluidVelocityBoundaryName = @(x,y) (x-xmin).*(x-xmax).*(y-ymin).*(y-ymax);
fluidPressureBoundaryName = @(x,y) ones(size(x)); % ones since no BC
[dof_u, nt_u] = bc_array_generator_2d(uFem, fluidBoundaryName, fluidVelocityBoundaryName);
[dof_p, nt_p] = bc_array_generator_2d(pFem, fluidBoundaryName, fluidPressureBoundaryName);

%% Set up the boundary conditions:
essbc1Name = @(x,y,t) (1 - cos(t).*cos(t)).*((x >= xmin) & (x <= xmax) & (y == ymax) & (t <= pi/2)) ...
    + 1.*((x >= xmin) & (x <= xmax) & (y == ymax) & (t > pi/2)); % lid driven (time ramp)
essbc2Name = @(x,y) zeros(size(x));                      % y component always zero
u1essential = zeros(length(uFem.point),1);
u2essential = zeros(length(uFem.point),1);
I = find(nt_u == 0);
u1essential(I) = essbc1Name(uFem.point(1,I)', uFem.point(2,I)', 0); % t = 0
u2essential(I) = essbc2Name(uFem.point(1,I), uFem.point(2,I));

%% Prepare the matrix blocks:
aName = @(x,y) ones(size(x));  nQuadraturePoint = 7;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0;
Mhuu_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, MeshFluid, ...
    uFem, xDerivative1,yDerivative1, uFem,xDerivative2,yDerivative2,nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Mhuu_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, MeshFluid, ...
    uFem, xDerivative1,yDerivative1, uFem,xDerivative2,yDerivative2,nQuadraturePoint);

xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Mhup_hx0 = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, MeshFluid, ...
    uFem, xDerivative1,yDerivative1, pFem ,xDerivative2,yDerivative2,nQuadraturePoint);
Mhpu_h0x = Mhup_hx0';
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 0;
Mhup_hy0 = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, MeshFluid, ...
    uFem, xDerivative1,yDerivative1, pFem ,xDerivative2,yDerivative2,nQuadraturePoint);
Mhpu_h0y = Mhup_hy0';

Muu_hxx = Mhuu_hxx(dof_u, dof_u);
Muu_hyy = Mhuu_hyy(dof_u, dof_u);
Mup_hx0 = Mhup_hx0(dof_u, dof_p); Mpu_h0x = Mup_hx0';
Mup_hy0 = Mhup_hy0(dof_u, dof_p); Mpu_h0y = Mup_hy0';

%% Assemble the blocks to form the matrix for the FE system:
nu = 1; rho = 1; % fluid parameters
SP0_u = sparse(length(dof_u), length(dof_u));
SP0_p = sparse(length(dof_p), length(dof_p));
M_tmp = [nu*(Muu_hxx + Muu_hyy), SP0_u, -(1/rho)*Mup_hx0; ...
    SP0_u, nu*(Muu_hxx + Muu_hyy), -(1/rho)*Mup_hy0; ...
    -(1/rho)*Mpu_h0x, -(1/rho)*Mpu_h0y, SP0_p ];

%% Append the gamma 1 vector (lagrange multiplier) to make system consistent
% (one alternative is to instead introduce a penalty parameter)
oneName = @(x,y) ones(size(x));
xDerivative = 0; yDerivative = 0;
oneh_h = load_vector_assembler_2d_lagrange_tri_global(oneName, MeshFluid, pFem, xDerivative, ...
    yDerivative, nQuadraturePoint);
one_vec = oneh_h(dof_p);
tmp_vec = [zeros(2*length(dof_u),1); one_vec];
M = [M_tmp, tmp_vec; tmp_vec', 0];

%% Prepare the vector blocks:
f1Name = @(x,y) zeros(size(x));  % no body forces
f2Name = @(x,y) zeros(size(x));  % no body forces
xDerivative = 0; yDerivative = 0;
f1h_h = load_vector_assembler_2d_lagrange_tri_global(f1Name, MeshFluid, uFem, xDerivative, ...
    yDerivative, nQuadraturePoint);
f1_h = f1h_h(dof_u);
f2h_h = load_vector_assembler_2d_lagrange_tri_global(f2Name, MeshFluid, uFem, xDerivative, ...
    yDerivative, nQuadraturePoint);
f2_h = f2h_h(dof_u);

tmp = -nu*Mhuu_hxx*u1essential; bc_u1_eMhxx = tmp(dof_u);
tmp = -nu*Mhuu_hyy*u1essential; bc_u1_eMhyy = tmp(dof_u);
tmp = -nu*Mhuu_hxx*u2essential; bc_u2_eMhxx = tmp(dof_u);
tmp = -nu*Mhuu_hyy*u2essential; bc_u2_eMhyy = tmp(dof_u);
tmp = -(1/rho)*Mhpu_h0x*u1essential; bc_u1_eMh0x = tmp(dof_p);
tmp = -(1/rho)*Mhpu_h0y*u2essential; bc_u2_eMh0y = tmp(dof_p);

%% Assemble the right hand side vector using the vector blocks:
rhs = [f1_h + bc_u1_eMhxx + bc_u1_eMhyy; ...
    f2_h + bc_u2_eMhxx + bc_u2_eMhyy; ...
    -bc_u1_eMh0x - bc_u2_eMh0y; ...
    0];

%% Solve the FE system:
u1u2p = M\rhs;
u1 = u1u2p(1:length(dof_u));
u1_fe = u1essential; u1_fe(dof_u) = u1; % u1_fe is FE solution in x direction
u2 = u1u2p(length(dof_u) + 1:(2*length(dof_u)));
u2_fe = u2essential; u2_fe(dof_u) = u2; % u2_fe is FE solution in x direction
p = u1u2p((2*length(dof_u))+1:(end-1));
p_fe = p; % p_fe is FE solution for pressure

%% START BEAM PROBLEM:
%% Make the beam Mesh:
domainBeam = [0,1]; nEleBeam = nx; %meshes must coinc
MeshBeam = mesh_generator_1d(domainBeam,nEleBeam);

%% Make the Fem data structure using C1 elements:
iDegree = 3; hDegree = 1;
wFem = femherm1d(MeshBeam, iDegree);

%% Set up the boundary conditions by making dof_u and nt arrays:
iBoundaryCondition = [0, 0, 0, 0]; % clamp at both ends of beam
[dof_w, nt_w] = bc_array_generator_1d(wFem, domainBeam, iBoundaryCondition);
wFem.dof_w = dof_w; wFem.nt = nt_w;

%% Make the time-independant A_22 matrix:
aName = @(x) ones(size(x));
iDerivative = 2; nGaussPoint = 5;
Ah_22 = biformglobherm1d(aName,...
    MeshBeam, wFem, iDerivative, wFem, iDerivative, nGaussPoint);
A_22 = Ah_22(dof_w, dof_w);

%% Make the time-independant A_00 matrix:
aName = @(x) ones(size(x));
iDerivative = 0; nGaussPoint = 5;
Ah_00 = biformglobherm1d(aName,...
    MeshBeam, wFem, iDerivative, wFem, iDerivative, nGaussPoint);
A_00 = Ah_00(dof_w, dof_w);

%% Set beam boundary conditions:
essName = @(x) [0, 0, 0, 0]; % first two zeros are beam positions, second two are spatial derivatives
I = find(nt_w == 0); 
w_e = zeros(length(wFem.point), 1); % essDttName possibly unneeded?
essDttName = @(x) [0, 0, 0, 0]; % first two zeros are beam positions, second two are spatial derivatives
w_edtt = zeros(length(wFem.point), 1);

%% Make the time-independant b_A22 vector (fixed beam boundaries)
w_e(I) = essName(wFem.point(I));
tmp = -Ah_22*w_e;
b_A22 = tmp(dof_w);

%% Make the time-independant b_A00 vector (fixed beam boundaries)
w_e(I) = essDttName(wFem.point(I));
tmp = -Ah_00*w_edtt;
b_A00 = tmp(dof_w);

%% Set the time domain and time-step size (tau):
tspan = [0, 5]; tau = 1/40;
Ntime = (tspan(2) - tspan(1))/tau; % number of time steps

%% Set up arrays that contain time and solution values:
t = zeros(Ntime + 1, 1); t(1) = tspan(1);
w = zeros(Ntime + 1, length(dof_w)); %y-displacement!!!  
wt = zeros(Ntime + 1, length(dof_w)); % y-velocity!!!
w_fe = zeros(length(t), length(wFem.point)); %y-displacement
wt_fe = zeros(length(t), length(wFem.point)); %y-velocity

%% Set the initial conditions:
w0Name = @(x) zeros(size(x)); % beam starts at 0 position
tmp = w0Name(wFem.point); w0 = tmp(dof_w);
w(1,:) = w0'; 
y1 = w0'; % y1 = u(t) = 0 beam starts at reference position
y2 = zeros(size(y1)); % y2 = u_t(t) = 0 since no-flow at start (time-ramp)


%% Loop over each time step
for kt = 1:Ntime
    %% Set the current time
    th = (kt - 1/2)*tau; % t_{k - 1/2}
    t(kt + 1) = kt*tau; %t(1) = t_0 = 0, ..., t(m-1) = t_m = t_end
    kt
    
    %% make q_0 vector
    %% note: for linformglobherm1d_fsi to work, nx must equal nEleBeam!!!
    iDerivative = 0; % pass in w(k,:) height of beam at the k-th (prior) time step
    qh_0 = linformglobherm1d_fsi(p_fe, w_fe(kt,:), MeshFluid, MeshBeam, pFem, wFem, iDerivative, nGaussPoint);
    q_0 = qh_0(dof_w);    % is there a way to get the vector straight from the pressure?...
                          % perhaps I should have interpolated the pressure
                          % somehow?
    
    %% Crank-Nicolson time stepping
    y2temp = (A_00/tau + A_22/4)\( -0.5*A_22*(y1 + tau*y2/2 + y1) + q_0 + b_A22 + b_A00 );
    y1 = (y1/tau + y2/2 + y2temp/2)*tau; % y1 is the beam displacement not on boundary (update the mesh)
    y2 = y2temp; % y2 is the beam velocity not on boundary (update the fluid interface by coinciding)
    
    %% Save beam disp in w and velocity in w_t (kt+1 new)
    w(kt,:) = y1'; 
    wt(kt,:) = y2'; 
    
    %% Update Finite Element Solution by enforcing boundary conditions (kt old)
    w_e = zeros(1, length(wFem.point));     
    I = find(nt_w == 0); 
    w_e(I) = essName(wFem.point(I)); % normally time dependent, so in loop (wasted here)
    tmp = w_e; tmp(dof_w) = w(kt,:);  
    w_fe(kt,:) = tmp; % displacement saved in w_fe
    
    wt_e = zeros(1, length(wFem.point));
    wt_e(I) = essName(wFem.point(I)); % (generally time dependent but wasted computation here)
    tmp = wt_e; tmp(dof_w) = wt(kt,:); 
    wt_fe(kt,:) = tmp; % velocity saved in wt_fe
    
    %% Update the fluid Mesh:
    MeshFluid = elastic_mesh_update(MeshFluid, MeshBeam, wFem, dof_u, nt_u, w_fe(kt,:), kt);
    
    %% Update the fluid spaces:
    uFem = fem_generator_2d_lagrange(MeshFluid,uDegree);
    pFem = fem_generator_2d_lagrange(MeshFluid,pDegree);

    %mesh_viewer_2d(MeshFluid, MeshFluid.node, 0, 0, 0)
    %pause
    
    %% Update fluid velocity boundary conditions:
    I = find(nt_u == 0);
    u1essential(I) = essbc1Name(uFem.point(1,I)', uFem.point(2,I)', t(kt+1)); 
    if kt == 1
        u2essential(I) = essbc2Name(uFem.point(1,I), uFem.point(2,I));
    else
        u2essential(I) = essbc2Name_time(uFem.point(1,I), uFem.point(2,I), wt_fe(kt,:), MeshBeam, wFem);
    end
    
    %% Remake the matrix blocks for the fuid problem:
    aName = @(x,y) ones(size(x));  nQuadraturePoint = 7;
    xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0;
    Mhuu_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, MeshFluid, ...
        uFem, xDerivative1,yDerivative1, uFem,xDerivative2,yDerivative2,nQuadraturePoint);
    xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
    Mhuu_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, MeshFluid, ...
        uFem, xDerivative1,yDerivative1, uFem,xDerivative2,yDerivative2,nQuadraturePoint);
    
    xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
    Mhup_hx0 = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, MeshFluid, ...
        uFem, xDerivative1,yDerivative1, pFem ,xDerivative2,yDerivative2,nQuadraturePoint);
    Mhpu_h0x = Mhup_hx0';
    xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 0;
    Mhup_hy0 = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, MeshFluid, ...
        uFem, xDerivative1,yDerivative1, pFem ,xDerivative2,yDerivative2,nQuadraturePoint);
    Mhpu_h0y = Mhup_hy0';
    
    Muu_hxx = Mhuu_hxx(dof_u, dof_u);
    Muu_hyy = Mhuu_hyy(dof_u, dof_u);
    Mup_hx0 = Mhup_hx0(dof_u, dof_p); Mpu_h0x = Mup_hx0';
    Mup_hy0 = Mhup_hy0(dof_u, dof_p); Mpu_h0y = Mup_hy0';
    
    %% Assemble the blocks to form the matrix for the fluid system:
    nu = 1; rho = 1; % fluid parameters
    SP0_u = sparse(length(dof_u), length(dof_u));
    SP0_p = sparse(length(dof_p), length(dof_p));
    M_tmp = [nu*(Muu_hxx + Muu_hyy), SP0_u, -(1/rho)*Mup_hx0; ...
        SP0_u, nu*(Muu_hxx + Muu_hyy), -(1/rho)*Mup_hy0; ...
        -(1/rho)*Mpu_h0x, -(1/rho)*Mpu_h0y, SP0_p ];
    
    
    %% Remake vector blocks for fluid problem
    tmp = -nu*Mhuu_hxx*u1essential; bc_u1_eMhxx = tmp(dof_u);
    tmp = -nu*Mhuu_hyy*u1essential; bc_u1_eMhyy = tmp(dof_u);
    tmp = -nu*Mhuu_hxx*u2essential; bc_u2_eMhxx = tmp(dof_u);
    tmp = -nu*Mhuu_hyy*u2essential; bc_u2_eMhyy = tmp(dof_u);
    tmp = -(1/rho)*Mhpu_h0x*u1essential; bc_u1_eMh0x = tmp(dof_p);
    tmp = -(1/rho)*Mhpu_h0y*u2essential; bc_u2_eMh0y = tmp(dof_p);

    %% Assemble the right hand side vector using the vector blocks for fluid problem:
    rhs = [f1_h + bc_u1_eMhxx + bc_u1_eMhyy; ...
        f2_h + bc_u2_eMhxx + bc_u2_eMhyy; ...
        -bc_u1_eMh0x - bc_u2_eMh0y; ...
        0];
    
    %% Solve the Fluid system:
    u1u2p = M\rhs;
    u1 = u1u2p(1:length(dof_u));
    u1_fe = u1essential; u1_fe(dof_u) = u1; % u1_fe is FE solution in x direction
    u2 = u1u2p(length(dof_u) + 1:(2*length(dof_u)));
    u2_fe = u2essential; u2_fe(dof_u) = u2; % u2_fe is FE solution in x direction
    p = u1u2p((2*length(dof_u))+1:(end-1));
    p_fe = p; % p_fe is FE solution for pressure
end

figure(1)
quiver(uFem.point(1, :), uFem.point(2, :), u1_fe', u2_fe', ...
'Color', [1, 0, 0])
xlabel('x'), ylabel('y'), title('velocity field'), grid on
axis equal; figure(2)
% plot the pressure function
%x = pFem.point(1,:); y = pFem.point(2,:); z = p_fe;
%scatter3(x,y,z,100,z,'.')
%a = -37.5; b = 30; view(a,b)
%colormap(jet), colorbar
%xlabel('x'), ylabel('y'), zlabel('p'), title('pressure field')
figure(3)
mesh_viewer_2d(MeshFluid, MeshFluid.node, 0, 0, 0)

    
    
    
    
    
    
    
    



