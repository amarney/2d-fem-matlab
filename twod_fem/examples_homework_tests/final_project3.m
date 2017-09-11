%% NUMBER 31:
% c
% Prepare the Mesh and Fem structures
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 40; ny = 40;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,nx);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 1; Fem = fem_generator_2d_lagrange(Mesh, iDegree);

% Prepare dof_u and the essential BC
boundName = @(x,y) x.*(x-xmax).*y.*(y-ymax);
boundDirName = @(x,y) (x).*(x-xmax).*(y).*(y-ymax); % all dirichlet
[dof_u, nt] = bc_array_generator_2d(Fem, boundName, boundDirName);
g_D = @(x,y) x.^2 + y.^3 .* (x - sin(pi*x)) + x.* y.* (-54 + 19*y.^2);
u_e = zeros(length(Fem.point),1);
I = find(nt == 0); u_e(I) = g_D(Fem.point(1,I), Fem.point(2,I));

% Prepare the matrices
aName = @(x,y) 1 + x + y./2;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0; nQuadraturePoint = 7;
Ah_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Ah_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
A_hxx = Ah_hxx(dof_u, dof_u); A_hyy = Ah_hyy(dof_u, dof_u);

bName = @(x,y) x/9;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0; nQuadraturePoint = 7;
Bh_hx0 = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 0; nQuadraturePoint = 7;
Bh_hy0 = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
B_hx0 = Bh_hx0(dof_u, dof_u); B_hy0 = Bh_hy0(dof_u, dof_u);

cName = @(x,y) 2*ones(size(x));
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Ch_h00 = stiffness_matrix_assembler_2d_lagrange_tri_global(cName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
C_h00 = Ch_h00(dof_u, dof_u);

% prepare b/c vectors associated with matrices
tmp = -Ah_hxx*u_e; bc_eAxx = tmp(dof_u);
tmp = -Ah_hyy*u_e; bc_eAyy = tmp(dof_u);
tmp = -Bh_hx0*u_e; bc_eBx0 = tmp(dof_u);
tmp = -Bh_hy0*u_e; bc_eBy0 = tmp(dof_u);
tmp = -Ch_h00*u_e; bc_eC00 = tmp(dof_u);

% prepare the load vectors
fName = @(x,y) (1/9)*(-18 + x.^2 .*(20 + 3*y.^2) + x.*y.*(-54 + 19*y.^2));
xDerivative = 0; yDerivative = 0;
fh_h0 = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
f_h0 = fh_h0(dof_u);

% solve the fe system to obtain the fe solution
u = (A_hxx + A_hyy + B_hx0 + B_hy0 + C_h00)\(f_h0 + bc_eAxx + bc_eAyy + bc_eBx0 + bc_eBy0 + bc_eC00);
u_fe = u_e; u_fe(dof_u) = u;

% post process
pt.x = pi/4; pt.y = pi/6;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:, Mesh.element(:, k));
    point_is_in = ptInTriangle(pt, element);
    if (point_is_in)
        break;
    end
end 
u_fe_loc = u_fe(Fem.T(:,k));
u_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 0, 0); %xDeriv = 0, yDeriv = 0
ux_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 1, 0); %xDeriv = 1, yDeriv = 0
uy_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
element, iDegree, 0, 1); %xDeriv = 0, yDeriv = 1







%% PART D
% Prepare the Mesh and Fem structures
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 40; ny = 40;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,nx);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 2; Fem = fem_generator_2d_lagrange(Mesh, iDegree);

% Prepare dof_u and the essential BC
boundName = @(x,y) x.*(x-xmax).*y.*(y-ymax);
boundDirName = @(x,y) (x).*(x-xmax).*(y).*(y-ymax); % all dirichlet
[dof_u, nt] = bc_array_generator_2d(Fem, boundName, boundDirName);
g_D = @(x,y) x.^2 + y.^3 .* (x - sin(pi*x)) + x.* y.* (-54 + 19*y.^2);
u_e = zeros(length(Fem.point),1);
I = find(nt == 0); u_e(I) = g_D(Fem.point(1,I), Fem.point(2,I));

% Prepare the matrices
aName = @(x,y) 1 + x + y./2;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0; nQuadraturePoint = 7;
Ah_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Ah_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
A_hxx = Ah_hxx(dof_u, dof_u); A_hyy = Ah_hyy(dof_u, dof_u);

bName = @(x,y) x/9;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0; nQuadraturePoint = 7;
Bh_hx0 = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 0; nQuadraturePoint = 7;
Bh_hy0 = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
B_hx0 = Bh_hx0(dof_u, dof_u); B_hy0 = Bh_hy0(dof_u, dof_u);

cName = @(x,y) 2*ones(size(x));
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Ch_h00 = stiffness_matrix_assembler_2d_lagrange_tri_global(cName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
C_h00 = Ch_h00(dof_u, dof_u);

% prepare b/c vectors associated with matrices
tmp = -Ah_hxx*u_e; bc_eAxx = tmp(dof_u);
tmp = -Ah_hyy*u_e; bc_eAyy = tmp(dof_u);
tmp = -Bh_hx0*u_e; bc_eBx0 = tmp(dof_u);
tmp = -Bh_hy0*u_e; bc_eBy0 = tmp(dof_u);
tmp = -Ch_h00*u_e; bc_eC00 = tmp(dof_u);

% prepare the load vectors
fName = @(x,y) (1/9)*(-18 + x.^2 .*(20 + 3*y.^2) + x.*y.*(-54 + 19*y.^2));
xDerivative = 0; yDerivative = 0;
fh_h0 = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
f_h0 = fh_h0(dof_u);

% solve the fe system to obtain the fe solution
u = (A_hxx + A_hyy + B_hx0 + B_hy0 + C_h00)\(f_h0 + bc_eAxx + bc_eAyy + bc_eBx0 + bc_eBy0 + bc_eC00);
u_fe = u_e; u_fe(dof_u) = u;

% post process
pt.x = pi/4; pt.y = pi/6;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:, Mesh.element(:, k));
    point_is_in = ptInTriangle(pt, element);
    if (point_is_in)
        break;
    end
end 

u_fe_loc = u_fe(Fem.T(:,k));
u_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 0, 0); %xDeriv = 0, yDeriv = 0
ux_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 1, 0); %xDeriv = 1, yDeriv = 0
uy_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
element, iDegree, 0, 1); %xDeriv = 0, yDeriv = 1

