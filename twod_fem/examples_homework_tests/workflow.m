%   FE TEST
% prepare the mesh and fem space
xmin = 0; xmax = 1; ymin = 0; ymax = 1.5; nx = 20; ny = 20;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,ny);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh,1,1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 2; Fem = fem_generator_2d_lagrange(Mesh,iDegree);
Fem1 = Fem; Fem2 = Fem;

%prepare dof_u and ess bc
bName = @(x,y) x.*(x-xmax).*y.*(y-ymax);
bdName = @(x,y) x.*(x-xmax).*y.*(y-ymax);
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);
g_D = @(x, y) exp(x).*exp(y-1);
u_e = zeros(length(Fem.point), 1);
I = find(nt == 0); u_e(I) = g_D(Fem.point(1, I), Fem.point(2, I));

% prepare the matrices
aName = @(x, y) 1 + cos(x).*exp(y);
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0; nQuadraturePoint = 7;
Sh_hx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem1, ... 
    xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Sh_hy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem1, ... 
    xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint);
S_hx = Sh_hx(dof_u, dof_u); S_hy = Sh_hy(dof_u, dof_u);

cName = @(x, y) 1 + sin(x).*exp(y);
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Mh_h = stiffness_matrix_assembler_2d_lagrange_tri_global(cName, Mesh, Fem1, ... 
    xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint);
M_h = Mh_h(dof_u, dof_u);

% prepare b/c vectors associated with matrices
tmp = -Sh_hx*u_e; bc_eSx = tmp(dof_u);
tmp = -Sh_hy*u_e; bc_eSy = tmp(dof_u);
tmp = -Mh_h*u_e; bc_eM = tmp(dof_u);

% prepare the load vectors
fName = @(x, y) exp(x + y - 1).*(-1 - 3*exp(y).*cos(x) ...
+ 2*exp(y).*sin(x));
xDerivative = 0; yDerivative = 0;
fh_h = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
f_h = fh_h(dof_u);

% solve the fe system to obtain the fe solution
u = (S_hx + S_hy + M_h)\(f_h + bc_eSx + bc_eSy + bc_eM);
u_fe = u_e; u_fe(dof_u) = u;

% post process
u = @(x, y) exp(x).*exp(y-1);
u_dx = @(x, y) exp(x).*exp(y-1);
u_dy = @(x, y) exp(x).*exp(y-1);

err_l2_20 = sqrt(error_2d_lagrange_global(u, u_fe, Mesh, Fem, 0, 0, nQuadraturePoint));
err_H1_x = error_2d_lagrange_global(u_dx, u_fe, Mesh, Fem, 1, 0, nQuadraturePoint);
err_H1_y = error_2d_lagrange_global(u_dy, u_fe, Mesh, Fem, 0, 1, nQuadraturePoint);
err_H1_20 = sqrt(err_H1_x + err_H1_y);






