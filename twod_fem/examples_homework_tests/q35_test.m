%   q 35
% Prepare the Mesh and Fem structures
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 40; ny = 40;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,nx);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 2; Fem = fem_generator_2d_lagrange(Mesh, iDegree);

% Prepare dof_u and the essential BC
bName = @(x,y) x.*(x-xmax).*y.*(y-ymax);
bdName = @(x,y) x.*y;
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);
g_D = @(x,y) 6*ones(size(x));
u_e = zeros(length(Fem.point),1);
I = find(nt == 0); u_e(I) = g_D(Fem.point(1,I), Fem.point(2,I));

 % prepare natural BC
[Kt, et] = elem_edge_types_2D_tri(Mesh, bdName);
Mesh.Kt = Kt; Mesh.et = et;
g_N = @(x,y) zeros(size(x));
% u_n = zeros(length(Fem.point),1);                             not need?
edgeI = find(Mesh.et == -1); I = unique(Mesh.edge(:,edgeI));
u_n(I) = g_N(Fem.point(1,I), Fem.point(2,I));

% prepare the matrices
aName = @(x,y) ones(size(x)); nQuadraturePoint = 7;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0;
Ah_hx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Ah_hy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
A_hx = Ah_hx(dof_u, dof_u); A_hy = Ah_hy(dof_u, dof_u);

cName = @(x,y) 2*ones(size(x));
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Ch_h = stiffness_matrix_assembler_2d_lagrange_tri_global(cName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
C_h = Ch_h(dof_u, dof_u);

% prepare dirichlet b/c vectors associated with matrices
tmp = -Ah_hx*u_e; bc_eAx = tmp(dof_u);
tmp = -Ah_hy*u_e; bc_eAy = tmp(dof_u);
tmp = -Ch_h*u_e; bc_eC = tmp(dof_u);

% prepare neumann b/c vector 
nGaussPoint = 5;
bch_N = neumann_bc_vector_assembler_2d_lagrange_tri_global(g_N, ...
    Mesh, Fem, nGaussPoint);
bc_N = bch_N(dof_u);

% prepare the load vector
fName = @(x,y) 12*ones(size(x));
xDerivative = 0; yDerivative = 0;
fh_h = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
f_h = fh_h(dof_u);

% Solve the fe system to obtain the fe solution
u = (A_hx + A_hy + C_h)\(f_h + bc_N + bc_eAx + bc_eAy + bc_eC);
u_fe = u_e; u_fe(dof_u) = u;

pt.x = 0.51; pt.y = 0.51;
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



