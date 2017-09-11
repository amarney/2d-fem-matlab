% Q 30
% PART C:
% Prepare the Mesh and Fem structures
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 40; ny = 40;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,nx);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 1; Fem = fem_generator_2d_lagrange(Mesh, iDegree);

% prepare dof_u and the essential BC
bName = @(x,y) x.*(x - xmax).*y.*(y - ymax);
bdName = @(x,y) (x+5).*(x-xmax+5).*(y+5).*(y-ymax+5); % no dirichlet
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);
g_D = @(x,y) zeros(size(x));
u_e = zeros(length(Fem.point),1);
I = find(nt == 0); u_e(I) = g_D(Fem.point(1,I), Fem.point(2,I));

% prepare the matrices
aName = @(x,y) 1 + x + y./2;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0; nQuadraturePoint = 7;
Ah_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Ah_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
A_hxx = Ah_hxx(dof_u, dof_u); A_hyy = Ah_hyy(dof_u, dof_u);

cName = @(x,y) 2*ones(size(x));
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Ch_h00 = stiffness_matrix_assembler_2d_lagrange_tri_global(cName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
C_h00 = Ch_h00(dof_u, dof_u);

% prepare the boundary vectors for the matrices (none)

% prepare the load vector
fName = @(x,y) pi*cos(2*pi*y).*sin(pi*x) + ...
    .5*cos(pi*x).*((4 + 5*pi^2 *(2 + 2*x + y)).*cos(2*pi*y) + 2*pi*sin(2*pi*y));
xDerivative = 0; yDerivative = 0;
fh_h = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
f_h = fh_h(dof_u);

% solve the fe system
u = (A_hxx + A_hyy + C_h00)\(f_h);
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
    element, iDegree, 0, 0);
ux_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 1, 0);
uy_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 0, 1);


% Q 30
% PART C:
% Prepare the Mesh and Fem structures
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 40; ny = 40;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,nx);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 1; Fem = fem_generator_2d_lagrange(Mesh, iDegree);

% prepare dof_u and the essential BC
bName = @(x,y) x.*(x - xmax).*y.*(y - ymax);
bdName = @(x,y) (x+5).*(x-xmax+5).*(y+5).*(y-ymax+5); % no dirichlet
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);
g_D = @(x,y) zeros(size(x));
u_e = zeros(length(Fem.point),1);
I = find(nt == 0); u_e(I) = g_D(Fem.point(1,I), Fem.point(2,I));

% prepare the matrices
aName = @(x,y) 1 + x + y./2;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0; nQuadraturePoint = 7;
Ah_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Ah_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
A_hxx = Ah_hxx(dof_u, dof_u); A_hyy = Ah_hyy(dof_u, dof_u);

cName = @(x,y) 2*ones(size(x));
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Ch_h00 = stiffness_matrix_assembler_2d_lagrange_tri_global(cName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
C_h00 = Ch_h00(dof_u, dof_u);

% prepare the boundary vectors for the matrices (none)

% prepare the load vector
fName = @(x,y) pi*cos(2*pi*y).*sin(pi*x) + ...
    .5*cos(pi*x).*((4 + 5*pi^2 *(2 + 2*x + y)).*cos(2*pi*y) + 2*pi*sin(2*pi*y));
xDerivative = 0; yDerivative = 0;
fh_h = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
f_h = fh_h(dof_u);

% solve the fe system
u = (A_hxx + A_hyy + C_h00)\(f_h);
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
    element, iDegree, 0, 0);
ux_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 1, 0);
uy_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 0, 1);



% Q 30
% PART C:
% Prepare the Mesh and Fem structures
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 40; ny = 40;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,nx);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 1; Fem = fem_generator_2d_lagrange(Mesh, iDegree);

% prepare dof_u and the essential BC
bName = @(x,y) x.*(x - xmax).*y.*(y - ymax);
bdName = @(x,y) (x+5).*(x-xmax+5).*(y+5).*(y-ymax+5); % no dirichlet
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);
g_D = @(x,y) zeros(size(x));
u_e = zeros(length(Fem.point),1);
I = find(nt == 0); u_e(I) = g_D(Fem.point(1,I), Fem.point(2,I));

% prepare the matrices
aName = @(x,y) 1 + x + y./2;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0; nQuadraturePoint = 7;
Ah_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Ah_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
A_hxx = Ah_hxx(dof_u, dof_u); A_hyy = Ah_hyy(dof_u, dof_u);

cName = @(x,y) 2*ones(size(x));
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Ch_h00 = stiffness_matrix_assembler_2d_lagrange_tri_global(cName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
C_h00 = Ch_h00(dof_u, dof_u);

% prepare the boundary vectors for the matrices (none)

% prepare the load vector
fName = @(x,y) pi*cos(2*pi*y).*sin(pi*x) + ...
    .5*cos(pi*x).*((4 + 5*pi^2 *(2 + 2*x + y)).*cos(2*pi*y) + 2*pi*sin(2*pi*y));
xDerivative = 0; yDerivative = 0;
fh_h = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
f_h = fh_h(dof_u);

% solve the fe system
u = (A_hxx + A_hyy + C_h00)\(f_h);
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
    element, iDegree, 0, 0);
ux_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 1, 0);
uy_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 0, 1);

% PART D:
% Prepare the Mesh and Fem structures
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 40; ny = 40;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,nx);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 2; Fem = fem_generator_2d_lagrange(Mesh, iDegree);

% prepare dof_u and the essential BC
bName = @(x,y) x.*(x - xmax).*y.*(y - ymax);
bdName = @(x,y) (x+5).*(x-xmax+5).*(y+5).*(y-ymax+5); % no dirichlet
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);
g_D = @(x,y) zeros(size(x));
u_e = zeros(length(Fem.point),1);
I = find(nt == 0); u_e(I) = g_D(Fem.point(1,I), Fem.point(2,I));

% prepare the matrices
aName = @(x,y) 1 + x + y./2;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0; nQuadraturePoint = 7;
Ah_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Ah_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
A_hxx = Ah_hxx(dof_u, dof_u); A_hyy = Ah_hyy(dof_u, dof_u);

cName = @(x,y) 2*ones(size(x));
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Ch_h00 = stiffness_matrix_assembler_2d_lagrange_tri_global(cName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
C_h00 = Ch_h00(dof_u, dof_u);

% prepare the boundary vectors for the matrices (none)

% prepare the load vector
fName = @(x,y) pi*cos(2*pi*y).*sin(pi*x) + ...
    .5*cos(pi*x).*((4 + 5*pi^2 *(2 + 2*x + y)).*cos(2*pi*y) + 2*pi*sin(2*pi*y));
xDerivative = 0; yDerivative = 0;
fh_h = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
f_h = fh_h(dof_u);

% solve the fe system
u = (A_hxx + A_hyy + C_h00)\(f_h);
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
    element, iDegree, 0, 0);
ux_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 1, 0);
uy_fe_p = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u_fe_loc, ...
    element, iDegree, 0, 1);
