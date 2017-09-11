% FEM TEST 2:
% (c) Angelo Marney 
% solution: u(x,y) = e^x e^y
%           u_x = u_xx = u_y = u_yy = e^x e^y
% coeffs:    I = I_x = I_y = I_xx = I_yy = I_z = e^x e^y e^z
%           so I/I_x = I/I_y = I/I_xfemx = I/I_yy = I/I_z = 1
%           k = -2 e^x e^y 
% pde:         
%           (I/I_z) u_xx + (I/I_z) u_yy = -k
% weak: 
%   \int I/Iz ux vx + \int I/Iz uy vy = \int -kv
%   I_xx + I_yy = -k_0

%% Make a data structure for the mesh:
xMin = 0; xMax = 1; yMin = 0; yMax = 1; nx = 10; ny = 10;
Mesh = mesh_generator_2d_triangle(xMin, xMax, yMin, yMax, nx, ny);


%% Add edge data to mesh:
[edge, M_e, M_ne] = mesh_edge_generator(Mesh,1,1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;

%% Make a data structure for finite element shape function access:
iDegree = 2; Fem = fem_generator_2d_lagrange(Mesh, iDegree);
Fem1 = Fem; Fem2 = Fem;

%% prepare dof_u and ess bc:
bName = @(x,y) x.*(x - xMax).*y.*(y - yMax);
bdName = @(x,y) x.*(x - xMax).*y.*(y-  yMax);
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);
g_D = @(x,y) exp(x).*exp(y);
u_e = zeros(length(Fem.point), 1);
I = find(nt == 0); u_e(I) = g_D(Fem.point(1, I), Fem.point(2, I));

%% prepare the matrices:
IName = @(x,y) -1*ones(size(x)); nQuadraturePoint = 7;

xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0;
Ih_xx = stiffness_matrix_assembler_2d_lagrange_tri_global(IName, Mesh, Fem1, ... 
    xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint);
I_xx = Ih_xx(dof_u, dof_u);

xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Ih_yy = stiffness_matrix_assembler_2d_lagrange_tri_global(IName, Mesh, Fem1, ... 
    xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint);
I_yy = Ih_yy(dof_u, dof_u);

%% prepare the -k_0 vector
kName = @(x,y) exp(x).*exp(y) + exp(x).* exp(y); % remove the negative sign, since - is in pde
xDerivative = 0; yDerivative = 0;
kh_0 = load_vector_assembler_2d_lagrange_tri_global(kName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
k_0 = kh_0(dof_u);

%% prepare the b/c vectors associated with the matrices:
tmp = -Ih_xx*u_e; b_I_xx= tmp(dof_u);
tmp = -Ih_yy*u_e; b_I_yy= tmp(dof_u);

%% solve the system:
%  ( Ix_x0 + I_xx + Iy_y0 + I_yy)*u_fe = k_0 + b_Ix_x0 + b_I_xx + b_Iy_y0 +  b_Iy_y0
u = (I_xx + I_yy)\(k_0 + b_I_xx +  b_I_yy);
u_fe = u_e; u_fe(dof_u) = u;

%% error:
uName = @(x,y) exp(x).*exp(y); xDerivative = 0; yDerivative = 0; nQuadraturePoint = 7;
err = sqrt(error_2d_lagrange_global(uName, u_fe, Mesh, Fem, xDerivative, yDerivative, nQuadraturePoint))

pt.x = 0.5; pt.y = 0.5;
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


