% EX of local load vector usage
nQuadraturePoint = 7; element = [0 0.1 0.05; 1 1.1 1.2];
fName = @(x,y) sin(x) + cos(y);
iDegree = 1; xDerivative = 0; yDerivative = 0;
localLoadVector = load_vector_assembler_2d_lagrange_tri_local(fName, element, ...
    iDegree, xDerivative, yDerivative, nQuadraturePoint);

% EX of local stiffness matrix usage
nQuadraturePoint = 7; element = [0 0.1 0.05; 1 1.1 1.2];
aName = @(x,y) sin(x) + cos(y);

iDegree1 = 1; iDegree2 = 1;
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
M_1 = stiffness_matrix_assembler_2d_lagrange_tri_local(aName, element, iDegree1, xDerivative1, ...
    yDerivative1, iDegree2, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0;
M_2 = stiffness_matrix_assembler_2d_lagrange_tri_local(aName, element, iDegree1, xDerivative1, ...
    yDerivative1, iDegree2, xDerivative2, yDerivative2, nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
M_3 = stiffness_matrix_assembler_2d_lagrange_tri_local(aName, element, iDegree1, xDerivative1, ...
    yDerivative1, iDegree2, xDerivative2, yDerivative2, nQuadraturePoint);

clear
clc
% Q 26 part c
% i)
% - - - 
xMin = 0; xMax = 1; yMin = 0; yMax = 1; nx = 20; ny = 20;
Mesh = mesh_generator_2d_triangle(xMin, xMax, yMin, yMax, nx, ny);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1); 
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 1; Fem = fem_generator_2d_lagrange(Mesh, iDegree); 
fName = @(x,y) cos(x + 2*y); aName = @(x,y) 1;
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
nQuadraturePoint = 7;
M = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
v = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, Fem, ...
    xDerivative1, yDerivative1, nQuadraturePoint); %xderiv1, yderiv1 = 0
ufeGlobal = M\v; % L2 projection
realans = cos(Fem.point(1,:) + 2*Fem.point(2,:));
pt.x = pi/4; pt.y = pi/6;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:, Mesh.element(:, k));
    point_is_in = ptInTriangle(pt, element);
    if (point_is_in)
        break;
    end
end
ufeLocal = ufeGlobal(Fem.T(:,k));
ufe_at_pt1 = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, ufeLocal, ...
    element, iDegree, xDerivative1, yDerivative1); %xderiv1,yderiv1 = 0
% - - -

% ii)
xMin = 0; xMax = 1; yMin = 0; yMax = 1; nx = 20; ny = 20;
Mesh = mesh_generator_2d_triangle(xMin, xMax, yMin, yMax, nx, ny);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1); 
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 2; Fem = fem_generator_2d_lagrange(Mesh, iDegree); 
fName = @(x,y) cos(x + 2*y); aName = @(x,y) 1;
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
nQuadraturePoint = 7;
M = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);
v = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, Fem, ...
    xDerivative1, yDerivative1, nQuadraturePoint); %xderiv1, yderiv1 = 0
ufeGlobal = M\v; % L2 projection
pt.x = pi/4; pt.y = pi/6;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:, Mesh.element(:, k));
    point_is_in = ptInTriangle(pt, element);
    if (point_is_in)
        break;
    end
end
ufeLocal = ufeGlobal(Fem.T(:,k));
ufe_at_pt2 = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, ufeLocal, ...
    element, iDegree, xDerivative1, yDerivative1); %xderiv1,yderiv1 = 0
 %T = table(realans', L2Projection2);
 realans = cos(pi/4 + 2*pi/6);
 
 
 % iii)
xMin = 0; xMax = 1; yMin = 0; yMax = 1; %nx = 20; ny = 20;
uName = @(x,y) cos(x + 2*y); aName = @(x,y) 1;
udxName = @(x,y) - sin(x + 2*y); udyName = @(x,y) -2*sin(x + 2*y);
iDegree = 2; nQuadraturePoint =  7; xDerivative = 0; yDerivative = 0;
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
err1 = zeros(10,1);
for nx = 10:10:100
    ny = nx;
    Mesh = mesh_generator_2d_triangle(xMin, xMax, yMin, yMax, nx, ny);
    [edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1); 
    Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
    Fem = fem_generator_2d_lagrange(Mesh, iDegree); Fem1 = Fem; Fem2 = Fem;
    M = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem1, ... 
        xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint);
    v = load_vector_assembler_2d_lagrange_tri_global(uName, Mesh, Fem, xDerivative, ... 
        yDerivative, nQuadraturePoint);
    u_fe = M\v; 
    err1(nx/10) = sqrt(error_2d_lagrange_global(uName, u_fe, Mesh, Fem, xDerivative, yDerivative, nQuadraturePoint));
end

iDegree = 2; nQuadraturePoint =  7; 
err2 = zeros(10,1);
for nx = 10:10:100
    ny = nx;
    Mesh = mesh_generator_2d_triangle(xMin, xMax, yMin, yMax, nx, ny);
    [edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1); 
    Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
    Fem = fem_generator_2d_lagrange(Mesh, iDegree); Fem1 = Fem; Fem2 = Fem;
    M = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem1, ... 
        xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint);
    v = load_vector_assembler_2d_lagrange_tri_global(uName, Mesh, Fem, xDerivative, ... 
        yDerivative, nQuadraturePoint);
    u_fe = M\v; 
    err_H1_x = error_2d_lagrange_global(udxName, u_fe, Mesh, Fem, 1, 0, nQuadraturePoint);
    err_H1_y = error_2d_lagrange_global(udyName, u_fe, Mesh, Fem, 0, 1, nQuadraturePoint);
    err2(nx/10) = sqrt(err_H1_x + err_H1_y);
end

% iv
nEle = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
h = 1./nEle;
rc = polyfit(log(h), log(err1)',1);
r_L2 = rc(1); C_L2 = exp(rc(2));

rc = polyfit(log(h), log(err2)',1);
r_H1 = rc(1); C_H1 = exp(rc(2));


% 27
load('mesh_for_test_FE_dof_u_2d_tri.mat')
Mesh.node = mesh.p; Mesh.element = mesh.t;
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1); 
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 2; Fem = fem_generator_2d_lagrange(Mesh, iDegree);
bName = @(x,y) x.*(x-1).*y.*(y-1); %level set function for domain
bdName = @(x,y) (x-1).*y; %dirchlet level set
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);

nDirichlet = sum(nt == 0);  %known
nNeumann = sum(nt == -1);   %unknown (flux known)
nInterior = sum(nt == -200);%unknown
