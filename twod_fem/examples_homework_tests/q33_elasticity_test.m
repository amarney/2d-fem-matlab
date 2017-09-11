% q33 linear elasticity TEST
% solution: u_1(x,y) = exp(x), u_2(x,y) = exp(y)

% prep solution domain, mesh, and fem
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 10; ny = 10;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,ny);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh,1,1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 2; Fem = fem_generator_2d_lagrange(Mesh,iDegree);

% prepare dof_u and boundary conditions
bName = @(x,y) x.*(x-xmax).*y.*(y-ymax);
bdName = @(x,y) x.*(x-xmax).*y.*(y-ymax);
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);

g1_D = @(x,y) exp(x);
g2_D = @(x,y) exp(y);
u1_e = zeros(length(Fem.point),1);
u2_e = zeros(length(Fem.point),1);
I = find(nt == 0);
u1_e(I) = g1_D(Fem.point(1,I)', Fem.point(2,I)');
u2_e(I) = g2_D(Fem.point(1,I), Fem.point(2,I));

% set parameters
lambda = 2; mu = 3;

% prepare the matrix blocks
% note Mhuu_hyx = Mhuu_hxy??
aName = @(x,y) ones(size(x)); nQuadraturePoint = 7;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0;
Mhuu_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName,Mesh,...
    Fem,xDerivative1,yDerivative1,Fem,xDerivative2,yDerivative2,nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Mhuu_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName,Mesh,...
    Fem,xDerivative1,yDerivative1,Fem,xDerivative2,yDerivative2,nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 1; yDerivative2 = 0;
Mhuu_hyx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName,Mesh,...
    Fem,xDerivative1,yDerivative1,Fem,xDerivative2,yDerivative2,nQuadraturePoint);
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 1;
Mhuu_hxy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName,Mesh,...
    Fem,xDerivative1,yDerivative1,Fem,xDerivative2,yDerivative2,nQuadraturePoint);

Muu_hxx = Mhuu_hxx(dof_u, dof_u);
Muu_hyy = Mhuu_hyy(dof_u, dof_u);
Muu_hyx = Mhuu_hyx(dof_u, dof_u);
Muu_hxy = Mhuu_hxy(dof_u, dof_u);

% form the matrix for the FE system:
%SP0_u = sparse(length(dof_u), length(dof_u));
M = [mu*(Muu_hxx + Muu_hyy) + (lambda + mu)*Muu_hxx, (lambda + mu)*Muu_hyx; ...
    (lambda + mu)*Muu_hxy, mu*(Muu_hxx + Muu_hyy) + (lambda + mu)*Muu_hyy];

% prepare the vector blocks:
f1Name = @(x,y) -8*exp(x);
f2Name = @(x,y) -8*exp(y);
xDerivative = 0; yDerivative = 0;
f1h_h = load_vector_assembler_2d_lagrange_tri_global(f1Name, Mesh, Fem, xDerivative, ...
    yDerivative, nQuadraturePoint);
f1_h = f1h_h(dof_u);
f2h_h = load_vector_assembler_2d_lagrange_tri_global(f2Name, Mesh, Fem, xDerivative, ...
    yDerivative, nQuadraturePoint);
f2_h = f2h_h(dof_u);

% prepare the boundary vector blocks
tmp = -mu*Mhuu_hxx*u1_e; bc_u1_eMhxx_mu = tmp(dof_u);
tmp = -(lambda + mu)*Mhuu_hxx*u1_e; bc_u1_eMhxx_lambda_mu = tmp(dof_u);
tmp = -mu*Mhuu_hyy*u1_e; bc_u1_eMhyy_mu = tmp(dof_u);
tmp = -(lambda + mu)*Mhuu_hyx*u2_e; bc_u2_eMhyx_lambda_mu = tmp(dof_u);

tmp = -mu*Mhuu_hxx*u2_e; bc_u2_eMhxx_mu = tmp(dof_u);
tmp = -mu*Mhuu_hyy*u2_e; bc_u2_eMhyy_mu = tmp(dof_u);
tmp = -(lambda + mu)*Mhuu_hyy*u2_e; bc_u2_eMhyy_lambda_mu = tmp(dof_u);
tmp = -(lambda + mu)*Mhuu_hxy*u1_e; bc_u1_eMhxy_lambda_mu = tmp(dof_u);

rhs = [f1_h + bc_u1_eMhxx_mu + bc_u1_eMhxx_lambda_mu + bc_u1_eMhyy_mu + bc_u2_eMhyx_lambda_mu; ...
    f2_h + bc_u2_eMhxx_mu + bc_u2_eMhyy_mu + bc_u2_eMhyy_lambda_mu + bc_u1_eMhxy_lambda_mu];

% solve the FE system
u1u2 = M\rhs;
u1 = u1u2(1:length(dof_u));
u1_fe = u1_e; u1_fe(dof_u) = u1;
u2 = u1u2(length(dof_u) + 1:(2*length(dof_u)));
u2_fe = u2_e; u2_fe(dof_u) = u2;


% post proces:
pt.x = 1; pt.y = 1;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:, Mesh.element(:, k));
    point_is_in = ptInTriangle(pt, element);
    if (point_is_in)
        break
    end
end
u1_fe_loc = u1_fe(Fem.T(:,k));
u1_fe_pt1 = evaluate_fe_function_2d_lagrange_tri(pt.x,pt.y,u1_fe_loc,element,iDegree,0,0);
u2_fe_loc = u2_fe(Fem.T(:,k));
u2_fe_pt1 = evaluate_fe_function_2d_lagrange_tri(pt.x,pt.y,u2_fe_loc,element,iDegree,0,0);

nQuadraturePoint = 7; 
H1_norm_u1 = 0; H1_norm_u2 = 0;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:, Mesh.element(:,k));
    u1_fe_loc = u1_fe(Fem.T(:,k));
    u2_fe_loc = u2_fe(Fem.T(:,k));
    qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
    qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
    xDerivative = 0; yDerivative = 0; 
    u1 = evaluate_fe_function_2d_lagrange_tri(qPoint(1,:), qPoint(2,:), ...
        u1_fe_loc, element, iDegree, xDerivative, yDerivative);
    u2 = evaluate_fe_function_2d_lagrange_tri(qPoint(1,:), qPoint(2,:), ...
        u2_fe_loc, element, iDegree, xDerivative, yDerivative);    
    xDerivative = 1; yDerivative = 0;
    u1_x = evaluate_fe_function_2d_lagrange_tri(qPoint(1,:), qPoint(2,:), ...
        u1_fe_loc, element, iDegree, xDerivative, yDerivative);
    u2_x = evaluate_fe_function_2d_lagrange_tri(qPoint(1,:), qPoint(2,:), ...
        u2_fe_loc, element, iDegree, xDerivative, yDerivative);        
    xDerivative = 0; yDerivative = 1;
    u1_y = evaluate_fe_function_2d_lagrange_tri(qPoint(1,:), qPoint(2,:), ...
        u1_fe_loc, element, iDegree, xDerivative, yDerivative);
    u2_y = evaluate_fe_function_2d_lagrange_tri(qPoint(1,:), qPoint(2,:), ...
        u2_fe_loc, element, iDegree, xDerivative, yDerivative);        
    integrand1 = u1.^2 + u1_x.^2 + u1_y.^2; integrand2 = u2.^2 + u2_x.^2 + u2_y.^2;
    H1_norm_u1 = H1_norm_u1 + sum(qWeight.*integrand1);
    H1_norm_u2 = H1_norm_u2 + sum(qWeight.*integrand2);
end

u1Name = @(x,y) exp(x); xDerivative = 0; yDerivative = 0; nQuadraturePoint = 7;
errorGlobal = error_2d_lagrange_global(u1Name, u1_fe, Mesh, Fem, xDerivative, yDerivative, nQuadraturePoint)

u2Name = @(x,y) exp(y); xDerivative = 0; yDerivative = 0; nQuadraturePoint = 7;
errorGlobal = error_2d_lagrange_global(u2Name, u2_fe, Mesh, Fem, xDerivative, yDerivative, nQuadraturePoint)

