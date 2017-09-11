%% STOKES FLUID Q32

% prepare solution domain, mesh, and fem
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 5; ny = 5;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,ny);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh,1,1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;

% TAYLOR HOOD PAIR:
uDegree = 2; pDegree = 1;
uFem = fem_generator_2d_lagrange(Mesh,uDegree);
pFem = fem_generator_2d_lagrange(Mesh,pDegree);

% prepare dof_u and boundary conditions
bName = @(x,y) (x-xmin).*(x-xmax).*(y-ymin).*(y-ymax);
bduName = @(x,y) (x-xmin).*(x-xmax).*(y-ymin).*(y-ymax);
bdpName = @(x,y) ones(size(x));
[dof_u_u, nt_u] = bc_array_generator_2d(uFem, bName, bduName);
[dof_u_p, nt_p] = bc_array_generator_2d(pFem, bName, bdpName);

g1_D = @(x,y) ...
    1.*((x >= -1) & (x <= 1) & (y == 1)); %...
    %+ 0.*((x < -1) | (x > 1) | (y ~= 1));
g2_D = @(x,y) zeros(size(x));
u1_e = zeros(length(uFem.point),1);
u2_e = zeros(length(uFem.point),1);
I = find(nt_u == 0);
u1_e(I) = g1_D(uFem.point(1,I)', uFem.point(2,I)');
u2_e(I) = g2_D(uFem.point(1,I), uFem.point(2,I));

% prepare the matrix blocks
aName = @(x,y) ones(size(x));  nQuadraturePoint = 7;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0;
Mhuu_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, ...
    uFem, xDerivative1,yDerivative1, uFem,xDerivative2,yDerivative2,nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Mhuu_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, ...
    uFem, xDerivative1,yDerivative1, uFem,xDerivative2,yDerivative2,nQuadraturePoint);

xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Mhup_hx0 = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, ...
    uFem, xDerivative1,yDerivative1, pFem ,xDerivative2,yDerivative2,nQuadraturePoint);
Mhpu_h0x = Mhup_hx0';
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 0;
Mhup_hy0 = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, ...
    uFem, xDerivative1,yDerivative1, pFem ,xDerivative2,yDerivative2,nQuadraturePoint);
Mhpu_h0y = Mhup_hy0';

Muu_hxx = Mhuu_hxx(dof_u_u, dof_u_u);
Muu_hyy = Mhuu_hyy(dof_u_u, dof_u_u);
Mup_hx0 = Mhup_hx0(dof_u_u, dof_u_p); Mpu_h0x = Mup_hx0';
Mup_hy0 = Mhup_hy0(dof_u_u, dof_u_p); Mpu_h0y = Mup_hy0';

% Form the matrix for the FE system:
nu = 1; rho = 1;
SP0_u = sparse(length(dof_u_u), length(dof_u_u));
SP0_p = sparse(length(dof_u_p), length(dof_u_p));
M_tmp = [nu*(Muu_hxx + Muu_hyy), SP0_u, -(1/rho)*Mup_hx0; ...
SP0_u, nu*(Muu_hxx + Muu_hyy), -(1/rho)*Mup_hy0; ...
-(1/rho)*Mpu_h0x, -(1/rho)*Mpu_h0y, SP0_p ];


% append gamma = 1 vector
fName = @(x,y) ones(size(x));
xDerivative = 0; yDerivative = 0;
fh_h = load_vector_assembler_2d_lagrange_tri_global(fName, Mesh, pFem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
one_vec = fh_h(dof_u_p);

tmp_vec = [zeros(2*length(dof_u_u),1); one_vec];
M = [M_tmp, tmp_vec; tmp_vec', 0];

% Prepare vector blocks:
f1Name = @(x,y) zeros(size(x));  % SHOULD BE ZEROS?
f2Name = @(x,y) zeros(size(x));  % SHOULD BE ZEROS?
xDerivative = 0; yDerivative = 0;
f1h_h = load_vector_assembler_2d_lagrange_tri_global(f1Name, Mesh, uFem, xDerivative, ...
    yDerivative, nQuadraturePoint);
f1_h = f1h_h(dof_u_u);
f2h_h = load_vector_assembler_2d_lagrange_tri_global(f2Name, Mesh, uFem, xDerivative, ...
    yDerivative, nQuadraturePoint);
f2_h = f2h_h(dof_u_u);

tmp = -nu*Mhuu_hxx*u1_e; bc_u1_eMhxx = tmp(dof_u_u);
tmp = -nu*Mhuu_hyy*u1_e; bc_u1_eMhyy = tmp(dof_u_u);
tmp = -nu*Mhuu_hxx*u2_e; bc_u2_eMhxx = tmp(dof_u_u);
tmp = -nu*Mhuu_hyy*u2_e; bc_u2_eMhyy = tmp(dof_u_u);
tmp = -(1/rho)*Mhpu_h0x*u1_e; bc_u1_eMh0x = tmp(dof_u_p);
tmp = -(1/rho)*Mhpu_h0y*u2_e; bc_u2_eMh0y = tmp(dof_u_p);

rhs = [f1_h + bc_u1_eMhxx + bc_u1_eMhyy; ...
    f2_h + bc_u2_eMhxx + bc_u2_eMhyy; ...
    -bc_u1_eMh0x - bc_u2_eMh0y; ...
    0];

% Solve the FE system:
u1u2p = M\rhs;
u1 = u1u2p(1:length(dof_u_u));
u1_fe = u1_e; u1_fe(dof_u_u) = u1;
u2 = u1u2p(length(dof_u_u) + 1:(2*length(dof_u_u)));
u2_fe = u2_e; u2_fe(dof_u_u) = u2;
p = u1u2p((2*length(dof_u_u))+1:(end-1));
p_fe = p;

% visualize FE solution: plot of velocity field
figure(1)
quiver(uFem.point(1, :), uFem.point(2, :), u1_fe', u2_fe', ...
'Color', [1, 0, 0])
xlabel('x'), ylabel('y'), title('velocity field')
axis equal; axis([xmin xmax ymin ymax]); figure(2)
% plot the pressure function
x = pFem.point(1,:); y = pFem.point(2,:); z = p_fe;
scatter3(x,y,z,100,z,'.')
a = -37.5; b = 30; view(a,b)
colormap(jet), colorbar
xlabel('x'), ylabel('y'), zlabel('p'), title('pressure field')

% evaluate the pressure at at the points (-.9,.9), (.9,.9)
% pt.x = -.9; pt.y = .9;
% for k = 1:size(Mesh.element,2)
%     element = Mesh.node(:, Mesh.element(:, k));
%     point_is_in = ptInTriangle(pt, element);
%     if (point_is_in)
%         break
%     end
% end
% pfeLocal = p_fe(pFem.T(:,k));
% xDerivative1 = 0; yDerivative1 = 0;
% pfe_at_pt1 = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, pfeLocal, ...
%     element, pDegree, xDerivative1, yDerivative1);
% 
% pt.x = .9; pt.y = .9;
% for k = 1:size(Mesh.element,2)
%     element = Mesh.node(:, Mesh.element(:, k));
%     point_is_in = ptInTriangle(pt, element);
%     if (point_is_in)
%         break
%     end
% end
% pfeLocal = p_fe(pFem.T(:,k));
% xDerivative1 = 0; yDerivative1 = 0;
% pfe_at_pt2 = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, pfeLocal, ...
%     element, pDegree, xDerivative1, yDerivative1);
% 
% flag = 0; tol = 10e-3;
% for k1 = -0.2:0.001:0.2
%     for k2 = 0.2:0.001:0.9
%         pt.x = k1; pt.y = k2;
%         for k = 1:size(Mesh.element,2)
%             element = Mesh.node(:, Mesh.element(:, k));
%             point_is_in = ptInTriangle(pt, element);
%             if (point_is_in)
%                 break
%             end
%         end
%         u1_fe_loc = u1_fe(uFem.T(:,k));
%         check_1 = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u1_fe_loc, ...
%             element, uDegree, 0, 0);
%         u2_fe_loc = u2_fe(uFem.T(:,k));
%         check_2 = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, u2_fe_loc, ...
%             element, uDegree, 0, 0);        
%         if ((abs(check_1) < tol) && (abs(check_2) < tol))
%             flag = 1;
%             break
%         end
%     end
%     if (flag == 1)
%         break
%     end
% end
% 
% 
% 


