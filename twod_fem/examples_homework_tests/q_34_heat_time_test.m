%q 34 time-stepping heat equation
% prepare the mesh and fem space
xmin = 0; xmax = 1; ymin = 0; ymax = 1; nx = 40; ny = 40;
Mesh = mesh_generator_2d_triangle(xmin,xmax,ymin,ymax,nx,ny);
[edge, M_e, M_ne] = mesh_edge_generator(Mesh,1,1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
iDegree = 1; Fem = fem_generator_2d_lagrange(Mesh,iDegree);


% set up time-stepping and initial condition
tspan = [0,1]; tau = 1/40; Ntime = (tspan(2) - tspan(1))/tau;

% prepare dof_u and ess bc and function handles
bName = @(x,y) x.*(x-xmax).*y.*(y-ymax);
bdName = @(x,y) x.*(x-xmax).*y.*(y-ymax);
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);
Fem.dof_u = dof_u; Fem.nt = nt;

aName = @(x,y) 1 + x + cos(y);
fName = @(x,y,t) cos(x + y/2 + t).*(5*(2+x)/4 + 3*cos(y)/4) + ...
    0.5*cos(x + 3*y/2 + t);
gD_name = @(x,y,t) cos(x).*cos(y/2 + t) - sin(x).*cos(y/2).*sin(t) + ...
    sin(x).*sin(y/2).*(x.*y - x - y + 1 - cos(t));
gDdt_name = @(x,y,t) -cos(x).*sin(y/2 + t) - sin(x).*cos(y/2).*cos(t) + ...
    sin(x).*sin(y/2).*sin(t); % time derivative of dirch
u0_name = @(x,y) cos(x + y./2); % u0 - initial condition

% prepare the matrices and initial condition
nQuadraturePoint = 3; 

one_fun = @(x,y) ones(size(x));
xDerivative1 = 0; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 0;
Mh_h00 = stiffness_matrix_assembler_2d_lagrange_tri_global(one_fun, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);

xDerivative1 = 1; xDerivative2 = 1;% using aName 
Mh_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);

xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1; % using aName
Mh_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName, Mesh, Fem, ... 
    xDerivative1, yDerivative1, Fem, xDerivative2, yDerivative2, nQuadraturePoint);

M_h00 = Mh_h00(dof_u, dof_u); 
M_hxx = Mh_hxx(dof_u, dof_u); M_hyy = Mh_hyy(dof_u, dof_u);

M_left = (M_h00 + (tau/2)*(M_hxx + M_hyy)); 
M_right = (M_h00 - (tau/2)*(M_hxx + M_hyy)); 

t = zeros(Ntime + 1, 1); t(1) = tspan(1);
u = zeros(Ntime + 1, length(dof_u));
tmp = u0_name(Fem.point(1,:),Fem.point(2,:)); u0 = tmp(dof_u)';
u(1,:) = u0';

%Solve for the unknown at each time step:
for k = 1:Ntime
    th = (k-1/2)*tau; t(k+1) = k*tau;
    
    xDerivative = 0; yDerivative = 0;
    fh_h = load_vector_assembler_2d_lagrange_tri_global_t(fName, Mesh, Fem, xDerivative, ... 
        yDerivative, nQuadraturePoint,th);
    f_h = fh_h(dof_u);
    
    I = find(nt == 0);
    u_e = zeros(length(Fem.point),1); u_e(I) = gD_name(Fem.point(1,I), Fem.point(2,I), th);
    u_edt = zeros(length(Fem.point),1); u_edt(I) = gDdt_name(Fem.point(1,I), Fem.point(2,I), th);
    
    tmp = -Mh_h00*u_edt; bc_eMh00 = tmp(dof_u);
    tmp = -Mh_hxx*u_e; bc_eMhxx = tmp(dof_u);
    tmp = -Mh_hyy*u_e; bc_eMhyy = tmp(dof_u);
    
    u_tmp = M_right*u0 + tau*(f_h + bc_eMh00 + bc_eMhxx + bc_eMhyy);
    u0 = M_left\u_tmp;
    u(k+1,:) = u0';
end

% Form the FE solution:
u_fe = zeros(length(t), length(Fem.point));
I = find(Fem.nt == 0);
for k = 1:length(t)
    u_e = zeros(1,length(Fem.point));
    u_e(I) = gD_name(Fem.point(1,I), Fem.point(2,I), th);
    tmp = u_e; tmp(dof_u) = u(k,:);
    u_fe(k,:) = tmp;
end

%%% POST PROCESS
% check pi/4, pi/6, 0.5
% pt.x = pi/4; pt.y = pi/6; pt.t = 0;
% for k = 1:size(Mesh.element,2)
%     element = Mesh.node(:, Mesh.element(:, k));
%     point_is_in = ptInTriangle(pt, element);
%     if (point_is_in)
%         break;
%     end
% end 
% time_index = find(t == pt.t);
% u_time0 = u_fe(time_index,:);
% ufeLocal = u_time0(Fem.T(:,k));
% ufe_at_pt0 = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, ufeLocal, ...
%     element, iDegree, 0, 0); %xDeriv = 0, yDeriv = 0

% 
% 
% pt.x = pi/4; pt.y = pi/6; pt.t = 0;
% for k = 1:size(Mesh.element,2)
%     element = Mesh.node(:, Mesh.element(:, k));
%     point_is_in = ptInTriangle(pt, element);
%     if (point_is_in)
%         break;
%     end
% end 
% time_index = find(t == pt.t);
% u_time2 = u_fe(time_index,:);
% ufeLocal = u_time2(Fem.T(:,k));
% ufe_at_pt2 = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, ufeLocal, ...
%     element, iDegree, 0, 0); %xDeriv = 0, yDeriv = 0


pt.x = pi/4; pt.y = pi/6; pt.t = 0.5;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:, Mesh.element(:, k));
    point_is_in = ptInTriangle(pt, element);
    if (point_is_in)
        break;
    end
end %ufeGlobal = M\v = u
% element = correct element, it contains pt.x and pt.y
% now need to find time slice index
time_index = find(t == pt.t);
u_time1 = u_fe(time_index,:);
ufeLocal = u_time1(Fem.T(:,k));
ufe_at_pt1 = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, ufeLocal, ...
    element, iDegree, 0, 0); %xDeriv = 0, yDeriv = 0

pt.x = pi/4; pt.y = pi/6; pt.t = 1;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:, Mesh.element(:, k));
    point_is_in = ptInTriangle(pt, element);
    if (point_is_in)
        break;
    end
end 
time_index = find(t == pt.t);
u_time2 = u_fe(time_index,:);
ufeLocal = u_time2(Fem.T(:,k));
ufe_at_pt2 = evaluate_fe_function_2d_lagrange_tri(pt.x, pt.y, ufeLocal, ...
    element, iDegree, 0, 0); %xDeriv = 0, yDeriv = 0

% figure(1);
% trisurf(Fem.T', Fem.point(1,:), Fem.point(2,:), u_time0, ...
%     'facecolor','interp','EdgeColor','none')

figure(2);
trisurf(Fem.T', Fem.point(1,:), Fem.point(2,:), u_time1, ...
    'facecolor','interp','EdgeColor','none')
%zlim([-40,75]), 
view(-136.5,36)
xlabel('x'), ylabel('y'), zlabel('u(X,0.5)')
title('finite element solution at t = 0.5')

figure(3); 
trisurf(Fem.T', Fem.point(1,:), Fem.point(2,:), u_time2, ...
    'facecolor','interp','EdgeColor','none')
%zlim([-40,75]), view(-26.5,45)
view(-136.5,36)
xlabel('x'), ylabel('y'), zlabel('u(X,1)')
title('finite element solution at t = 1')






