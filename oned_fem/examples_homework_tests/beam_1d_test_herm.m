%beam_1d_timeless
% cantilevered beam: u(0) = 0, u'(0) = pi
%                    u''(1) =0, u'''(1) = pi^3
% solution: u(x) = sin(pi*x);


% make mesh
domain = [0,1];
nElement = 4;
Mesh = mesh_generator_1d(domain, nElement);

% make fem
   iDegree = 3;
   Fem = femherm1d(Mesh, iDegree);
   hDegree = 1; % h degree matches Fem.degree
%hDegree = 1;
%Fem = fem_generator_1d_lagrange(Mesh, hDegree); % equivalent to femherm1d

% make dof_u and nt
iBoundaryCondition = [0,-1,0,-1];
[dof_u, nt] = bc_array_generator_1d(Fem, domain, iBoundaryCondition);
Fem.dof_u = dof_u; Fem.nt = nt;

% make A matrix
aName = @(x) ones(size(x));
iDerivative = 2;
nGaussPoint = 5;
Ah_22 = biformglobherm1d(aName,...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
A_22 = Ah_22(dof_u, dof_u);

% make q vector
qName = @(x) pi^4 *sin(pi*x);
iDerivative = 0;
nGaussPoint = 5;
qh_0 = linformglobherm1d(qName, ...
    Mesh, Fem, iDerivative, nGaussPoint);
q_0 = qh_0(dof_u);

% make b_A22 vector
essName = @(x) [0, pi];
I = find(nt == 0); 
u_e = zeros(length(Fem.point), 1);
u_e(I) = essName(Fem.point(I));
tmp = -Ah_22*u_e;
b_A22 = tmp(dof_u);

% make b_n1 and b_n2 vectors
natName = @(x) [pi^3, 0];
u_n = zeros(length(Fem.point), 1);
I = find(nt == -1); 
u_n(I) = natName(Fem.point(I));
b_n = u_n(dof_u);

% Solve the system of equations and assemble u_fe = uGlobal
uh_fe = A_22\(q_0 + b_A22 - b_n);
u_fe = u_e;
u_fe(dof_u) = uh_fe;

% evaluate FE soution at 1
tx = 1;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    if (tx - element(1))*(tx - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        ufe_tx1 = evalfeherm1d(tx, u_fe_loc, element, hDegree, 0);
        break
    end
end

figure(1), hold on
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    uLocal = u_fe(Fem.T(:,k));
    x = linspace(element(1), element(2), 11);
    plot(x, evalfeherm1d(x, uLocal, element, hDegree, 0));
end

x = linspace(0,1,nElement*11);
plot(x,sin(pi*x))

uName = @(x) sin(pi*x);
derivOrder = 0; nGaussPoint = 5;
error = errorglobherm1d(uName, u_fe, Mesh, Fem, derivOrder, nGaussPoint)

