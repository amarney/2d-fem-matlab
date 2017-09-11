%beam_1d_timeless
% sol: u(x) = sin(4*pi*x)
% u'(x) = 4*pi*cos(4*pi*x)

% make mesh
domain = [0,1];
nElement = 200;
Mesh = mesh_generator_1d(domain, nElement);

% make fem
iDegree = 3;
Fem = femherm1d(Mesh, iDegree);
hDegree = 1; 

% make dof_u and nt
iBoundaryCondition = [0,0, 0, 0];
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
qName = @(x) 256*pi^4*sin(4*pi*x);
iDerivative = 0;
nGaussPoint = 5;
qh_0 = linformglobherm1d(qName, ...
    Mesh, Fem, iDerivative, nGaussPoint);
q_0 = qh_0(dof_u);

% make b_A22 vector
essName = @(x) [0, 0, 4*pi, 4*pi];
I = find(nt == 0); %I = I(1:2); % only u, not u'
u_e = zeros(length(Fem.point), 1);
u_e(I) = essName(Fem.point(I));

tmp = -Ah_22*u_e;
b_A22 = tmp(dof_u);

% Solve the system of equations and assemble u_fe = uGlobal
uh_fe = A_22\(q_0 + b_A22);
u_fe = u_e;
u_fe(dof_u) = uh_fe;


% post processing: plot 
figure(1), hold on
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    uLocal = u_fe(Fem.T(:,k));
    x = linspace(element(1), element(2), 11);
    plot(x, evalfeherm1d(x, uLocal, element, hDegree, 0));
end

% post processing: compute error
uName = @(x) sin(4*pi*x);
derivOrder = 0; nGaussPoint = 5;
u0Error = errorglobherm1d(uName, u_fe, Mesh, Fem, derivOrder, nGaussPoint)

