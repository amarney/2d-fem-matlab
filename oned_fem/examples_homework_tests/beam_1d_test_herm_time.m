% solution: u(x,t) = 0.01*t^2 * sin(pi*x)
% problem: u_{xxxx} = u_{tt} + f(x,t)
% u_{xxxx} = pi^4 * (0.01)*t^2 * sin(pi*x)
% u_{t} = 0.02*t*sin(pi*x)
% u_{tt} = 0.02*sin(pi*x)
% so f(x,t) = 0.01 * sin(pi*x) * (pi^4 * t^2 + 2)
% initial conditions: u(x,0) = 0
%                     u_t(x,0) = 0
% boundary conditions: u(0,t) = 0 (dir)
%                      u(1,t) = 0 (dir)
%                      u'(0,t) = (pi*t^2)/100 * cos(pi*0) = (pi*t^2)/100 (dir)
%                      u'(1,t) = -(pi*t^2)/100 (dir)

% make mesh:
domain = [0,1]; nEle = 5;
Mesh = mesh_generator_1d(domain, nEle);

% use C1 elements:
iDegree = 3;
Fem = femherm1d(Mesh, iDegree);
hDegree = 1;

% set bc type (make dof_u and nt):
iBoundaryCondition = [0, 0, 0, 0];
[dof_u, nt] = bc_array_generator_1d(Fem, domain, iBoundaryCondition);
Fem.dof_u = dof_u; Fem.nt = nt;

% make A_22 matrix (no-time)
aName = @(x) ones(size(x));
iDerivative = 2; nGaussPoint = 5;
Ah_22 = biformglobherm1d(aName,...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
A_22 = Ah_22(dof_u, dof_u);

% make A_00 matrix (no-time)
aName = @(x) ones(size(x));
iDerivative = 0; nGaussPoint = 5;
Ah_00 = biformglobherm1d(aName,...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
A_00 = Ah_00(dof_u, dof_u);

% set time:
tspan = [0, 4.5]; tau = 1/40;
Ntime = (tspan(2) - tspan(1))/tau;
t = zeros(Ntime + 1, 1); t(1) = tspan(1);
u = zeros(Ntime + 1, length(dof_u));
fName = @(x,t) pi^4 * 0.01 * sin(pi*x) .* (t.^2 + 2);

% set initial conditions:
u0Name = @(x) zeros(size(x)); % works because u(x,0) = dudx(x,0) = 0;
tmp = u0Name(Fem.point); u0 = tmp(dof_u);
u(1,:) = u0'; % row 1: first timeslice
              % row 2: second timeslice, ...
              % contains unknown nodes (interior pts for ess-problem)
y1 = u0'; % since y1 = u(t)

% normally do ut0name
y2 = zeros(size(y1)); % since y2 = dudt(t), but dudt(0) = 0 in this problem
              %y1 = u(1:length(u0)/2);
              %y2 = u((length(u0)/2 + 1):length(u0));
%uu = [u u];
              
              
% set boundary conditions:
essName = @(x,t) [0, 0, 0.01*pi*t.^2, -0.01*pi*t.^2];
I = find(nt == 0);
u_e = zeros(length(Fem.point), 1);

essDttName = @(x,t) [0, 0, 0.01*pi, -0.01*pi];
u_edtt = zeros(length(Fem.point), 1);

% loop over each time step
for k = 1:Ntime
    th = (k - 1/2)*tau; % t_{k - 1/2}
    t(k + 1) = k*tau; %t(1) = t_0, t(2) = t_1, ... t(m-1) = t_m = t_end
    
    % make f_0 vector
    iDerivative = 0;
    fh_0 = linformglobherm1d_t(fName, Mesh, Fem, iDerivative, nGaussPoint, th);
    f_0 = fh_0(dof_u);
    
    % make b_A22 vector
    u_e(I) = essName(Fem.point(I), th); 
    tmp = -Ah_22*u_e;
    b_A22 = tmp(dof_u);
    
    % make b_A00 vector
    u_edtt(I) = essDttName(Fem.point(I), th);
    tmp = -Ah_00*u_edtt;
    b_A00 = tmp(dof_u);
    
    y2temp = (A_00/tau + A_22/4)\( -0.5*A_22*(y1 + tau*y2/2 + y1) + f_0 + b_A22 + b_A00 );
    y1 = (y1/tau + y2/2 + y2temp/2)*tau;
    y2 = y2temp;
    
    u(k+1,:) = y1';
end

% form the finite element solution:
u_fe = zeros(length(t), length(Fem.point));
I = find(nt == 0);
for k = 1:length(t)
    u_e = zeros(1, length(Fem.point));
    u_e(I) = essName(Fem.point(I), t(k));
    tmp = u_e; tmp(dof_u) = u(k,:);
    u_fe(k,:) = tmp;
end

% post processing:

% set solution for error calculation:
uName = @(x,t)  0.01*t.^2 .* sin(pi*x);
derivOrder = 0; nGaussPoint = 5; time = 91;
error = errorglobherm1d(uName, u_fe(time,:), Mesh, Fem, derivOrder, nGaussPoint, t(time));

% make a movie
figure(3); clf;
jump = 1; frame_count = 1;
for k = 1:jump:length(t)
    plot(Mesh.node, u_fe(k,1:length(Mesh.node)), 'r o', 'markers', 2)
    hold on
    plot(Mesh.node, 0.01*(t(k))^2*sin(pi*Mesh.node), 'b')
    hold off
    axis([-0.1, 1.1, -0.1, 1.1])
    xlabel('x'); ylabel('u(x,t)');
    text(0.8, 1, ['t = ', num2str(t(k))]);
    M(frame_count) = getframe;
    frame_count = frame_count + 1;
end
