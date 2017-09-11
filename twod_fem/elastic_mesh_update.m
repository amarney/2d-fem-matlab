function MeshNew = elastic_mesh_update(MeshFluid, MeshBeam, wFem, dof_u, nt_u, beamDisplacement, kt)
%% ELASTICITY MESH UPDATER:
%% Make mFem:
mDegree = 2; mFem = fem_generator_2d_lagrange(MeshFluid,mDegree); % match uFem

%% bottom boundary:
% beamDisplacement = w_fe(kt,:) (wGlobal)

% wInterpolate is the bottom boundary

% all i need to do is evalute beamDisplacement at mFem.points
%               well, the mFem.points on the bottom boundary
% [I,J] = find(M_e == 1) -> gives nodes of boundary edges
% e55 = M_ne(I(5), J(5)) -> gives edge connecting I(5) and J(5)
% [I(5); J(5)] == edge(:, e55)
% [MeshFluid.node(:,I(5)), MeshFluid.node(:,J(5))] -> coordinates of boundary nodes
% [MeshFluid.node(:,I); MeshFluid.node(:,J)] -> coordinates of boundary
%

%% Locate the boundary nodes to make dof_m, nt_m arrays
% assume old information carries over
dof_m = dof_u; nt_m = nt_u;

%% Set up the boundary conditions (using old mesh):
% this is where the most work went
% need to change depending on problem
bm1Name = @(x,y) zeros(size(x)); %zero x displacement on all boundaries
m1_e = zeros(length(mFem.point),1);
m2_e = zeros(length(mFem.point),1);
I = find(nt_m == 0);
m1_e(I) = bm1Name(mFem.point(1,I)', mFem.point(2,I)');
if kt == 1 % it starts off as 1D since it ramps, need to change dep on problem
    m2_e(I) = bm1Name(mFem.point(1,I)', mFem.point(2,I)'); % zeros
else
    m2_e(I) = bm2Name(mFem.point(1,I)', mFem.point(2,I)', beamDisplacement, MeshBeam, wFem);
end

%% Set material (mesh scaling) parameters
lambda = 1; mu = 1;

%% prepare the matrix blocks
% note that Mhuu_hyx = Mhuu_hxy is probably true, but lets be safe
aName = @(x,y) ones(size(x)); nQuadraturePoint = 7;
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0;
Mhuu_hxx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName,MeshFluid,...
    mFem,xDerivative1,yDerivative1,mFem,xDerivative2,yDerivative2,nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Mhuu_hyy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName,MeshFluid,...
    mFem,xDerivative1,yDerivative1,mFem,xDerivative2,yDerivative2,nQuadraturePoint);
xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 1; yDerivative2 = 0;
Mhuu_hyx = stiffness_matrix_assembler_2d_lagrange_tri_global(aName,MeshFluid,...
    mFem,xDerivative1,yDerivative1,mFem,xDerivative2,yDerivative2,nQuadraturePoint);
xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 0; yDerivative2 = 1;
Mhuu_hxy = stiffness_matrix_assembler_2d_lagrange_tri_global(aName,MeshFluid,...
    mFem,xDerivative1,yDerivative1,mFem,xDerivative2,yDerivative2,nQuadraturePoint);

Muu_hxx = Mhuu_hxx(dof_m, dof_m);
Muu_hyy = Mhuu_hyy(dof_m, dof_m);
Muu_hyx = Mhuu_hyx(dof_m, dof_m);
Muu_hxy = Mhuu_hxy(dof_m, dof_m);

% form the matrix for the FE system:
%SP0_u = sparse(length(dof_u), length(dof_u));
M = [mu*(Muu_hxx + Muu_hyy) + (lambda + mu)*Muu_hxx, (lambda + mu)*Muu_hyx; ...
    (lambda + mu)*Muu_hxy, mu*(Muu_hxx + Muu_hyy) + (lambda + mu)*Muu_hyy];

% prepare the vector blocks:
f1Name = @(x,y) zeros(size(x));
f2Name = @(x,y) zeros(size(y));
xDerivative = 0; yDerivative = 0;
f1h_h = load_vector_assembler_2d_lagrange_tri_global(f1Name, MeshFluid, mFem, xDerivative, ...
    yDerivative, nQuadraturePoint);
f1_h = f1h_h(dof_m);
f2h_h = load_vector_assembler_2d_lagrange_tri_global(f2Name, MeshFluid, mFem, xDerivative, ...
    yDerivative, nQuadraturePoint);
f2_h = f2h_h(dof_m);

% prepare the boundary vector blocks
tmp = -mu*Mhuu_hxx*m1_e; bc_u1_eMhxx_mu = tmp(dof_m);
tmp = -(lambda + mu)*Mhuu_hxx*m1_e; bc_u1_eMhxx_lambda_mu = tmp(dof_m);
tmp = -mu*Mhuu_hyy*m1_e; bc_u1_eMhyy_mu = tmp(dof_m);
tmp = -(lambda + mu)*Mhuu_hyx*m2_e; bc_u2_eMhyx_lambda_mu = tmp(dof_m);

tmp = -mu*Mhuu_hxx*m2_e; bc_u2_eMhxx_mu = tmp(dof_m);
tmp = -mu*Mhuu_hyy*m2_e; bc_u2_eMhyy_mu = tmp(dof_m);
tmp = -(lambda + mu)*Mhuu_hyy*m2_e; bc_u2_eMhyy_lambda_mu = tmp(dof_m);
tmp = -(lambda + mu)*Mhuu_hxy*m1_e; bc_u1_eMhxy_lambda_mu = tmp(dof_m);

rhs = [f1_h + bc_u1_eMhxx_mu + bc_u1_eMhxx_lambda_mu + bc_u1_eMhyy_mu + bc_u2_eMhyx_lambda_mu; ...
    f2_h + bc_u2_eMhxx_mu + bc_u2_eMhyy_mu + bc_u2_eMhyy_lambda_mu + bc_u1_eMhxy_lambda_mu];

% solve the FE system
m1m2 = M\rhs;
m1 = m1m2(1:length(dof_m));
m1_fe = m1_e; m1_fe(dof_m) = m1; % m1_fe is displacement in x direction
m2 = m1m2(length(dof_m) + 1:(2*length(dof_m)));
m2_fe = m2_e; m2_fe(dof_m) = m2; % m2_fe is displacement in y direction

%% Note: Fem.point index numbering matches Mesh.node index (for first
% length(Mesh.node) = n_nodes, so:
new_x = MeshFluid.node(1,:)+ m1_fe(1:length(MeshFluid.node))';
new_y =  MeshFluid.node(2,:)+ m2_fe(1:length(MeshFluid.node))'; % new_x and y are row vectors

MeshNew.node = [new_x; new_y];
MeshNew.element = MeshFluid.element; % element connectivity stays the same
MeshNew.edge = MeshFluid.edge; % edge connectivity stays the same
MeshNew.M_e = MeshFluid.M_e; % edge connevitivity stays the same
MeshNew.M_ne = MeshFluid.M_ne;% edge connevitivity stays the same

%% ISOPARAMETRIC ELEMENTS?
end

%% boundary condition sub/local function
function bm2Name = bm2Name(x,y,beamDisplacement,MeshBeam,wFem)

% reso = 20;
% wInterpolate = zeros(1,(reso+1)*size(MeshBeam.element,2));
% xInterpolate = wInterpolate;
% for k = 1:size(MeshBeam.element,2)
%     elementBeam = MeshBeam.node(MeshBeam.element(:,k));
%     wLocal = beamDisplacement(wFem.T(:,k));
%     xl = linspace(elementBeam(1), elementBeam(2), reso+1);
%     xInterpolate(((k-1)*(reso+1)+1):(k*(reso+1))) = xl;
%     wInterpolate(((k-1)*(reso+1)+1):(k*(reso+1))) = evalfeherm1d(xl, wLocal, ...
%         elementBeam, 1, 0); % wInterp = u1*N1 + u2*N2 + u1'*N3 + u2'*N4
% end
%
%
% % match = false(1, numel(xInterpolate));
% % val     = xInterpolate(1);       % Remember the first element
% % index = 1;
% % for k = 2:numel(xInterpolate)  % Loop over elements
% %   if xInterpolate(k) == val      % If current value equals stored value
% %     match(index) = true;  % Set flag, that a(index) is repeated
% %   else              % Value has changed
% %     val     = xInterpolate(k);   % Remeber new value and its index
% %     index = k;
% %   end
% % end
% % result = find(match) + 1;
%
% [CxInterpolate,IC,~] = unique(xInterpolate);
% CwInterpolate = wInterpolate(IC);

%u = unique(xInterpolate);
%h = histc(xInterpolate, u);
%u(h>1); find repeated values
%repeated, indicesfind(xInterpolate == u(h>1));


%beps = 10e-8;

for entry = 1:length(y)
    if ((x(entry) >= 0) && (x(entry) <= 1) && (y(entry) == 1)) % top %y == 1
        bm2Name(entry) = 0;
    elseif ((x(entry) == 0) && (y(entry) >= 0) && (y(entry) <= 1)) % left
        bm2Name(entry) = 0;
    elseif ((x(entry) == 1) && (y(entry) >= 0) && (y(entry) <= 1)) % right
        bm2Name(entry) = 0;
    else % else, on bottom
        for k = 1:size(MeshBeam.element,2)
            elementBeam = MeshBeam.node(MeshBeam.element(:,k));
            if (x(entry) - elementBeam(1))*(x(entry) - elementBeam(2)) <= 0
                wLocal = beamDisplacement(wFem.T(:,k)); 
                bm2Name(entry) = evalfeherm1d(x(entry), wLocal, elementBeam, 0, 0);
            end
        end
        %bm2Name = interp1(CxInterpolate, CwInterpolate, x, 'pchip');
        %bm2Name = interp2(CxInterpolate,CwInterpolate,CwInterpolate, x,y);
    end
    
end


% g2_D = @(x,y) 0.*((x >= 0) & (x <= 1) & (y == 1)) + ...
%     0.*((x == 0) & (y >= 0) & (y <= 1)) + ...
%     0.*((x == 1) & (y >= 0) & (y <= 1)) + ...
%     interp1(xInterpolate, wInterpolate, x, 'pchip').*((x >= 0) & (x <= 1) & (y <= );
end


