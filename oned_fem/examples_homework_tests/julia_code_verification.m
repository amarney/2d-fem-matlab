% make a data structure for the mesh
domain = [0,1];
nElement = 20;
Mesh = mesh_generator_1d(domain, nElement);

% make a data structure to access the fe shape
iDegree = 2;
Fem = fem_generator_1d_lagrange(Mesh, iDegree);

% use galerkin
uFem = Fem;
vFem = Fem;

% make bc data structs
iBoundaryCondition = [-1, 0];
[uDof, nodeType] = bc_array_generator_1d(Fem, domain, iBoundaryCondition);

% Assemble f_0
fName = @(x) exp(x).*(sin(x) - cos(x) - 1 + 1./(x + 1));
%fName = @(x) exp(x).*(sin(x) - cos(x) + x);
iDerivative = 0;
nGaussPoint = 5;
fh_0 = load_vector_assembler_1d_lagrange_global(fName, ... 
    Mesh, Fem, iDerivative, nGaussPoint);

f_0 = fh_0(uDof);

% Assemble A_11
aName = @(x) 1 + cos(x);
iDerivative = 1;
Ah = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

A_11 = Ah(uDof, uDof);

% Assemble C_00
cName = @(x) 1./(x + 1);
%cName = @(x) x + ones(size(x));
iDerivative = 0;
Ch = stiffness_matrix_assembler_1d_lagrange_global(cName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

C_00 = Ch(uDof, uDof);

% Assemble b_eA and b_eC
essName = @(x) exp(1);
uEss = zeros(length(Fem.point),1);
indEss = find(nodeType == 0);
uEss(indEss) = essName(Fem.point(indEss));
temp = -Ah*uEss; bc_eA = temp(uDof);
temp = -Ch*uEss; bc_eC = temp(uDof);

% assemble b_n
natName = @(x) -2*ones(size(x));
uNat = zeros(length(Fem.point),1);
indNat = find(nodeType == -1);
uNat(indNat) = natName(Fem.point(indNat));
b_n = uNat(uDof);
b_n(1) = b_n(1);

% solve system
uhGlobal = (A_11 + C_00)\(f_0 + bc_eA + bc_eC + b_n);
uGlobal = uEss;
uGlobal(uDof) = uhGlobal;

tx = 0.5;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    if (tx - element(1))*(tx - element(2)) <= 0
        uLocal = uGlobal(Fem.T(:,k)); % we have found the correct element
        utx = evaluate_fe_function_1d_lagrange(tx, uLocal, element, iDegree, 0);
    end
end






