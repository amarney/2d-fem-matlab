% issues:
%   problem 19, error for derivative doesn't decrease
%   problem 17, ufe_tx3 and ufe_tx4 are far off (also has derivatives...)



% Problem 15:
% -------------------------------------------------------------------------
% Weak Form:
%   Find u \in S such that
%   
%   \int_{0}^{1} (2+x) u'(x) v'(x) dx + \int_{0}^{1} (1 + x^2) u(x) v(x) dx
%       = \int_{0}^{1} (-3(1+x) cos(x) - (-4 + 2x + x^3)sin(1-x)) v(x) dx
%   
%   for all v \in T where
%   T = {w | w \in H1(0,1), w(0) = 0}
%   S = {w | w \in H1(0,1), w(0) = sin(1)}
% -------------------------------------------------------------------------

% using linear finite element functions:
% step 1: prepare the BVP domain, boundary index, and function handles
domain = [0,1]; iBoundaryCondition = [0,-1]; %left is essential
aName = @(x) 2 + x; 
cName = @(x) 1 + x.^2;
fName = @(x) -3.*(1+x).*cos(1-x) - (-4 + 2.*x + x.^3).*sin(1-x);
eName = @(x) sin(1);
natName = @(x) 0;

% step 2: prepare the finite element trial space, solution set, and 
%         BC data structures 
nElement = 20; Mesh = mesh_generator_1d(domain, nElement);
iDegree = 1; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
[uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition);

% step 3: set the number of gauss points and derivative index, so that we  
%         can assemble global load vector
nGaussPoint = 3; iDerivative = 0;
f = load_vector_assembler_1d_lagrange_global(fName, ... 
    Mesh, Fem, iDerivative, nGaussPoint);


% step 4: similiarly assembler the global stiffness matrices
iDerivative = 1;
S = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
iDerivative = 0;
M = stiffness_matrix_assembler_1d_lagrange_global(cName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

% step 5: extract the vector/matrices which are not associated with the
%         essential BC and form vectors to take care of the essential BC
fh = f(uDof);
uEssential = zeros(length(Fem.point),1); ind = find(nodeType == 0); % find essential
uEssential(ind) = eName(Fem.point(ind));
Sh = S(uDof, uDof); temp = -S*uEssential; bc_eS = temp(uDof);
Mh = M(uDof, uDof); temp = -M*uEssential; bc_eM = temp(uDof);

% step 6: form the vector to take care of the natural BC
uNatural = zeros(length(Fem.point),1); ind = find(nodeType == -1); % find natural
uNatural(ind) = natName(Fem.point(ind));
bc_n = uNatural(uDof); bc_n(1) = -bc_n(1);

% step 7: solve (Sh + Mh)u = (fh + bc_n + bc_eS + bc_eM)
%         for the coefficients u
u = (Sh + Mh)\(fh + bc_n + bc_eS + bc_eM);

% step 8: enforce the essential boundary condition to get the final 
%         coefficients of the L2 projection onto our finite element space
u_fe = uEssential; u_fe(uDof) = u;

% ------------
%   notes:
%   > L2 projection can be evaluated locally using
%     evaluate_fe_function_1d_lagrange
%   > it would be nice to initial iDegree = [0, 1, 2, 3]
%     so that iDegree(0) = 0, but Matlab won't do that.
% ------------

% step 9: postprocessing (evaluating u_fe(pi/5)
tx = [0,1/5,1/4,1/3,1/2,pi/5,1];
for entry = 1:length(tx)
    for k = 1:size(Mesh.element,2)
        element = Mesh.node(Mesh.element(:,k));
        if (tx(entry) - element(1))*(tx(entry) - element(2)) <= 0
            u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
            ufe_tx1(entry) = evaluate_fe_function_1d_lagrange(tx(entry), u_fe_loc, element, iDegree, 0);
        end
    end
end

% using quadratic finite element functions:
% step 1: prepare the BVP domain, boundary index, and function handles
domain = [0,1]; iBoundaryCondition = [0,-1]; %left is essential
aName = @(x) 2 + x; 
cName = @(x) 1 + x.^2;
fName = @(x) -3.*(1+x).*cos(1-x) - (-4 + 2.*x + x.^3).*sin(1-x);
eName = @(x) sin(1);
natName = @(x) 0;

% step 2: prepare the finite element trial space, solution set, and 
%         BC data structures 
nElement = 20; Mesh = mesh_generator_1d(domain, nElement);
iDegree = 2; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
[uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition);

% step 3: set the number of gauss points and derivative index, so that we  
%         can assemble global load vector
nGaussPoint = 3; iDerivative = 0;
f = load_vector_assembler_1d_lagrange_global(fName, ... 
    Mesh, Fem, iDerivative, nGaussPoint);

% step 4: similiarly assembler the global stiffness matrices
iDerivative = 1;
S = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
iDerivative = 0;
M = stiffness_matrix_assembler_1d_lagrange_global(cName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

% step 5: extract the vector/matrices which are not associated with the
%         essential BC and form vectors to take care of the essential BC
fh = f(uDof);
uEssential = zeros(length(Fem.point),1); ind = find(nodeType == 0); % find essential
uEssential(ind) = eName(Fem.point(ind));
Sh = S(uDof, uDof); temp = -S*uEssential; bc_eS = temp(uDof);
Mh = M(uDof, uDof); temp = -M*uEssential; bc_eM = temp(uDof);

% step 6: form the vector to take care of the natural BC
uNatural = zeros(length(Fem.point),1); ind = find(nodeType == -1); % find natural
uNatural(ind) = natName(Fem.point(ind));
bc_n = uNatural(uDof); bc_n(1) = -bc_n(1);

% step 7: solve (Sh + Mh)u = (fh + bc_n + bc_eS + bc_eM)
%         for the coefficients u
u = (Sh + Mh)\(fh + bc_n + bc_eS + bc_eM);

% step 8: enforce the essential boundary condition to get the final 
%         coefficients of the L2 projection onto our finite element space
u_fe = uEssential; u_fe(uDof) = u;

% step 9: postprocessing (evaluating u_fe(pi/5))
tx = pi/5;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    if (tx - element(1))*(tx - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        ufe_tx2 = evaluate_fe_function_1d_lagrange(tx, u_fe_loc, element, iDegree, 0);
    end
end

% Problem 16:
% -------------------------------------------------------------------------
% Weak Form:
%   Find u \in S such that
%   
%   \int_{0}^{1} (2+cos(x)) u'(x) v'(x) dx + \int_{0}^{1} (1 + x^2) u(x) v(x) dx
%       = \int_{0}^{1} (-e^(x^2) * (3 + 7x^2 + (2+4x^2)cos(x) - 2x sin(x))) v(x) dx
%   
%   for all v \in T where
%   T = {w | w \in H1(0,1), w(1) = 0}
%   S = {w | w \in H1(0,1), w(1) = e^1}
% -------------------------------------------------------------------------

% using linear finite element functions:
% step 1: prepare the BVP domain, boundary index, and function handles
domain = [0,1]; iBoundaryCondition = [-1,0]; % right is essential
aName = @(x) 2 + cos(x); 
cName = @(x) 1 + x.^2;
fName = @(x) -exp(x.^2).*(3 + 7*x.^2 + (2+4*x.^2).*cos(x) - 2*x.*sin(x));
eName = @(x) exp(x);
natName = @(x) 0;

% step 2: prepare the finite element trial space, solution set, and 
%         BC data structures 
nElement = 20; Mesh = mesh_generator_1d(domain, nElement);
iDegree = 1; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
[uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition);

% step 3: set the number of gauss points and derivative index, so that we  
%         can assemble global load vector
nGaussPoint = 3; iDerivative = 0;
f = load_vector_assembler_1d_lagrange_global(fName, ... 
    Mesh, Fem, iDerivative, nGaussPoint);

% step 4: similiarly assembler the global stiffness matrices
iDerivative = 1;
S = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
iDerivative = 0;
M = stiffness_matrix_assembler_1d_lagrange_global(cName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

% step 5: extract the vector/matrices which are not associated with the
%         essential BC and form vectors to take care of the essential BC
fh = f(uDof);
uEssential = zeros(length(Fem.point),1); ind = find(nodeType == 0); % find essential
uEssential(ind) = eName(Fem.point(ind));
Sh = S(uDof, uDof); temp = -S*uEssential; bc_eS = temp(uDof);
Mh = M(uDof, uDof); temp = -M*uEssential; bc_eM = temp(uDof);

% step 6: form the vector to take care of the natural BC
uNatural = zeros(length(Fem.point),1); ind = find(nodeType == -1); % find natural
uNatural(ind) = natName(Fem.point(ind));
bc_n = uNatural(uDof); bc_n(1) = -bc_n(1);

% step 7: solve (Sh + Mh)u = (fh + bc_n + bc_eS + bc_eM)
%         for the coefficients u
u = (Sh + Mh)\(fh + bc_n + bc_eS + bc_eM);

% step 8: enforce the essential boundary condition to get the final 
%         coefficients of the L2 projection onto our finite element space
u_fe = uEssential; u_fe(uDof) = u;

% step 9: postprocessing (integrating u_fe(x))
gaussPoint = gauss_node_generator_1d_local([0,1],3);
gaussWeight = gauss_weight_generator_1d_local([0,1],3);
u_fe_gaussPoint = zeros(3,1);
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    if (gaussPoint(1) - element(1))*(gaussPoint(1) - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        u_fe_gaussPoint(1) = evaluate_fe_function_1d_lagrange(gaussPoint(1), u_fe_loc, element, iDegree, 0);
    elseif (gaussPoint(2) - element(1))*(gaussPoint(2) - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        u_fe_gaussPoint(2) = evaluate_fe_function_1d_lagrange(gaussPoint(2), u_fe_loc, element, iDegree, 0);
    elseif (gaussPoint(3) - element(1))*(gaussPoint(3) - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        u_fe_gaussPoint(3) = evaluate_fe_function_1d_lagrange(gaussPoint(3), u_fe_loc, element, iDegree, 0);
    end
end
integral1 = gaussWeight*u_fe_gaussPoint;

% using quadratic finite element functions:
% step 1: prepare the BVP domain, boundary index, and function handles
domain = [0,1]; iBoundaryCondition = [-1,0]; % right is essential
aName = @(x) 2 + cos(x); 
cName = @(x) 1 + x.^2;
fName = @(x) -exp(x.^2).*(3 + 7*x.^2 + (2+4*x.^2).*cos(x) - 2*x.*sin(x));
eName = @(x) exp(x);
natName = @(x) 0;

% step 2: prepare the finite element trial space, solution set, and 
%         BC data structures 
nElement = 20; Mesh = mesh_generator_1d(domain, nElement);
iDegree = 2; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
[uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition);

% step 3: set the number of gauss points and derivative index, so that we  
%         can assemble global load vector
nGaussPoint = 3; iDerivative = 0;
f = load_vector_assembler_1d_lagrange_global(fName, ... 
    Mesh, Fem, iDerivative, nGaussPoint);

% step 4: similiarly assembler the global stiffness matrices
iDerivative = 1;
S = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
iDerivative = 0;
M = stiffness_matrix_assembler_1d_lagrange_global(cName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

% step 5: extract the vector/matrices which are not associated with the
%         essential BC and form vectors to take care of the essential BC
fh = f(uDof);
uEssential = zeros(length(Fem.point),1); ind = find(nodeType == 0); % find essential
uEssential(ind) = eName(Fem.point(ind));
Sh = S(uDof, uDof); temp = -S*uEssential; bc_eS = temp(uDof);
Mh = M(uDof, uDof); temp = -M*uEssential; bc_eM = temp(uDof);

% step 6: form the vector to take care of the natural BC
uNatural = zeros(length(Fem.point),1); ind = find(nodeType == -1); % find natural
uNatural(ind) = natName(Fem.point(ind));
bc_n = uNatural(uDof); bc_n(1) = -bc_n(1);

% step 7: solve (Sh + Mh)u = (fh + bc_n + bc_eS + bc_eM)
%         for the coefficients u
u = (Sh + Mh)\(fh + bc_n + bc_eS + bc_eM);

% step 8: enforce the essential boundary condition to get the final 
%         coefficients of the L2 projection onto our finite element space
u_fe = uEssential; u_fe(uDof) = u;

% step 9: postprocessing (evaluating u_fe(pi/5))
gaussPoint = gauss_node_generator_1d_local([0,1],3);
gaussWeight = gauss_weight_generator_1d_local([0,1],3);
u_fe_gaussPoint = zeros(3,1);
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    if (gaussPoint(1) - element(1))*(gaussPoint(1) - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        u_fe_gaussPoint(1) = evaluate_fe_function_1d_lagrange(gaussPoint(1), u_fe_loc, element, iDegree, 0);
    elseif (gaussPoint(2) - element(1))*(gaussPoint(2) - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        u_fe_gaussPoint(2) = evaluate_fe_function_1d_lagrange(gaussPoint(2), u_fe_loc, element, iDegree, 0);
    elseif (gaussPoint(3) - element(1))*(gaussPoint(3) - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        u_fe_gaussPoint(3) = evaluate_fe_function_1d_lagrange(gaussPoint(3), u_fe_loc, element, iDegree, 0);
    end
end
integral2 = gaussWeight*u_fe_gaussPoint;

% Problem 17:
% -------------------------------------------------------------------------
% Weak Form:
%   Find u \in S such that
%   
%   \int_{0}^{1} (2+cos(x)) u'(x) v'(x) dx + \int_{0}^{1} (1 + x^2) u(x) v(x) dx
%       = \int_{0}^{1} e^x * (-6 -8*x + x^2 + x^4 - (2+x)^2 cos(x) + (2+2x+x^2) sin(x)) v(x) dx
%   
%   for all v \in T where
%   T = {w | w \in H1(0,1), w(0) = 0, w(1) = 0}
%   S = {w | w \in H1(0,1), w(0) = 2, w(1) = 3*e}
% -------------------------------------------------------------------------

% using linear finite element functions:
% step 1: prepare the BVP domain, boundary index, and function handles
domain = [0,1]; iBoundaryCondition = [0,0]; % both are essential
aName = @(x) 2 + cos(x); 
cName = @(x) 1 + x.^2;
fName = @(x) exp(x).*(-6 - 8*x + x.^2 + x.^4 - (2+x).^2 .*cos(x) + (2+2*x+x.^2).*sin(x));
eName = @(x) ...
    2*(x < 1) ...
    + 3*exp(x).*(x >= 1);
natName = @(x) 0;

% step 2: prepare the finite element trial space, solution set, and 
%         BC data structures 
nElement = 20; Mesh = mesh_generator_1d(domain, nElement);
iDegree = 1; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
[uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition);

% step 3: set the number of gauss points and derivative index, so that we  
%         can assemble global load vector
nGaussPoint = 3; iDerivative = 0;
f = load_vector_assembler_1d_lagrange_global(fName, ... 
    Mesh, Fem, iDerivative, nGaussPoint);

% step 4: similiarly assembler the global stiffness matrices
iDerivative = 1;
S = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
iDerivative = 0;
M = stiffness_matrix_assembler_1d_lagrange_global(cName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

% step 5: extract the vector/matrices which are not associated with the
%         essential BC and form vectors to take care of the essential BC
fh = f(uDof);
uEssential = zeros(length(Fem.point),1); ind = find(nodeType == 0); % find essential
uEssential(ind) = eName(Fem.point(ind));
Sh = S(uDof, uDof); temp = -S*uEssential; bc_eS = temp(uDof);
Mh = M(uDof, uDof); temp = -M*uEssential; bc_eM = temp(uDof);

% step 6: form the vector to take care of the natural BC
uNatural = zeros(length(Fem.point),1); ind = find(nodeType == -1); % find natural
uNatural(ind) = natName(Fem.point(ind));
bc_n = uNatural(uDof); bc_n(1) = -bc_n(1);

% step 7: solve (Sh + Mh)u = (fh + bc_n + bc_eS + bc_eM)
%         for the coefficients u
u = (Sh + Mh)\(fh + bc_n + bc_eS + bc_eM);

% step 8: enforce the essential boundary condition to get the final 
%         coefficients of the L2 projection onto our finite element space
u_fe = uEssential; u_fe(uDof) = u;

% step 9: postprocessing (evaluating u_fe(pi/4)')
tx = pi/4;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    if (tx - element(1))*(tx - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        ufe_tx3 = evaluate_fe_function_1d_lagrange(tx, u_fe_loc, element, iDegree, 1);
    end
end

% using quadratic finite element functions:
% step 1: prepare the BVP domain, boundary index, and function handles
domain = [0,1]; iBoundaryCondition = [0,0]; % both are essentials
aName = @(x) 2 + cos(x); 
cName = @(x) 1 + x.^2;
fName = @(x) exp(x).*(-6 - 8*x + x.^2 + x.^4 - (2+x).^2 .*cos(x) + (2+2*x+x.^2).*sin(x));
eName = @(x) ...
    2*(x < 1) ...
    + 3*exp(x).*(x >= 1);
natName = @(x) 0;

% step 2: prepare the finite element trial space, solution set, and 
%         BC data structures 
nElement = 20; Mesh = mesh_generator_1d(domain, nElement);
iDegree = 2; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
[uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition);

% step 3: set the number of gauss points and derivative index, so that we  
%         can assemble global load vector
nGaussPoint = 3; iDerivative = 0;
f = load_vector_assembler_1d_lagrange_global(fName, ... 
    Mesh, Fem, iDerivative, nGaussPoint);

% step 4: similiarly assembler the global stiffness matrices
iDerivative = 1;
S = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
iDerivative = 0;
M = stiffness_matrix_assembler_1d_lagrange_global(cName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

% step 5: extract the vector/matrices which are not associated with the
%         essential BC and form vectors to take care of the essential BC
fh = f(uDof);
uEssential = zeros(length(Fem.point),1); ind = find(nodeType == 0); % find essential
uEssential(ind) = eName(Fem.point(ind));
Sh = S(uDof, uDof); temp = -S*uEssential; bc_eS = temp(uDof);
Mh = M(uDof, uDof); temp = -M*uEssential; bc_eM = temp(uDof);

% step 6: form the vector to take care of the natural BC
uNatural = zeros(length(Fem.point),1); ind = find(nodeType == -1); % find natural
uNatural(ind) = natName(Fem.point(ind));
bc_n = uNatural(uDof); bc_n(1) = -bc_n(1);

% step 7: solve (Sh + Mh)u = (fh + bc_n + bc_eS + bc_eM)
%         for the coefficients u
u = (Sh + Mh)\(fh + bc_n + bc_eS + bc_eM);

% step 8: enforce the essential boundary condition to get the final 
%         coefficients of the L2 projection onto our finite element space
u_fe = uEssential; u_fe(uDof) = u;

% step 9: postprocessing (evaluating u_fe(pi/4)')
tx = pi/4;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    if (tx - element(1))*(tx - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        ufe_tx4 = evaluate_fe_function_1d_lagrange(tx, u_fe_loc, element, iDegree, 1);
    end
end

% Problem 18:
% -------------------------------------------------------------------------
% Weak Form:
%   Find u \in S such that
%   
%   \int_{0}^{1} (2+x) u'(x) v'(x) dx + \int_{0}^{1} (1 + x) u'(x) v(x) dx + \int_{0}^{1} (1 + x^2) u(x) v(x) dx
%       = \int_{0}^{1} e^x * (-6 -8*x + x^2 + x^4 - (2+x)^2 cos(x) + (2+2x+x^2) sin(x)) v(x) dx
%   
%   for all v \in T where
%   T = {w | w \in H1(0,1), w(0) = 0, w(1) = 0}
%   S = {w | w \in H1(0,1), w(0) = 1, w(1) = 4*cos(1)}
% -------------------------------------------------------------------------

% using linear finite element functions:
% step 1: prepare the BVP domain, boundary index, and function handles
domain = [0,1]; iBoundaryCondition = [0,0]; % both are essential
aName = @(x) 2 + x;
bName = @(x) 1 + x;
cName = @(x) 1 + x.^2;
fName = @(x) (3 + 13*x + 4*x.^2 + 3*x.^3).*cos(x) + (12 + 5*x - 3*x.^2).*sin(x);
eName = @(x) ...
    1*(x < 1) ...
    + 4*cos(x).*(x >= 1);
natName = @(x) 0;

% step 2: prepare the finite element trial space, solution set, and 
%         BC data structures 
nElement = 20; Mesh = mesh_generator_1d(domain, nElement);
iDegree = 1; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
[uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition);

% step 3: set the number of gauss points and derivative index, so that we  
%         can assemble global load vector
nGaussPoint = 3; iDerivative = 0;
f = load_vector_assembler_1d_lagrange_global(fName, ... 
    Mesh, Fem, iDerivative, nGaussPoint);

% step 4: similiarly assembler the global stiffness matrices
iDerivative = 1;
S = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
iDerivative = [0,1];
N = stiffness_matrix_assembler_1d_lagrange_global(bName, ...
    Mesh, Fem, iDerivative(1), Fem, iDerivative(2), nGaussPoint);
iDerivative = 0;
M = stiffness_matrix_assembler_1d_lagrange_global(cName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

% step 5: extract the vector/matrices which are not associated with the
%         essential BC and form vectors to take care of the essential BC
fh = f(uDof);
uEssential = zeros(length(Fem.point),1); ind = find(nodeType == 0); % find essential
uEssential(ind) = eName(Fem.point(ind));
Sh = S(uDof, uDof); temp = -S*uEssential; bc_eS = temp(uDof);
Nh = N(uDof, uDof); temp = -N*uEssential; bc_eN = temp(uDof);
Mh = M(uDof, uDof); temp = -M*uEssential; bc_eM = temp(uDof);

% step 6: form the vector to take care of the natural BC
uNatural = zeros(length(Fem.point),1); ind = find(nodeType == -1); % find natural
uNatural(ind) = natName(Fem.point(ind));
bc_n = uNatural(uDof); bc_n(1) = -bc_n(1);

% step 7: solve Au = b  for the coefficients u
u = (Sh + Nh + Mh)\(fh + bc_n + bc_eS + bc_eN + bc_eM);

% step 8: enforce the essential boundary condition to get the final 
%         coefficients of the L2 projection onto our finite element space
u_fe = uEssential; u_fe(uDof) = u;

% step 9: postprocessing (evaluating u_fe(pi/7))
tx = pi/7;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    if (tx - element(1))*(tx - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        ufe_tx5 = evaluate_fe_function_1d_lagrange(tx, u_fe_loc, element, iDegree, 0);
    end
end

% using quadratic finite element functions:
% step 1: prepare the BVP domain, boundary index, and function handles
domain = [0,1]; iBoundaryCondition = [0,0]; % both are essential
aName = @(x) 2 + x;
bName = @(x) 1 + x;
cName = @(x) 1 + x.^2;
fName = @(x) (3 + 13*x + 4*x.^2 + 3*x.^3).*cos(x) + (12 + 5*x - 3*x.^2).*sin(x);
eName = @(x) ...
    1*(x < 1) ...
    + 4*cos(x).*(x >= 1);
natName = @(x) 0;

% step 2: prepare the finite element trial space, solution set, and 
%         BC data structures 
nElement = 20; Mesh = mesh_generator_1d(domain, nElement);
iDegree = 2; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
[uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition);

% step 3: set the number of gauss points and derivative index, so that we  
%         can assemble global load vector
nGaussPoint = 3; iDerivative = 0;
f = load_vector_assembler_1d_lagrange_global(fName, ... 
    Mesh, Fem, iDerivative, nGaussPoint);

% step 4: similiarly assembler the global stiffness matrices
iDerivative = 1;
S = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
iDerivative = [0,1];
N = stiffness_matrix_assembler_1d_lagrange_global(bName, ...
    Mesh, Fem, iDerivative(1), Fem, iDerivative(2), nGaussPoint);
iDerivative = 0;
M = stiffness_matrix_assembler_1d_lagrange_global(cName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

% step 5: extract the vector/matrices which are not associated with the
%         essential BC and form vectors to take care of the essential BC
fh = f(uDof);
uEssential = zeros(length(Fem.point),1); ind = find(nodeType == 0); % find essential
uEssential(ind) = eName(Fem.point(ind));
Sh = S(uDof, uDof); temp = -S*uEssential; bc_eS = temp(uDof);
Nh = N(uDof, uDof); temp = -N*uEssential; bc_eN = temp(uDof);
Mh = M(uDof, uDof); temp = -M*uEssential; bc_eM = temp(uDof);

% step 6: form the vector to take care of the natural BC
uNatural = zeros(length(Fem.point),1); ind = find(nodeType == -1); % find natural
uNatural(ind) = natName(Fem.point(ind));
bc_n = uNatural(uDof); bc_n(1) = -bc_n(1);

% step 7: solve Au = b  for the coefficients u
u = (Sh + Nh + Mh)\(fh + bc_n + bc_eS + bc_eN + bc_eM);

% step 8: enforce the essential boundary condition to get the final 
%         coefficients of the L2 projection onto our finite element space
u_fe = uEssential; u_fe(uDof) = u;

% step 9: postprocessing (evaluating u_fe(pi/7))
tx = pi/7;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    if (tx - element(1))*(tx - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        ufe_tx6 = evaluate_fe_function_1d_lagrange(tx, u_fe_loc, element, iDegree, 0);
    end
end

% Problem 19:
% -------------------------------------------------------------------------
% Weak Form:
%   Find u \in S such that
%   
%   \int_{0}^{1} (2+e^x) u'(x) v'(x) dx + \int_{0}^{1} (1 + x^2) u(x) v(x) dx
%       = \int_{0}^{1} [(e^x (x-1) + x(3 + x^2)) *cos(x) + (4 + e^x (2+x))*sin(x)] v(x) dx
%   
%   for all v \in T where
%   T = {w | w \in H1(0,1), w(0) = 0, w(1) = 0}
%   S = {w | w \in H1(0,1), w(0) = 0, w(1) = 0}
% -------------------------------------------------------------------------

% using cubic finite element functions:
% step 1: prepare the BVP domain, boundary index, and function handles
domain = [0,1]; iBoundaryCondition = [-1,-1]; % none are essential
aName = @(x) 2 + exp(x); 
cName = @(x) 1 + x.^2; 
fName = @(x) (exp(x).*(x-1) + x.*(x.^2 + 3)).*cos(x) + (4 + exp(x).*(x + 2)).*sin(x); 
eName = @(x) 0; 
natName = @(x) ...
    3.*(x < 1) ...
    + (exp(x)+2).*(cos(x) - sin(x)).*(x >= 1);
uName = @(x) x.*cos(x);
uPrimeName = @(x) cos(x) - x.*sin(x);

% Run all the steps in a loop
errorL2 = zeros(1,7);
errorH1 = zeros(1,7);
NN = [10, 20, 30, 40, 50, 60, 70]; % for mesh size
for count = 1:7
    nElement = 10*count;
    Mesh = mesh_generator_1d(domain, nElement);
    iDegree = 3; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
    [uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition);
    nGaussPoint = 4; iDerivative = 0;
    f = load_vector_assembler_1d_lagrange_global(fName, ... 
        Mesh, Fem, iDerivative, nGaussPoint);
    iDerivative = 1;
    S = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
        Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
    iDerivative = 0;
    M = stiffness_matrix_assembler_1d_lagrange_global(cName, ...
        Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);
    fh = f(uDof);
    uEssential = zeros(length(Fem.point),1); ind = find(nodeType == 0); % find essential
    uEssential(ind) = eName(Fem.point(ind));
    Sh = S(uDof, uDof); temp = -S*uEssential; bc_eS = temp(uDof);
    Mh = M(uDof, uDof); temp = -M*uEssential; bc_eM = temp(uDof);
    uNatural = zeros(length(Fem.point),1); ind = find(nodeType == -1); % find natural
    uNatural(ind) = natName(Fem.point(ind));
    bc_n = uNatural(uDof); bc_n(1) = -bc_n(1);
    u = (Sh + Mh)\(fh + bc_n + bc_eS + bc_eM);
    u_fe = uEssential; u_fe(uDof) = u;
    
    iDerivative = 0; nGaussPoint = 4;
    errorL2(count) = sqrt(error_1d_lagrange_global(uName, u_fe, Mesh, Fem, iDerivative, nGaussPoint));
    iDerivative = 1; nGaussPoint = 4;
    errorH1(count) = sqrt(error_1d_lagrange_global(uPrimeName, u_fe, Mesh, Fem, iDerivative, nGaussPoint));
end
h = 1./NN;
rc = polyfit(log(h), log(errorL2), 1);
r_L2 = rc(1);
C_L2 = rc(2);
rc = polyfit(log(h), log(errorH1), 1);
r_H1 = rc(1);
C_H1 = rc(2);

clear

% Problem 19:
% -------------------------------------------------------------------------
% Weak Form:
%   Find u \in S such that
%   
%   \int_{0}^{1} (2+cos(x)) u'(x) v'(x) dx 
%       = \int_{0}^{1} e^x(sin(x)(1+x) - (2+x)(2+cos(x))) v(x) dx
%   
%   for all v \in T where
%   T = {w | w \in H1(0,1), w(0) = 0, w(1) = 0}
%   S = T
% -------------------------------------------------------------------------

% using linear finite element functions:
% step 1: prepare the BVP domain, boundary index, and function handles
domain = [0,1]; iBoundaryCondition = [-1,-1]; % both are essential
aName = @(x) 2 + cos(x);
fName = @(x) exp(x).*(sin(x).*(1+x) - (2+x).*(2 + cos(x)));
eName = @(x) 0;
natName = @(x) ...
    3*(x < 1) ...
    + 2*exp(x).*(2+cos(x)).*(x >= 1);


% step 2: prepare the finite element trial space, solution set, and 
%         BC data structures 
nElement = 20; Mesh = mesh_generator_1d(domain, nElement);
iDegree = 1; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
[uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition);

% step 3: set the number of gauss points and derivative index, so that we  
%         can assemble global load vector
nGaussPoint = 3; iDerivative = 0;
f = load_vector_assembler_1d_lagrange_global(fName, ... 
    Mesh, Fem, iDerivative, nGaussPoint);

% step 4: similiarly assembler the global stiffness matrices
iDerivative = 1;
S = stiffness_matrix_assembler_1d_lagrange_global(aName, ...
    Mesh, Fem, iDerivative, Fem, iDerivative, nGaussPoint);

% step 5: extract the vector/matrices which are not associated with the
%         essential BC and form vectors to take care of the essential BC
fh = f(uDof);
uEssential = zeros(length(Fem.point),1); ind = find(nodeType == 0); % find essential
uEssential(ind) = eName(Fem.point(ind));
Sh = S(uDof, uDof); temp = -S*uEssential; bc_eS = temp(uDof);

% step 6: form the vector to take care of the natural BC
uNatural = zeros(length(Fem.point),1); ind = find(nodeType == -1); % find natural
uNatural(ind) = natName(Fem.point(ind));
bc_n = uNatural(uDof); bc_n(1) = -bc_n(1);

% step 7: solve Au = b  for the coefficients u
u = (Sh)\(fh + bc_n + bc_eS);

% step 8: enforce the essential boundary condition to get the final 
%         coefficients of the L2 projection onto our finite element space
u_fe = uEssential; u_fe(uDof) = u;

% step 9: postprocessing (evaluating u_fe(pi/6))
tx = pi/6;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    if (tx - element(1))*(tx - element(2)) <= 0
        u_fe_loc = u_fe(Fem.T(:,k)); % we have found the correct element
        ufe_final = evaluate_fe_function_1d_lagrange(tx, u_fe_loc, element, iDegree, 0);
    end
end