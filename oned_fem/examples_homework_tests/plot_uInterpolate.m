% interpolate sin(4*pi*x)
% used to test evalfeherm1d

% Make the mesh:
domain = [0,1]; nElement = 20;
Mesh = mesh_generator_1d(domain, nElement);

% Form a data structure to access the finite element basis:
iDegree = 3; Fem = femherm1d(Mesh, iDegree); hDegree = 1;

% Create an array with the first 100 entries corresponding to the values
% of sin(4*pi*x) evaluated at the Mesh nodes:
uGlobal = sin(4*pi*Fem.point(1:nElement+1)); iDerivative = 0;

% Extend the array so that the last 100 entries correspond to the values of
% 2*pi*cos(2*pi*x) evaluated at the Mesh nodes:
uGlobal(nElement+2:2*(nElement+1)) = 4*pi*cos(4*pi*Fem.point(nElement+2:2*(nElement+1)));

% Set the resolution and prepare the figure:
reso = 20; figure(1); clf; hold on;
uInterpolate = zeros(1,(reso+1)*size(Mesh.element,2));

r = linspace(-1,1,11);
for k = 1:size(Mesh.element,2) % Loop over each element in the mesh
    element = Mesh.node(Mesh.element(:,k)); % Extract element endpoints
    uLocal = uGlobal(Fem.T(:,k)); % uLocal(1) = u1, uLocal(2) = u2, uLocal(3) = u1', uLocal(4) = u2'
    len = element(2) - element(1);
    c0 = len/2;
    c1 = (element(2) + element(1))/2; % rule = 11
    x_g = c0*r + c1;
    %x = linspace(element(1), element(2), reso+1); % x contains the points we are evaluating the shape functions at
    %uInterpolate(((k-1)*(reso+1)+1):(k*(reso+1))) = evalfeherm1d(x, uLocal, ...
    %    element, hDegree, 0); % Evaluate shape functions at x and scale by corresponding coefficients in uLocal 
    plot(x_g,evalfeherm1d(x_g, uLocal, ...
        element, hDegree, 0)); % Plot the result
end
hold off


% domain = [-1,1]; nElement = 1;
% Mesh = mesh_generator_1d(domain, nElement);
% iDegree = 3; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
% uGlobal = sin(2*pi*Fem.point); iDerivative = 0;
% reso = 20; figure(2); hold on;
% uInterpolate = zeros(1,(reso+1)*size(Mesh.element,2));
% for k = 1:size(Mesh.element,2) %loop over element
%     element = Mesh.node(Mesh.element(:,k));
%     uLocal = uGlobal(Fem.T(:,k));
%     x = linspace(element(1), element(2), reso+1);
%     uInterpolate(((k-1)*(reso+1)+1):(k*(reso+1))) = evaluate_fe_function_1d_lagrange(x, uLocal, ...
%         element, iDegree, iDerivative);
%     plot(x,uInterpolate(((k-1)*(reso+1)+1):(k*(reso+1))));
%     pause
% end