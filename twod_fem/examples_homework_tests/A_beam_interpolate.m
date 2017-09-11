% Make the mesh:
domain = [2,3]; nElement = 5;
Mesh = mesh_generator_1d(domain, nElement);

% Form a data structure to access the finite element basis:
iDegree = 3; Fem = femherm1d(Mesh, iDegree); hDegree = 1;

% Create an array with the first 100 entries corresponding to the values
% of sin(2*pi*x) evaluated at the Mesh nodes:
uGlobal = sin(2*pi*Fem.point(1:nElement+1)); iDerivative = 0;

% Extend the array so that the last 100 entries correspond to the values of
% 2*pi*cos(2*pi*x) evaluated at the Mesh nodes:
uGlobal(nElement+2:2*(nElement+1)) = 2*pi*cos(2*pi*Fem.point(nElement+2:2*(nElement+1)));

% Set the resolution and prepare the figure:
reso = 20; figure(1); clf; hold on;
uInterpolate = zeros(1,(reso+1)*size(Mesh.element,2));
xInterpolate = uInterpolate;

for k = 1:size(Mesh.element,2) % Loop over each element in the mesh
    element = Mesh.node(Mesh.element(:,k)); % Extract element endpoints
    uLocal = uGlobal(Fem.T(:,k)); % uLocal(1) = u1, uLocal(2) = u2, uLocal(3) = u1', uLocal(4) = u2'
    x = linspace(element(1), element(2), reso+1); % x contains the points we are evaluating the shape functions at
    xInterpolate(((k-1)*(reso+1)+1):(k*(reso+1))) = x;
    uInterpolate(((k-1)*(reso+1)+1):(k*(reso+1))) = evalfeherm1d(x, uLocal, ...
        element, hDegree, iDerivative); % uInterpolate = u1*N1 + u2*N2 + u1'*N3 + u2'*N4
    plot(x,uInterpolate(((k-1)*(reso+1)+1):(k*(reso+1)))); % Plot the result
end
hold off
