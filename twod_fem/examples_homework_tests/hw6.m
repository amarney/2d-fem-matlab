% % Q 22:
 xMin = 0; xMax = 1; yMin = 0.5; yMax = 1.1; nx = 5; ny = 4;
 Mesh = mesh_generator_2d_triangle(xMin,xMax,yMin,yMax,nx,ny);
 mesh_viewer_2d(Mesh, Mesh.node, 1, 1, 1)

% % Q 23
load('mesh_rec_for_FE_mesh_viewer_2D.mat')
Mesh.node = mesh.p; Mesh.element = mesh.t;
k = 15;
element = Mesh.node(:, Mesh.element(:,k));
centerCoordinate = mean(element');
x = centerCoordinate(1); y = centerCoordinate(2);

iDegree = 1;
f = zeros(6,2);
f(1,1) = shape_function_generator_2d_lagrange_triangle_local(x, y, element, ...
    iDegree, 0, 0, 2);
f(2,1) = shape_function_generator_2d_lagrange_triangle_local(x, y, element, ...
    iDegree, 1, 0, 2);
f(3,1) = shape_function_generator_2d_lagrange_triangle_local(x, y, element, ...
    iDegree, 0, 1, 2);
f(4,1) = NaN;
f(5,1) = NaN;
f(6,1) = NaN;

iDegree = 2;
f(1,2) = shape_function_generator_2d_lagrange_triangle_local(x, y, element, ...
    iDegree, 0, 0, 2);
f(2,2) = shape_function_generator_2d_lagrange_triangle_local(x, y, element, ...
    iDegree, 1, 0, 2);
f(3,2) = shape_function_generator_2d_lagrange_triangle_local(x, y, element, ...
    iDegree, 0, 1, 2);
f(4,2) = shape_function_generator_2d_lagrange_triangle_local(x, y, element, ...
    iDegree, 0, 0, 6);
f(5,2) = shape_function_generator_2d_lagrange_triangle_local(x, y, element, ...
    iDegree, 1, 0, 6);
f(6,2) = shape_function_generator_2d_lagrange_triangle_local(x, y, element, ...
    iDegree, 0, 1, 6)

% % Q24
load('mesh_rec_for_FE_mesh_viewer_2D.mat')
Mesh.node = mesh.p; Mesh.element = mesh.t; iDegree = 2;
[edge, M_e, M_ne] = mesh_edge_generator(Mesh, 1, 1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;
Fem = fem_generator_2d_lagrange(Mesh, iDegree);
mesh_viewer_2d(Mesh, Fem.point, 1, 1, 0)

% % Q 25
f = @(x,y) cos(x + y.^2);
k = 14; element = Mesh.node(:,Mesh.element(:,k)); nQuadraturePoint = 3;
qNode = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
f_val = f(qNode(1,:), qNode(2,:));
Int1 = sum(qWeight.*f_val)
nQuadraturePoint = 7;
qNode = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
f_val = f(qNode(1,:), qNode(2,:));
Int2 = sum(qWeight.*f_val)



%f = @(x, y) sin(x) + cos(y);
%element = [0 1 0; 0 0 1]; nQuadraturePoint = 7;
%qnodes = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
%qweights = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
%f_val = f(qnodes(1,:), qnodes(2,:));
%Int = sum(qweights.*f_val)

