format long e
% FINAL PROJECT
%   Q.28.
%   a)
%%%%   1.0210077216510462818098013
%   1.352182434943101024230872
%   b)
f = @(x,y) exp(x-y/2);
xMin = 0; xMax = 1; yMin = 0; yMax = 1; 
error = zeros(3,1); 
tot_err1 = zeros(8,1); tot_err2 = zeros(8,1); tot_err3 = zeros(8,1);
for nx = 5:1:12
    ny = nx;
    Mesh = mesh_generator_2d_triangle(xMin, xMax, yMin, yMax, nx, ny);
	nQuadraturePoint = 1; integral_value = 0;
    for k = 1:size(Mesh.element,2) % loop over elements
        element = Mesh.node(:,Mesh.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
    end
    error(1) = abs(1.3521824349431010242308724 -integral_value);
	nQuadraturePoint = 3; integral_value = 0;
    for k = 1:size(Mesh.element,2) % loop over elements
        element = Mesh.node(:,Mesh.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
    end 
    error(2) = abs(1.3521824349431010242308724 -integral_value);
	nQuadraturePoint = 7; integral_value = 0;
    for k = 1:size(Mesh.element,2) % loop over elements
        element = Mesh.node(:,Mesh.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
    end  
    nx, error(3) = abs(1.3521824349431010242308724 -integral_value)
    tot_err1(nx-4) = error(1);  tot_err2(nx-4) = error(2); tot_err3(nx-4) = error(3);
end
% c)
%tot_err3(8) = 3.661338147750001e-16;
nEle = [5 6 7 8 9 10 11 12];
h = 1./nEle;
rc = polyfit(log(h), log(tot_err1)',1);
r1 = rc(1); c1 = exp(rc(2));
rc = polyfit(log(h), log(tot_err2)',1);
r3 = rc(1); c3 = exp(rc(2));
rc = polyfit(log(h), log(tot_err3)',1);
r7 = rc(1); c7 = exp(rc(2));
% d) no they dont, r7 should be 8, but its 6.
%   why not: NOT SURE, COME BACK LATER
clc, clear

% Q.29
% a)
% 3.658709154598590166803084
% b)
f = @(x,y) exp(x-y/2);
load('quadrature_on_circle_5_fields.mat')
mesh_1.node = mesh_1.p; mesh_1.element = mesh_1.t; mesh_1.edge = mesh_1.e;
nQuadraturePoint = 1; integral_value = 0; err1 = zeros(6,1); err3 = zeros(6,1);
for k = 1:size(mesh_1.element,2) % loop over elements
        element = mesh_1.node(:,mesh_1.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err1(1) = abs(3.658709154598590166803084 - integral_value);
nQuadraturePoint = 3; integral_value = 0;
for k = 1:size(mesh_1.element,2) % loop over elements
        element = mesh_1.node(:,mesh_1.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err3(1) = abs(3.658709154598590166803084 - integral_value);
mesh_2.node = mesh_2.p; mesh_2.element = mesh_2.t; mesh_2.edge = mesh_2.e;
nQuadraturePoint = 1; integral_value = 0;
for k = 1:size(mesh_2.element,2) % loop over elements
        element = mesh_2.node(:,mesh_2.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err1(2) = abs(3.658709154598590166803084 - integral_value);
nQuadraturePoint = 3; integral_value = 0;
for k = 1:size(mesh_2.element,2) % loop over elements
        element = mesh_2.node(:,mesh_2.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err3(2) = abs(3.658709154598590166803084 - integral_value);
mesh_3.node = mesh_3.p; mesh_3.element = mesh_3.t; mesh_3.edge = mesh_3.e;
nQuadraturePoint = 1; integral_value = 0;
for k = 1:size(mesh_3.element,2) % loop over elements
        element = mesh_3.node(:,mesh_3.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err1(3) = abs(3.658709154598590166803084 - integral_value);
nQuadraturePoint = 3; integral_value = 0;
for k = 1:size(mesh_3.element,2) % loop over elements
        element = mesh_3.node(:,mesh_3.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err3(3) = abs(3.658709154598590166803084 - integral_value);
mesh_4.node = mesh_4.p; mesh_4.element = mesh_4.t; mesh_4.edge = mesh_4.e;
nQuadraturePoint = 1; integral_value = 0;
for k = 1:size(mesh_4.element,2) % loop over elements
        element = mesh_4.node(:,mesh_4.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err1(4) = abs(3.658709154598590166803084 - integral_value);
nQuadraturePoint = 3; integral_value = 0;
for k = 1:size(mesh_4.element,2) % loop over elements
        element = mesh_4.node(:,mesh_4.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err3(4) = abs(3.658709154598590166803084 - integral_value);
mesh_5.node = mesh_5.p; mesh_5.element = mesh_5.t; mesh_5.edge = mesh_5.e;
nQuadraturePoint = 1; integral_value = 0;
for k = 1:size(mesh_5.element,2) % loop over elements
        element = mesh_5.node(:,mesh_5.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err1(5) = abs(3.658709154598590166803084 - integral_value);
nQuadraturePoint = 3; integral_value = 0;
for k = 1:size(mesh_5.element,2) % loop over elements
        element = mesh_5.node(:,mesh_5.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err3(5) = abs(3.658709154598590166803084 - integral_value);
mesh_6.node = mesh_6.p; mesh_6.element = mesh_6.t; mesh_6.edge = mesh_6.e;
nQuadraturePoint = 1; integral_value = 0;
for k = 1:size(mesh_6.element,2) % loop over elements
        element = mesh_6.node(:,mesh_6.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err1(6) = abs(3.658709154598590166803084 - integral_value)
nQuadraturePoint = 3; integral_value = 0;
for k = 1:size(mesh_6.element,2) % loop over elements
        element = mesh_6.node(:,mesh_6.element(:,k));
        qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
        qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
        f_val = f(qPoint(1,:), qPoint(2,:));
        integral_value = integral_value + sum(qWeight.*f_val);
end
err3(6) = abs(3.658709154598590166803084 - integral_value)

h_array = sparse(zeros(size(mesh_1.edge,2)));
for k = 1:size(mesh_1.edge,2) %loop over edge
    kEdge = mesh_1.node(:,mesh_1.edge(:,k));
    h_array(k) = norm(kEdge(:,2) - kEdge(:,1)); %kEdggeLength
end
h_1 = max(max(h_array'));

h_array = sparse(zeros(size(mesh_2.edge,2)));
for k = 1:size(mesh_2.edge,2) %loop over edge
    kEdge = mesh_2.node(:,mesh_2.edge(:,k));
    h_array(k) = norm(kEdge(:,2) - kEdge(:,1)); %kEdggeLength
end
h_2 = max(max(h_array'));

h_array = sparse(zeros(size(mesh_3.edge,2)));
for k = 1:size(mesh_3.edge,2) %loop over edge
    kEdge = mesh_3.node(:,mesh_3.edge(:,k));
    h_array(k) = norm(kEdge(:,2) - kEdge(:,1)); %kEdggeLength
end
h_3 = max(max(h_array'));

h_array = sparse(zeros(size(mesh_4.edge,2)));
for k = 1:size(mesh_4.edge,2) %loop over edge
    kEdge = mesh_4.node(:,mesh_4.edge(:,k));
    h_array(k) = norm(kEdge(:,2) - kEdge(:,1)); %kEdggeLength
end
h_4 = max(max(h_array'));

h = [h_1, h_2, h_3, h_4];
rc = polyfit(log(h), log(err1(1:4))',1);
r1 = rc(1); c1 = exp(rc(2));
rc = polyfit(log(h), log(err3(1:4))',1);
r3 = rc(1); c3 = exp(rc(2));
% h_array = sparse(zeros(size(mesh_5.edge,2)));
% for k = 1:size(mesh_5.edge,2) %loop over edge
%     kEdge = mesh_5.node(:,mesh_5.edge(:,k));
%     h_array(k) = norm(kEdge(:,2) - kEdge(:,1)); %kEdggeLength
% end
% h_5 = max(max(h_array'));
% 
% h_array = sparse(zeros(size(mesh_6.edge,2)));
% for k = 1:size(mesh_6.edge,2) %loop over edge
%     kEdge = mesh_6.node(:,mesh_6.edge(:,k));
%     h_array(k) = norm(kEdge(:,2) - kEdge(:,1)); %kEdggeLength
% end
% h_6 = max(max(h_array'));

%h = [h_1 h2 h3 h4 h5 h6]















