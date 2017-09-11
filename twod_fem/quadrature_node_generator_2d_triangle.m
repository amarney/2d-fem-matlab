function qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint)
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
if (nQuadraturePoint <= 0)
    fprintf('Error: please enter a positive value for n_g \n');
    return
elseif nQuadraturePoint == 1
    qPoint = [1/3, 1/3, 1/3];
elseif nQuadraturePoint == 3
    qPoint = ...
        [.5 .5 0; ...
        0 .5 .5; ... 
        .5 0 .5];
elseif nQuadraturePoint == 7
    s1 = (9+2*sqrt(15))/21; t1 = (6-sqrt(15))/21;
    s2 = (9-2*sqrt(15))/21; t2 = (6+sqrt(15))/21;
    qPoint = [s1 t1 t1; t1 s1 t1; t1 t1 s1; s2 t2 t2; ...
                  t2 s2 t2; t2 t2 s2; 1/3 1/3 1/3];
end

% map back
A1 = element(:,1); A2 = element(:,2); A3 = element(:,3); 
invL = [1 1 1; A1 A2 A3]; 
one_x_y_vector = invL*qPoint';
qPoint = one_x_y_vector(2:3,:);
end