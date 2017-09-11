function qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint)
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
A1 = element(:,1); A2 = element(:,2); A3 = element(:,3);
tri_mat = [A1' 1; A2' 1; A3' 1];
K = det(tri_mat)/2; %maybe remove the /2?
if (nQuadraturePoint <= 0)
    fprintf('Please enter a positive value for n_g. \n');
    return
elseif nQuadraturePoint == 1
    qWeight = K;
elseif nQuadraturePoint == 3
    qWeight = [K/3, K/3, K/3];
elseif nQuadraturePoint == 7
    c1 = (155 - sqrt(15))/1200; c2 = (155 + sqrt(15))/1200;
    c3 = 9/40;
    qWeight = [K*c1, K*c1, K*c1, K*c2, K*c2, K*c2, K*c3];
end
end
