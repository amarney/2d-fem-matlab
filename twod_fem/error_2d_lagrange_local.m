function errorLocal = error_2d_lagrange_local(uName, uLocal, element, iDegree, xDerivative, yDerivative, nQuadraturePoint)
%--------------------------------------------------------------------------
%error_2d_lagrange_local:
%--------------------------------------------------------------------------
qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
qWeight = quadrature_weight_generator_2d_triangle(element, nQuadraturePoint);
integrand = (feval(uName, qPoint(1,:),qPoint(2,:))-evaluate_fe_function_2d_lagrange_tri(qPoint(1,:), qPoint(2,:), ...
    uLocal, element, iDegree, xDerivative, yDerivative)).^2;
errorLocal = sum(qWeight.*integrand);
return;