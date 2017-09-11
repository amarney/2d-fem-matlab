function u = evaluate_fe_function_2d_lagrange_tri(x, y, uLocal, ...
    element, iDegree, xDerivative, yDerivative)
%--------------------------------------------------------------------------
%evaluate_fe_function_2d_lagrange:
%--------------------------------------------------------------------------
u = zeros(size(x)); ldof = (iDegree + 1)*(iDegree + 2)/2;
for iShape = 1:ldof 
    u = u + uLocal(iShape)*shape_function_generator_2d_lagrange_triangle_local(x, y, ...
        element, iDegree, xDerivative, yDerivative, iShape);
end
end