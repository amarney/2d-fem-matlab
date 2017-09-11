function loadVectorLocal = load_vector_assembler_2d_lagrange_tri_local_t(fName, element, iDegree, ...
    xDerivative, yDerivative, nQuadraturePoint,t)
%--------------------------------------------------------------------------
%load_vector_assembler_2d_lagrange_tri_local:
%--------------------------------------------------------------------------
qPoint = quadrature_node_generator_2d_triangle(element, nQuadraturePoint);
qWeight = quadrature_weight_generator_2d_triangle(element,nQuadraturePoint);
ldof = (iDegree + 1)*(iDegree + 2)/2;
loadVectorLocal = zeros(ldof, 1);
fVal = feval(fName,qPoint(1,:),qPoint(2,:),t);
for iShape = 1:ldof
    integrand = shape_function_generator_2d_lagrange_triangle_local(qPoint(1,:), ...
        qPoint(2,:), element, iDegree, xDerivative, yDerivative, iShape);    
    %integrand = feval(fName,qPoint(1,:),qPoint(2,:)).*shape_function_generator_2d_lagrange_triangle_local(qPoint(1,:), ...
    %    qPoint(2,:), element, iDegree, xDerivative, yDerivative, iShape);
    loadVectorLocal(iShape) = sum(qWeight.*fVal.*integrand); 
end