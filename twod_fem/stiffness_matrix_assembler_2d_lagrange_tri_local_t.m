function stiffnessMatrixLocal = stiffness_matrix_assembler_2d_lagrange_tri_local_t(aName, element, iDegree1, xDerivative1, ...
    yDerivative1, iDegree2, xDerivative2, yDerivative2, nQuadraturePoint,t)
%--------------------------------------------------------------------------
%stiffness_matrix_assembler_2d_lagrange_tri_local:
%--------------------------------------------------------------------------
qPoint = quadrature_node_generator_2d_triangle(element,nQuadraturePoint);
qWeight = quadrature_weight_generator_2d_triangle(element,nQuadraturePoint);
ldof1 = (iDegree1 + 1)*(iDegree1 + 2)/2; ldof2 = (iDegree2 + 1)*(iDegree2 + 2)/2;
stiffnessMatrixLocal = zeros(ldof1, ldof2);
if (nargin(aName) == 3)
    fVal = feval(aName,qPoint(1,:),qPoint(2,:),t);
end
if (nargin(aName) == 2)
    fVal = feval(aName,qPoint(1,:),qPoint(2,:));
end
for iShape1 = 1:ldof1
    for iShape2 = 1:ldof2
        integrand = shape_function_generator_2d_lagrange_triangle_local(qPoint(1,:), qPoint(2,:), element, iDegree1, xDerivative1, yDerivative1, iShape1).*...
            shape_function_generator_2d_lagrange_triangle_local(qPoint(1,:), qPoint(2,:), element, iDegree2, xDerivative2, yDerivative2, iShape2);
        %integrand = feval(aName,qPoint(1,:),qPoint(2,:)).*shape_function_generator_2d_lagrange_triangle_local(qPoint(1,:), qPoint(2,:), element, iDegree1, xDerivative1, yDerivative1, iShape1).*...
        %    shape_function_generator_2d_lagrange_triangle_local(qPoint(1,:), qPoint(2,:), element, iDegree2, xDerivative2, yDerivative2, iShape2);
        stiffnessMatrixLocal(iShape1,iShape2) = sum(qWeight.*fVal.*integrand);
    end
end
end









