function stiffnessMatrixGlobal = stiffness_matrix_assembler_2d_lagrange_tri_global_t(aName, Mesh, Fem1, ... 
    xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint,t)
%--------------------------------------------------------------------------
%stiffness_matrix_assembler_2d_lagrange_tri_global:  
% deriv1 is for u, deriv2 is for v (actually dont matter)
%--------------------------------------------------------------------------
iDegree1 = Fem1.degree;
iDegree2 = Fem2.degree;
stiffnessMatrixGlobal = sparse(size(Fem1.point, 2), size(Fem2.point, 2)); %sparse
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:,Mesh.element(:,k));
    stiffnessMatrixLocal = stiffness_matrix_assembler_2d_lagrange_tri_local_t(aName, element, ...
        iDegree1, xDerivative1, yDerivative1, iDegree2, xDerivative2, yDerivative2, nQuadraturePoint,t);
    stiffnessMatrixGlobal(Fem1.T(:,k), Fem2.T(:,k)) = stiffnessMatrixGlobal(Fem1.T(:,k), Fem2.T(:,k)) + stiffnessMatrixLocal;
end
return;
end