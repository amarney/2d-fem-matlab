function stiffnessMatrixGlobal = stiffness_matrix_assembler_1d_lagrange_global(aName, Mesh, Fem1, ... 
    iDerivative1, Fem2, iDerivative2, nGaussPoint)
%--------------------------------------------------------------------------
%stiffness_matrix_assembler_lagrange_global:
%   assembles the global stiffnes matrix. 
%
%   Parameters  :   aName - function handle for DE coefficients.
%                   Mesh - the Mesh data structure which contains node
%                          coordinates, the element connectivity array,
%                          and the edge connectivity array.
%                   Fem1 - the Fem data structure for u(x); contains
%                          x-coordinates of of the finite element space, the
%                          connectivity array, and the degree of the basis
%                          polynomials.
%                   Fem2 - the Fem data structure for v(x); contains
%                          x-coordinates of of the finite element space, the
%                          connectivity array, and the degree of the basis
%                          polynomials.
%                   iDerivative1 - derivative index for u(x)   
%                   iDerivative2 - derivative index for v(x)
%                   nGaussPoint - integer representing the number of gauss
%                                 to be generator.
%
%   Return      :   stiffnessMatrixGlobal - the global stiffness matrix 
%                                           associated with coefficients of
%                                           DE.
%                                     
%   
%   Example     :
%       See readme       
%--------------------------------------------------------------------------
iDegree1 = Fem1.degree;
iDegree2 = Fem2.degree;
stiffnessMatrixGlobal = sparse(size(Fem1.point, 2), size(Fem2.point, 2)); %sparse
for k = 1:size(Mesh.element,2)
    elem = Mesh.node(Mesh.element(:,k));
    stiffnessMatrixLocal = stiffness_matrix_assembler_1d_lagrange_local(aName, elem, ...
        iDegree1, iDerivative1, iDegree2, iDerivative2, nGaussPoint);
    stiffnessMatrixGlobal(Fem1.T(:,k), Fem2.T(:,k)) = stiffnessMatrixGlobal(Fem1.T(:,k), Fem2.T(:,k)) + stiffnessMatrixLocal;
end
return;
end