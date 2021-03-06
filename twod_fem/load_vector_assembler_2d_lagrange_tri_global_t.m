function loadVectorGlobal = load_vector_assembler_2d_lagrange_tri_global_t(fName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint,t)
%--------------------------------------------------------------------------
%load_vector_assembler_2d_lagrange_tri_global:
%   assembles the global load vector.      
%--------------------------------------------------------------------------
iDegree = Fem.degree;
loadVectorGlobal = zeros(size(Fem.point, 2), 1);
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:,Mesh.element(:,k));
    loadVectorLocal = load_vector_assembler_2d_lagrange_tri_local_t(fName, element, iDegree, ...
        xDerivative, yDerivative, nQuadraturePoint,t);
    loadVectorGlobal(Fem.T(:,k)) = loadVectorGlobal(Fem.T(:,k)) + loadVectorLocal;
end
return;