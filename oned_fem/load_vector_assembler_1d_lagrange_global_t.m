function loadVectorGlobal = load_vector_assembler_1d_lagrange_global_t(fName, Mesh, Fem, iDerivative, nGaussPoint, t)
%--------------------------------------------------------------------------
%load_vector_assembler_global_t:
%   assembles the global load vector. 
%
%   Parameters  :   fName - function handle for external forces.
%                   Mesh - the Mesh data structure which contains node
%                          coordinates, the element connectivity array,
%                          and the edge connectivity array.
%                   Fem - the Fem data structure which contains
%                         x-coordinates of of the finite element space, the
%                         connectivity array, and the degree of the basis
%                         polynomials.           
%                   element - array consisting of coordinates of nodes
%                             associated with a given element.
%                   iDerivative - index specifying which derivative of the
%                                 shape function is evaluated.                   
%                   nGaussPoint - integer representing the number of gauss
%                                  to be generator.
%
%   Return      :   loadVectorGlobal - the global load vector associated with
%                                     external forces of a DE.
%   
%   Example     :
%       See readme       
%--------------------------------------------------------------------------
iDegree = Fem.degree;
loadVectorGlobal = zeros(size(Fem.point, 2), 1);
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    loadVectorLocal = load_vector_assembler_1d_lagrange_local_t(fName, element, iDegree, ...
        iDerivative, nGaussPoint, t);
    loadVectorGlobal(Fem.T(:,k)) = loadVectorGlobal(Fem.T(:,k)) + loadVectorLocal;
end
return;