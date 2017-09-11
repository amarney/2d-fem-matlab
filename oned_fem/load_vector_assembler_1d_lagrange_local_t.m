function loadVectorLocal = load_vector_assembler_1d_lagrange_local_t(fName, element, iDegree, ...
iDerivative, nGaussPoint, t)
%--------------------------------------------------------------------------
%load_vector_assembler_local_t:
%   computes the local load vector on a particular element by integrating
%   finite element basis functions against the external force fName.
%   Designed to be used in ``load_vector_assembler_global''.  Essentially
%   computes a discrete linear form.
%
%   Parameters  :   fName - function handle for external forces.
%                   uLocal - array containing local values of uGlobal on 
%                            ``element'',where uGlobal contains the 
%                            coefficients of the finite element solution.
%                   element - array consisting of coordinates of nodes
%                             associated with a given element.
%                   iDegree - index specifying the degree of the
%                            Lagrange polynomial being evaluated.
%                   iDerivative - index specifying which derivative of the
%                                 shape function is evaluated.                   
%                   nGaussPoint - integer representing the number of gauss
%                                  to be generator.
%
%   Return      :   loadVectorLocal - the local load vector associated with
%                                     external forces of a DE.
%   
%   Example     :
%       See readme      
%--------------------------------------------------------------------------
gaussPoint = gauss_node_generator_1d_local(element,nGaussPoint);
gaussWeight = gauss_weight_generator_1d_local(element,nGaussPoint);
loadVectorLocal = zeros(iDegree + 1, 1);
for iShape = 1:(iDegree+1)
    integrand = feval(fName,gaussPoint,t).*shape_function_generator_1d_lagrange_local(gaussPoint, ...
        element, iDegree, iDerivative, iShape);
    loadVectorLocal(iShape) = sum(gaussWeight.*integrand); 
end