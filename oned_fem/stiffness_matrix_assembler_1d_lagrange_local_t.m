function stiffnessMatrixLocal = stiffness_matrix_assembler_1d_lagrange_local_t(aName, element, iDegree1, iDerivative1, ...
    iDegree2, iDerivative2, nGaussPoint, t)
%--------------------------------------------------------------------------
%stiffness_matrix_assembler_local_t:
%   computes the local stiffness matrix on a particular element by integrating
%   two finite element basis functions against the DE coefficient aName.
%   Designed to be used in ``stiffness_matrix_assembler_global''. Essentially
%   computes a discrete bilinear form.  
%
%   Parameters  :   aName - function handle for DE coefficients.
%                   element - array consisting of coordinates of nodes
%                             associated with a given element.
%                   iDegree1 - degree index for u(x).
%                   iDerivative1 - derivative index for u(x). 
%                   iDegree2 - degree index for v(x).
%                   iDerivative2 - derivative index for v(x).  
%                   nGaussPoint - integer representing the number of gauss
%                                  to be generator.
%
%   Return      :   stiffnessMatrixLocal - the local stiffness matrix 
%                                          associated with coefficients in
%                                          the DE.                                    
%   
%   Example     :
%       See readme      
%--------------------------------------------------------------------------
gaussPoint = gauss_node_generator_1d_local(element,nGaussPoint);
gaussWeight = gauss_weight_generator_1d_local(element,nGaussPoint);
stiffnessMatrixLocal = zeros(iDegree1 + 1, iDegree2 + 1);
for iShape1 = 1:(iDegree1+1)
    for iShape2 = 1:(iDegree2+1)
        integrand = feval(aName,gaussPoint,t).*shape_function_generator_1d_lagrange_local(gaussPoint,element,iDegree1,iDerivative1,iShape1).*...
            shape_function_generator_1d_lagrange_local(gaussPoint,element,iDegree2,iDerivative2,iShape2);
        stiffnessMatrixLocal(iShape1,iShape2) = sum(gaussWeight.*integrand);
    end
end
end









