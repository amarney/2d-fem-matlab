function errorLocal = error_1d_lagrange_local(uName, uLocal, element, iDegree, iDerivative, nGaussPoint)
%--------------------------------------------------------------------------
%error_1d_lagrange_local:
%   computes the local error on a given element using the L^2 norm,
%
%   norm(u - u_h)_{L^2(\K) =  \int_{K} (u^{(r)}(x) - u_h^{(r)}(x))^2 dx.   
%
%   Parameters  :   uName - function handle for true function.
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
%   Return      :   errorLocal - integer representing the local error on a
%                                 particular element.
%   
%   Example     :
%       To find the local error on the second element for u(x) = cos(x)
%       over a four-element mesh with domain [0,1] do
%
%       domain = [0,1]; nElement = 4;
%       Mesh = mesh_generator_1d(domain, nElement);
%       iDegree = 3; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
%       uName = @(x) cos(x); uGlobal = feval(uName, Fem.point);
%       iDerivative = 0; nGaussPoint = 5; k = 2;
%       element = Mesh.node(Mesh.element(:,k)); 
%       uLocal = uGlobal(Fem.T(:,k));
%       errorLocal = error_1d_lagrange_local(uName, uLocal, element, ...
%           iDegree, iDerivative, nGaussPoint)           
%--------------------------------------------------------------------------
gaussNode = gauss_node_generator_1d_local(element, nGaussPoint);
gaussWeight = gauss_weight_generator_1d_local(element, nGaussPoint);
integrand = (feval(uName, gaussNode)-evaluate_fe_function_1d_lagrange(gaussNode, ...
    uLocal, element, iDegree, iDerivative)).^2;
errorLocal = sum(gaussWeight.*integrand);
return;