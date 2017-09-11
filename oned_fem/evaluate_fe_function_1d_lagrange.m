function u = evaluate_fe_function_1d_lagrange(x, uLocal, ...
    element, iDegree, iDerivative)
%--------------------------------------------------------------------------
%evaluate_fe_function_1d_lagrange:
%   locally evaluates a finite element function on an element using
%
%   u^{(r)}(x) = \sum_{i=1}^{p+1} uaLocal(i) * L_{K_k, i)^{(r)}(x)   
%
%   or equivalently
%
%   u^{(iDerivative)}(x) = \sum_{iShape=1}^{iDegree+1} uLocal(iShape) *
%                          L_{K_k, iShape}^{(iDerivative)}(x).
%
%   The value L_{K_k, i}^{(r)}(x) is just the r-th derivative of the i-th 
%   local lagrangian shape function on element k. Used for interpolating a
%   function with finite element basis functions (FE interpolation).
%
%   Parameters  :   x - array consisting of coordinates where the finite
%                       element function is being evaluated.
%                   uLocal - array containing local values of uGlobal on 
%                            ``element'',where uGlobal contains the 
%                            coefficients of the finite element solution.
%                   element - array consisting of coordinates of nodes
%                             associated with a given element.
%                   iDegree - index specifying the degree of the
%                            Lagrange polynomial being evaluated.
%                   iDerivative - index specifying which derivative of the
%                                 shape function is evaluated.
%
%   Return      :   u - array consisting of finite element function values 
%                       at points in x.
%   
%   Example     :   
%       To interpolate u(x) = sin(2*pi*x) with finite element basis
%       functions over [2,3] and plot the interpolation do
%
%       domain = [2,3]; nElement = 5;
%       Mesh = mesh_generator_1d(domain, nElement);
%       iDegree = 2; Fem = fem_generator_1d_lagrange(Mesh,iDegree);
%       uGlobal = sin(2*pi*Fem.point); iDerivative = 0;
%       nPlot = 20; figure(1); clf; hold on;
%       uInterpolate = zeros(1,(nPlot+1)*size(Mesh.element,2));
%       for k = 1:size(Mesh.element,2) %loop over element
%           element = Mesh.node(Mesh.element(:,k));
%           uLocal = uGlobal(Fem.T(:,k));
%           x = linspace(element(1), element(2), nPlot+1);
%           uInterpolate(((k-1)*(nPlot+1)+1):(k*(nPlot+1))) = evaluate_fe_function_1d_lagrange(x, uLocal, ...
%           element, iDegree, iDerivative);
%           plot(x,uInterpolate(((k-1)*(nPlot+1)+1):(k*(nPlot+1))));
%           pause
%       end
%
%       To plot the absolute error between our interpolation and u(x),
%       before the for loop include
%
%       uTrue = zeros(1,(nPlot+1)*size(Mesh.element,2));
%
%       and within the for loop include
%       
%       uTrue(((k-1)*(nPlot+1)+1):(k*(nPlot+1))) = sin(2*pi*x);
%       plot(x, abs(uTrue(((k-1)*(nPlot+1)+1):(k*(nPlot+1))) - uInterpolate(((k-1)*(nPlot+1)+1):(k*(nPlot+1)))), ...
%              'b', 'LineWidth', 2);
%       axis([2, 3, 0, 1])
%
%--------------------------------------------------------------------------
u = zeros(size(x));
for iShape = 1:(iDegree+1) %loop over shape index
    u = u + uLocal(iShape)*shape_function_generator_1d_lagrange_local(x, ...
        element, iDegree, iDerivative, iShape);
end
end