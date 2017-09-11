function f = shape_function_generator_1d_lagrange_local(x, element, ...
    iDegree, iDerivative, iShape)
%--------------------------------------------------------------------------
%shape_function_generator_1d_lagrange_local:
%   returns the values (or derivative values) for the chosen local
%   finite-element shape function over the given element. Does so by 
%   performing an affine transformation from ``element'' to [-1,1], 
%   then calls ``shape_function_generator_1d_lagrange_reference'', and 
%   then doesa transformation back to original coordinates. 
%
%   Parameters  :   x - array consisting of coordinates where the shape
%                       function is being evaluated.
%                   element - array consisting of coordinates of nodes
%                             associated with a given element.
%                   iShape - index specifying which local
%                            shape function is evaluated.
%                   iDerivative - index specifying which derivative of the
%                                 shape function is evaluated.
%                   iDegree - index specifying the iDegree of the
%                            Lagrange polynomial being evaluated.
%
%   Return      :   f - array consisting of values (or derivative values)
%                       corresponding to the coordinates of t.
%   
%   Example     :
%       To make an example mesh with domain [0,3] and ten elements do
%       domain = [0,3]; nElement = 10;
%       Mesh = mesh_generator_1d(domain, nElement);
%       
%       To evaluate the first linear lagrange shape function  at x = 1.35, 
%       where x lies within the fifth element, do
%       x = 1.35; element = Mesh.node(Mesh.element(:,5));
%       iDegree = 1; iDerivative = 0; iShape = 1;
%       f = shape_function_generator_1d_lagrange_local(x, element, ...
%       iDegree, iDerivative, iShape)
%
%       To plot this shape function on the fifth element do
%       x = element(1):0.001:element(2);
%       f = shape_function_generator_1d_lagrange_local(x, element, ...
%       iDegree, iDerivative, iShape);
%       plot(x,f)
%
%       To plot all the shape functions for third degree polynomials on the
%       fifth element do
%       iDegree = 3; iDerivative = 0;
%       f1 = shape_function_generator_1d_lagrange_local(x, element, ...
%       iDegree, iDerivative, 1);
%       f2 = shape_function_generator_1d_lagrange_local(x, element, ...
%       iDegree, iDerivative, 2);
%       f3 = shape_function_generator_1d_lagrange_local(x, element, ...
%       iDegree, iDerivative, 3);
%       f4 = shape_function_generator_1d_lagrange_local(x, element, ...
%       iDegree, iDerivative, 4);
%       plot(x,f1,x,f2,x,f3,x,f4)
%
%   NOTE: Only goes up to third degree! 
%--------------------------------------------------------------------------
t = (2/(element(2)-element(1))).*(x-element(1)) - 1;
f = ((2/(element(2)-element(1)))^iDerivative)*shape_function_generator_1d_lagrange_reference(t, iDegree, iDerivative, iShape);
end