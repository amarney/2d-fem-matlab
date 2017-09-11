function gaussPoint = gauss_node_generator_1d_local(element, nGaussPoint)
%--------------------------------------------------------------------------
%gauss_node_generator_1d_local:
%   generates the gauss points locally on a given element. Generates  the
%   gauss points using ``gauss_node_generator_1d_reference'', and then maps
%   the points from the reference element to the given element.
%
%   Parameters  :   element - array consisting of coordinates of nodes
%                             associated with a given element.
%                   nGaussPoint - integer representing the number of gauss
%                                  to be generator.
%
%   Return      :   gaussPoint - array containing gaussian point values. 
%   
%   Example     :
%       To integrate cos(x)*(x^2 + 5x + 1) from 2 to 3 numerically do
%       nGaussPoint = 5; element = [2,3];
%       gaussPoint = gauss_node_generator_1d_local(element, nGaussPoint);
%       gaussWeight = gauss_weight_generator_1d_local(element, ...
%           nGaussPoint);
%       fName = @(x) cos(x).*(x.^2 + 5.*x + 1);
%       f = fName(gaussPoint);
%       intValue = sum(gaussWeight.*f)
%--------------------------------------------------------------------------
gaussPoint = gauss_node_generator_1d_reference(nGaussPoint);
gaussPoint = .5*(element(2)-element(1))*gaussPoint + .5*(element(1)+element(2));
end