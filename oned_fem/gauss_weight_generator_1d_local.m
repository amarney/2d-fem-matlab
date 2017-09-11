function gaussWeight = gauss_weight_generator_1d_local(element, nGaussPoint)
%--------------------------------------------------------------------------
%gauss_weight_generator_1d_local:
%   generates the local gauss weights on a given element. Generates  the
%   gauss weights using ``gauss_weight_generator_1d_reference'', then maps
%   the weights from the reference element to the given element.
%
%   Parameters  :   element - array consisting of coordinates of nodes
%                             associated with a given element.
%                   nGaussPoint - integer representing the number of gauss
%                                  to be generator.
%
%   Return      :   gaussWeight - array containing gaussian point values. 
%   
%   Example     :
%       To integrate x^2 from 0 to 1 numerically do
%       nGaussPoint = 2; element = [0,1];
%       gaussPoint = gauss_node_generator_1d_local(element, nGaussPoint);
%       gaussWeight = gauss_weight_generator_1d_local(element, ...
%           nGaussPoint);
%       fName = @(x) x.^2;
%       f = fName(gaussPoint);
%       intValue = sum(gaussWeight.*f)
%--------------------------------------------------------------------------
gaussWeight = gauss_weight_generator_1d_reference(nGaussPoint);
gaussWeight = .5*(element(2)-element(1))*gaussWeight;
end