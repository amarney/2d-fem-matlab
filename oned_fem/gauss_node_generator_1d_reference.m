function gaussPoint = gauss_node_generator_1d_reference(nGaussPoint)
%--------------------------------------------------------------------------
%gauss_node_generator_1d_reference:
%   generates the gauss nodes on the reference element [-1,1]. Values taken
%   from zeros of Legendre polynomial of degree nGaussPoints. 
%
%   Parameters  :   nGaussPoint - integer representing the number of gauss
%                                  to be generator.
%
%   Return      :   gaussPoint - array containing gaussian point values. 
%   
%   NOTE: Designed to be ``used in gauss_node_generator_1d_local''! 
%--------------------------------------------------------------------------
if (nGaussPoint <= 0)
    fprintf('Error: please enter a positive value for n_g \n');
    return
elseif (nGaussPoint == 1)
    gaussPoint = [0];
elseif (nGaussPoint == 2)
    gaussPoint = [-0.57735026918962576451, 0.57735026918962576451];
elseif (nGaussPoint == 3)
    gaussPoint = [-0.77459666924148337704, 0, 0.77459666924148337704];
elseif (nGaussPoint == 4)
    gaussPoint = [-0.86113631159405257522, -0.33998104358485626480, ...
0.33998104358485626480, 0.86113631159405257522];
elseif (nGaussPoint == 5)
    gaussPoint = [-0.90617984593866399280, -0.53846931010568309104, 0, ...
0.53846931010568309104, 0.90617984593866399280];
else
    fprintf('The maximum nGaussPoints value is 5, please try again. \n');
end    
end