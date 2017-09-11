function gaussWeight = gauss_weight_generator_1d_reference(nGaussPoint)
%--------------------------------------------------------------------------
%gauss_weight_generator_1d_reference:
%   generates the gauss weights on the reference element [-1,1]. Values
%   come from a formula which uses the gauss points.
%
%   Parameters  :   nGaussPoint - integer representing the number of gauss
%                                  to be generator.
%
%   Return      :   gaussWeight - array containing gaussian point values. 
%   
%   NOTE: Designed to be ``used in gauss_weight_generator_1d_local''! 
%--------------------------------------------------------------------------
if (nGaussPoint <= 0)
    fprintf('Please enter a positive value for n_g. \n');
    return
elseif (nGaussPoint == 1)
    gaussWeight = [2];
elseif (nGaussPoint == 2)
    gaussWeight = [1,1];
elseif (nGaussPoint == 3)
    gaussWeight = [0.5555555555555555556, 0.88888888888888888889, 0.5555555555555555556];
elseif (nGaussPoint == 4)
    gaussWeight = [0.3478548451374538574, 0.6521451548625461426, 0.6521451548625461426, ...
0.3478548451374538574];
elseif (nGaussPoint == 5)
    gaussWeight = [0.236926885056189088, 0.4786286704993664680, 0.56888888888888888889, ...
0.4786286704993664680, 0.236926885056189088];
else
    fprintf('The maximum nGaussPoint value is 5, please try again. \n');
end
end
    