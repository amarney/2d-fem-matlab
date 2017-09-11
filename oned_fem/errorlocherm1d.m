function errorLoc = errorlocherm1d(uName, uLocal, element, iDegree, iDerivative, nGaussPoint, t)
if nargin < 7
    gaussNode = gauss_node_generator_1d_local(element, nGaussPoint);
    gaussWeight = gauss_weight_generator_1d_local(element, nGaussPoint);
    integrand = (feval(uName, gaussNode)-evalfeherm1d(gaussNode, ...
        uLocal, element, iDegree, iDerivative)).^2;
    errorLoc = sum(gaussWeight.*integrand);
else
    gaussNode = gauss_node_generator_1d_local(element, nGaussPoint);
    gaussWeight = gauss_weight_generator_1d_local(element, nGaussPoint);
    integrand = (feval(uName, gaussNode, t)-evalfeherm1d(gaussNode, ...
        uLocal, element, iDegree, iDerivative)).^2;
    errorLoc = sum(gaussWeight.*integrand);
end
end