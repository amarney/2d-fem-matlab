function errorGlob = errorglobherm1d(uName, uGlobal, Mesh, Fem, iDerivative, nGaussPoint, t)
iDegree = Fem.degree;
errorGlob = 0;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    uLocal = uGlobal(Fem.T(:,k));
    if nargin < 7
        errorLocal = errorlocherm1d(uName, uLocal, ...
            element, iDegree, iDerivative, nGaussPoint);
    else
        errorLocal = errorlocherm1d(uName, uLocal, ...
            element, iDegree, iDerivative, nGaussPoint, t);
    end
    errorGlob = errorGlob + errorLocal;
end
return;