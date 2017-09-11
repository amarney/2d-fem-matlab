function errorGlobal = error_2d_lagrange_global(uName, uGlobal, Mesh, Fem, xDerivative, yDerivative, nQuadraturePoint)
iDegree = Fem.degree;
%--------------------------------------------------------------------------
%error_2d_lagrange_global:
%--------------------------------------------------------------------------
errorGlobal = 0;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:,Mesh.element(:,k));
    uLocal = uGlobal(Fem.T(:,k));
    errorLocal = error_2d_lagrange_local(uName, uLocal, ...
        element, iDegree, xDerivative, yDerivative, nQuadraturePoint);
    errorGlobal = errorGlobal + errorLocal;
end
return;