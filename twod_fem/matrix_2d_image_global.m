function matGlob = matrix_2d_image_global(aName, Mesh, Fem1, ... 
    xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint)
%--------------------------------------------------------------------------
% Generates global stiffness matrix. 
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
%--------------------------------------------------------------------------
iDegree1 = Fem1.degree;
iDegree2 = Fem2.degree;
matGlob = sparse(size(Fem1.point, 2), size(Fem2.point, 2)); %sparse
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:,Mesh.element(:,k));
    matLoc = matrix_2d_image_local(aName, element, ...
        iDegree1, xDerivative1, yDerivative1, iDegree2, xDerivative2, yDerivative2, nQuadraturePoint);
    matGlob(Fem1.T(:,k), Fem2.T(:,k)) = matGlob(Fem1.T(:,k), Fem2.T(:,k)) + matLoc;
end
return;
end