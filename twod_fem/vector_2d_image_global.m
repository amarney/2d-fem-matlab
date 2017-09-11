function vecGlob = vector_2d_image_global(fName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint)
%--------------------------------------------------------------------------
% Assembles the local load vector.
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
%--------------------------------------------------------------------------
iDegree = Fem.degree;
vecGlob = zeros(size(Fem.point, 2), 1);
for k = 1:size(Mesh.element,2)
    element = Mesh.node(:,Mesh.element(:,k));
    vecLoc = vector_2d_image_local(fName, element, iDegree, ...
        xDerivative, yDerivative, nQuadraturePoint);
    vecGlob(Fem.T(:,k)) = vecGlob(Fem.T(:,k)) + vecLoc;
end
return;