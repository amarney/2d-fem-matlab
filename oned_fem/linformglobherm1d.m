function vecGlob = linformglobherm1d(fName, Mesh, Fem, iDerivative, nGaussPoint)
% assemble load vector with hermite c1 elements
iDegree = Fem.degree;
vecGlob = zeros(size(Fem.point, 2), 1);
for k = 1:size(Mesh.element,2) %loop over all elements
    element = Mesh.node(Mesh.element(:,k));
    vecLoc = linformlocherm1d(fName, element, iDegree, ...
        iDerivative, nGaussPoint);
    vecGlob(Fem.T(:,k)) = vecGlob(Fem.T(:,k)) + vecLoc;
end
return;