function matGlob = biformglobherm1d(aName, Mesh, Fem1, iDerivative1, Fem2, iDerivative2, nGaussPoint)
% assemble stiffness matrix with hermite C1 elements
iDegree1 = Fem1.degree;
iDegree2 = Fem2.degree;
matGlob = sparse(size(Fem1.point, 2), size(Fem2.point, 2)); %sparse
for k = 1:size(Mesh.element,2) %loop over elements
    elem = Mesh.node(Mesh.element(:,k));
    matLoc = biformlocherm1d(aName, elem, ...
        iDegree1, iDerivative1, iDegree2, iDerivative2, nGaussPoint);
    matGlob(Fem1.T(:,k), Fem2.T(:,k)) = matGlob(Fem1.T(:,k), Fem2.T(:,k)) + matLoc;
end
return;
end