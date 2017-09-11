function vecLoc = linformlocherm1d_t(fName, element, iDegree, ...
    iDerivative, nGaussPoint,t)
gaussPoint = gauss_node_generator_1d_local(element,nGaussPoint);
gaussWeight = gauss_weight_generator_1d_local(element,nGaussPoint);
vecLoc = zeros(iDegree + 3, 1); % +1 to +3
jaco = (element(2) - element(1))/2;
for iShape = 1:(iDegree + 3) % +1 to +3
    if (iShape == 1) || (iShape == 2)
        integrand = feval(fName,gaussPoint,t).*shapeloc1dherm(gaussPoint, ...
            element, iDegree, iDerivative, iShape);
        vecLoc(iShape) = sum(gaussWeight.*integrand);
    end
    if (iShape == 3) || (iShape == 4)
        integrand = jaco*feval(fName,gaussPoint,t).*shapeloc1dherm(gaussPoint, ...
            element, iDegree, iDerivative, iShape);
        vecLoc(iShape) = sum(gaussWeight.*integrand);
    end
end