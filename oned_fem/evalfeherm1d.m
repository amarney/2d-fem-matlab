function u = evalfeherm1d(x, uLocal, ...
    element, iDegree, iDerivative)
% make uLocal a matrix, 2 by 2
% rows = iShape, cols = iiShape (derivs)
u = zeros(size(x));
jaco = (element(2) - element(1))/2;
for iShape = 1:(iDegree +3) %loop over shape index +1 -> +3
    if (iShape == 1) || (iShape == 2)
    u = u + uLocal(iShape)*shapeloc1dherm(x, ...
        element, iDegree, iDerivative, iShape);
    end    
    if (iShape == 3) || (iShape == 4)
        u = u + uLocal(iShape)*jaco*shapeloc1dherm(x, ...
        element, iDegree, iDerivative, iShape);
    end
end
end