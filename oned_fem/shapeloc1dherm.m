function f = shapeloc1dherm(x, element, iDegree, iDerivative, iShape)
% assume that iDegree = 1, hDegree = 3;
t = (2/(element(2)-element(1))).*(x-element(1)) - 1;
f = ((2/(element(2)-element(1)))^iDerivative)*shaperef1dherm(t, iDegree, iDerivative, iShape);
end