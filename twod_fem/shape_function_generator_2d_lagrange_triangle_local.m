function f = shape_function_generator_2d_lagrange_triangle_local(x, y, element, ...
    iDegree, xDerivative, yDerivative, iShape)
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
A1 = element(:,1); A2 = element(:,2); A3 = element(:,3);
B = [A2 - A1, A3 - A1]; invB = inv(B);
hx = invB(1,1)*(x - A1(1)) + invB(1,2)*(y - A1(2));
hy = invB(2,1)*(x - A1(1)) + invB(2,2)*(y - A1(2));
if xDerivative == 0 && yDerivative == 0
    f = shape_function_generator_2d_lagrange_triangle_reference(hx,hy,iDegree,xDerivative,yDerivative,iShape);
elseif xDerivative == 1 && yDerivative == 0
    f = shape_function_generator_2d_lagrange_triangle_reference(hx,hy,iDegree,1,0,iShape)*invB(1,1) + ...
        shape_function_generator_2d_lagrange_triangle_reference(hx,hy,iDegree,0,1,iShape)*invB(2,1);
elseif xDerivative == 0 && yDerivative == 1
    f = shape_function_generator_2d_lagrange_triangle_reference(hx,hy,iDegree,1,0,iShape)*invB(1,2) + ...
        shape_function_generator_2d_lagrange_triangle_reference(hx,hy,iDegree,0,1,iShape)*invB(2,2);
end 
end