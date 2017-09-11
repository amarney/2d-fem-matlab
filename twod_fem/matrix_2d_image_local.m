function matLoc = matrix_2d_image_local(aName, element, iDegree1, xDerivative1, ...
    yDerivative1, iDegree2, xDerivative2, yDerivative2, nQuadraturePoint)
%--------------------------------------------------------------------------
% Generates local stiffness matrix.
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
%--------------------------------------------------------------------------
qPoint = quadrature_node_generator_2d_triangle(element,nQuadraturePoint);
qWeight = quadrature_weight_generator_2d_triangle(element,nQuadraturePoint);
ldof1 = (iDegree1 + 1)*(iDegree1 + 2)/2; ldof2 = (iDegree2 + 1)*(iDegree2 + 2)/2;
matLoc = zeros(ldof1, ldof2);

% SET UP PARAMETERS FOR IMAGE/MATRIX/DATA: 
delta_x = 4.5e-6; delta_y = delta_x; nPixel = 90;
[X,Y] = meshgrid(0:delta_x:nPixel*delta_x, 0:delta_x:nPixel*delta_x);
fVal = interp2(X,Y,aName, qPoint(1,:), qPoint(2,:), 'cubic');


%fVal = feval(aName,qPoint(1,:),qPoint(2,:));
for iShape1 = 1:ldof1
    for iShape2 = 1:ldof2
        integrand = shape_function_generator_2d_lagrange_triangle_local(qPoint(1,:), qPoint(2,:), element, iDegree1, xDerivative1, yDerivative1, iShape1).*...
            shape_function_generator_2d_lagrange_triangle_local(qPoint(1,:), qPoint(2,:), element, iDegree2, xDerivative2, yDerivative2, iShape2);
        %integrand = feval(aName,qPoint(1,:),qPoint(2,:)).*shape_function_generator_2d_lagrange_triangle_local(qPoint(1,:), qPoint(2,:), element, iDegree1, xDerivative1, yDerivative1, iShape1).*...
        %    shape_function_generator_2d_lagrange_triangle_local(qPoint(1,:), qPoint(2,:), element, iDegree2, xDerivative2, yDerivative2, iShape2);
        matLoc(iShape1,iShape2) = sum(qWeight.*fVal.*integrand);
    end
end
end









