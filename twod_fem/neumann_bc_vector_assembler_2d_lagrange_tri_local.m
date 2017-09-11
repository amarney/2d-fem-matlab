function bc_Nloc = neumann_bc_vector_assembler_2d_lagrange_tri_local(gName, ...
    element, iDegree, xyEdgeNode, nQuadraturePoint)
edgeLength = norm(xyEdgeNode(:,2) - xyEdgeNode(:,1));
gaussPoint = gauss_node_generator_1d_reference(nQuadraturePoint); % [-1 1]
gaussWeight = gauss_weight_generator_1d_reference(nQuadraturePoint);
ldof = (iDegree + 1)*(iDegree + 2)/2;
bc_Nloc = zeros(ldof,1);
%f(X) = phi * g
%need F(t) = f(x(t),y(t)) replace x with x(t) and y with y(t)
%                         means evaluate gName at x(t) and y(t)
%    x(t) = (1-t)/2 *xyEdgeNode(1,1) + (1+t)/2 * xyEdgeNode(1,2)
%    y(t) = (1-t)/2 *xyEdgeNode(2,1) + (1+t)/2 * xyEdgeNode(2,2)
xt = xyEdgeNode(1,1).*(1 - gaussPoint)/2 + xyEdgeNode(1,2).*(1 + gaussPoint)/2;
yt = xyEdgeNode(2,1).*(1 - gaussPoint)/2 + xyEdgeNode(2,2).*(1 + gaussPoint)/2;
gNval = gName(xt,yt);
for iShape = 1:ldof
    integrand = gNval.*shape_function_generator_2d_lagrange_triangle_local(xt, ...
        yt, element, iDegree, 0, 0, iShape); %xderiv = 0, yderiv = 0
    bc_Nloc(iShape) = (edgeLength/2)*sum(gaussWeight.*integrand);
end
end