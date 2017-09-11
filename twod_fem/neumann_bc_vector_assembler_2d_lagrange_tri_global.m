function bc_N = neumann_bc_vector_assembler_2d_lagrange_tri_global(gName, ...
    Mesh, Fem, nQuadraturePoint)
iDegree = Fem.degree;
boundaryElement = find(Mesh.Kt == 1);
nBoundaryElement = length(boundaryElement);
bc_N = zeros(size(Fem.point,2),1); %maybe need to make smaller?....
for k = 1:nBoundaryElement % loop over boundary elements
    iElement = boundaryElement(k);
    iElementNode = Mesh.element(:,iElement);
    edge1 = Mesh.M_ne(iElementNode(1),iElementNode(2));
    edge2 = Mesh.M_ne(iElementNode(2),iElementNode(3));
    edge3 = Mesh.M_ne(iElementNode(3),iElementNode(1));
    et1 = Mesh.et(edge1); et2 = Mesh.et(edge2); et3 = Mesh.et(edge3);
    element = Mesh.node(:,Mesh.element(:,iElement));
    if et1 == -1 %NOTE, IT IS POSSIBLE TO HAVE 2 NEUMAN EDGES (BOTTOM RIGHT CORNER)
        %iEdgeNode = Mesh.edge(edge1); xyEdgeNode = Mesh.node(:,Mesh.edge(:,iEdgeNode));
        xyEdgeNode = Mesh.node(:,Mesh.edge(:,edge1));
        bc_Nloc = neumann_bc_vector_assembler_2d_lagrange_tri_local(gName, ...
            element, iDegree, xyEdgeNode, nQuadraturePoint);
    end
    if et2 == -1
        %iEdgeNode = Mesh.edge(edge2); xyEdgeNode = Mesh.node(:,Mesh.edge(:,iEdgeNode));
        xyEdgeNode = Mesh.node(:,Mesh.edge(:,edge2));
        bc_Nloc = neumann_bc_vector_assembler_2d_lagrange_tri_local(gName, ...
            element, iDegree, xyEdgeNode, nQuadraturePoint);
    end
    if et3 == -1
        %iEdgeNode = Mesh.edge(edge3); xyEdgeNode = Mesh.node(:,Mesh.edge(:,iEdgeNode));
        xyEdgeNode = Mesh.node(:,Mesh.edge(:,edge3));
        bc_Nloc = neumann_bc_vector_assembler_2d_lagrange_tri_local(gName, ...
            element, iDegree, xyEdgeNode, nQuadraturePoint);
    end
    if (et1 == -1 || et2 == -1 || et3 == -1)
        bc_N(Fem.T(:,iElement)) = bc_N(Fem.T(:,iElement)) + bc_Nloc;
    end
end
end
