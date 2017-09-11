function mesh_viewer_2d(Mesh, nodes, elementLabel, nodeLabel, edgeLabel)
%--------------------------------------------------------------------------
%mesh_viewer_2d:
%   displays a plot of the 2d mesh. 
%
%   Parameters  :   Mesh - the Mesh data structure which contains node
%                          coordinates, the element connectivity array,
%                          and the edge connectivity array.
%                   nodes - the nodes to be plotted, typically Mesh.node         
%                   elementLabel - 0 for no labels, 1 for labels
%                   nodeLabel - 0 for no labels, 1 for labels
%                   edgeLabel - 0 for no labels, 1 for labels
%
%   Return      :   displays plot
%   
%   Example     :
%       load('meshdata');
%       mesh_viewer_2d(Mesh, Mesh.node, 1, 1, 1)    
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
%--------------------------------------------------------------------------
triplot(Mesh.element', nodes(1,:)', nodes(2,:)','b');
hold on;
if (nodeLabel == 1)
    for k = 1:size(nodes,2)
        %caption = sprintf('%d',k);
        plot(nodes(1,k), nodes(2,k), 'k.');
        text(nodes(1,k), nodes(2,k), int2str(k));
    end
end
if (elementLabel == 1)
    for k = 1:size(Mesh.element,2)
        elementCoordinate = Mesh.node(:,Mesh.element(:,k));
        centerCoordinate = mean(elementCoordinate');
        textStruct = text(centerCoordinate(1), centerCoordinate(2), int2str(k));
        textStruct.Color = 'r';
    end
end
if (edgeLabel == 1)
    nConnectedComponent = 1; nBoundary = 1;
    [edge, M_e, M_ne] = mesh_edge_generator(Mesh, nConnectedComponent, nBoundary);
    for k = 1:size(edge,2)
        edgeCoordinate = Mesh.node(:,edge(:,k));
        centerCoordinate = mean(edgeCoordinate');
        text(centerCoordinate(1), centerCoordinate(2), int2str(k));
    end
end
title('Mesh')
end