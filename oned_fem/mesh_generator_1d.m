function Mesh = mesh_generator_1d(domain, nElement)
%--------------------------------------------------------------------------
%mesh_generator_1d:
%   generates an evenly-spaced one-dimensional Mesh data structure for use
%   in finite element code. 
%
%   Parameters  :   domain - a two-element array containing the
%                            coordinates of the endpoints.
%                   nElement - an integer representing the number of
%                               elements.
%
%   Return      :   Mesh.node - an array containing the nodes of the mesh.
%                   Mesh.element - an array containing element
%                                   connectivity information. The
%                                   column-index corresponds to the
%                                   element-index. The first (second) row
%                                   is the index for that element's left
%                                   (right) endpoint.
%                  Mesh.edge - an array containing edge connectivity
%                                information. The first row contains
%                                indices of the edges. The second (third)
%                                row contains the index of the element on
%                                the left (right) for a particular edge.
%                                An element-index of 0 means the edge is 
%                                at the boundary. Note that in the
%                                one-dimensional case the edges are just
%                                the nodes.
%   
%   Example     :
%       To generate a five-element mesh over the domain [0,1] do
%       domain = [0,1]; nElement = 5;
%       Mesh = mesh_generator_1d(domain, nElement);
%
%       To extract the coordinate of the sixth node do
%       Mesh.node(6)
%
%       To extract the nodal-indices associated with the fifth element do
%       Mesh.element(:,5)
%       
%       To extract the coordinates for the nodes of the fifth element do
%       element = Mesh.node(Mesh.element(:,5))
%
%       To find the element-indices associated with the sixth edge do
%       Mesh.edge(2:3,6)
%       
%--------------------------------------------------------------------------
x_l = domain(1); x_r = domain(2);
p = x_l:(x_r - x_l)/nElement:x_r;
t = zeros(2,nElement); t(:,1) = [1;2];
for i = 2:nElement
    t(:,i) = t(:,i-1) + 1;
end
e = zeros(3,nElement+1); e(1,:) = 1:nElement+1; e(2,:) = 0:nElement;
e(3,:) = [1:nElement,0];
Mesh = struct('node',p,'element',t,'edge',e);
end

