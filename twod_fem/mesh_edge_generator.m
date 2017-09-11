function [edge, M_e, M_ne] = mesh_edge_generator(Mesh, nConnectedComponent, nBoundary)
%--------------------------------------------------------------------------
%mesh_edge_generator:
%   generates the edge data structures for a given mesh. 
%
%   Parameters  :   Mesh - the Mesh data structure which contains node
%                          coordinates, the element connectivity array,
%                          and the edge connectivity array.
%                   nConnectedComponent - number of connected componenents      
%                   nBoundary - number of boundaries edges
%
%   Return      :   edge - array data structure where i-th column of `edge'
%                          contains the two indices for the nodes forming 
%                          the i-th edge.                     
%                   M_e - array data structure of size |N_h| by |N_h| which
%                         tells you if the i-j nodes form an edge shared by
%                         p elements. If M_e(i,j) = 1, then it is a
%                         boundary edge. If M_e(i,j) > 1, then it is an
%                         interior edge. If M_e(i,j) = 0 don't form an edge
%                   M_ne - array data structure which tells you if the i-j
%                          nodes form the p-th edge. If M_ne(i,j) = 0, then
%                          the i-j nodes don't form an edge. If M_ne = p,
%                          then  the i-j nodes form the p-th edge, p ~= 0.
%   Example     :
%       See readme.
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
%--------------------------------------------------------------------------
nEdge = length(Mesh.node) + length(Mesh.node) - (2*nConnectedComponent - nBoundary);
M_e = spalloc(length(Mesh.node), length(Mesh.node), nEdge);
for k = 1:size(Mesh.element,2)
    M_e(Mesh.element(:, k), Mesh.element(:, k)) = M_e(Mesh.element(:, k), Mesh.element(:, k)) ...
        + ones(3, 3);
end
[I, J, S] = find(M_e);
II = find(I<J); % indices of the 2 nodes on an edge is from small to large
edge = [I(II)'; J(II)'];
M_ne = sparse(I(II), J(II), 1:length(II), length(Mesh.node), length(Mesh.node));
M_ne = M_ne + M_ne';
end