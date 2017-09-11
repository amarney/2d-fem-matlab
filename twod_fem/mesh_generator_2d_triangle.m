function Mesh = mesh_generator_2d_triangle(xMin, xMax, yMin, yMax, nx, ny)
%--------------------------------------------------------------------------
%mesh_generator_2d_triangle:
%   generates a Cartesian triangular mesh on a rectangle domain. 
%
%   Parameters  :   xMin - the minimum x value
%                   xMax - the maximum x value
%                   yMin - the minimum y value
%                   yMax - the maximum y value
%                   nx - the number of elements in the x side
%                   ny -  the number of elements in the y side
%
%   Return      :   Mesh.node - an array containing the nodes of the mesh.
%                   Mesh.element - an array containing element
%                                   connectivity information
%   
%   Example     :
%       % domain = [0,1] \times [0.5,1.1] with 5 elements in x side
%       % and 4 elements in y side:
%       xMin = 0; xMax = 1; yMin = 0.5; yMax = 1.1; nx = 5; ny = 4;
%       figure(1); clf; elem_label = 1; node_label = 1;
%       mesh_viewer_2d(Mesh, Mesh.node, 1, 1, 0)
%
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
%--------------------------------------------------------------------------
xNode = xMin:(xMax - xMin)/nx:xMax;
yNode = yMin:(yMax - yMin)/ny:yMax;

nnx = max(size(xNode));
nny = max(size(yNode));

% form nodes
node = zeros(2,nnx*nny); nodeCount = 1;
for j = 1:nny
    for k = 1:nnx
        node(:, nodeCount) = [xNode(k),yNode(j)];
        nodeCount = nodeCount + 1;
    end
end

% form elements
t = zeros(3, 2*nx*ny);
% form the first two element in the first rectangle
t(:,1) = [1; nnx+2; nnx+1]; t(:,2) = [1; 2; nnx+2]; t_count = 3;
for j=1:ny
    for i=2:nx % count rectangles in x direction
        t(:,t_count:t_count+1) = t(:,t_count - 2:t_count - 1) + 1;
        t_count = t_count + 2;
    end
    if j < ny % move the next row of rectangles
        t(:,t_count:t_count+1) = t(:,t_count - 2*nx:t_count - 2*nx+1) + nnx;
        t_count = t_count + 2;
    end
end
Mesh = struct('node', node, 'element', t);
end