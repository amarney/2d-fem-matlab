function Fem = fem_generator_1d_lagrange(Mesh, iDegree)
%--------------------------------------------------------------------------
%fem_generator_1d_lagrange:
%   generates the Fem data structure, which is used to represent the finite
%   element space for use in finite element code.
%
%   Parameters  :   Mesh - the Mesh data structure which contains node
%                          coordinates, element connectivity information,
%                          and edge connectivity information.
%                   iDegree - degree associated with finite element space.
%
%   Return      :   Fem.point - an array containing the discretized
%                               x-coordinates of the finite element space.
%                   Fem.T - array containing connectivity information for
%                           the global degrees of freedom. The column-index
%                           corresponds to the element-index. The row-index
%                           corresponds to the local-shape-index (iShape).
%                           The entries tell you which coefficients of
%                           the local shape functions are nonzero. 
%                   Fem.degree - an integer giving the degree of
%                               polynomials used as basis functions for 
%                               the finite element space                             
%   
%   Example     :  
%       To create a sample Fem data structure do
%       domain = [0,1]; nElements = 5;
%       Mesh = mesh_generator_1d(domain,nElements);
%       degree = 2;
%       Fem = fem_generator_1d_lagrange(Mesh, iDegree);
%       
%       Given a finite element solution uGlobal, to extract uLocal for the
%       third element do
%       uLocal = uGlobal(Fem.T(:,3);
%--------------------------------------------------------------------------
nElement = size(Mesh.element,2);
t = zeros(iDegree+1,nElement); t(:,1) = (1:iDegree+1)';
for k=2:nElement
    t(:,k) = t(:,k-1)+iDegree;
end
p = zeros(1, (iDegree+1) + (nElement - 1)*iDegree);
k = 1;
element = Mesh.node(Mesh.element(:,k));
h = (element(2) - element(1))/iDegree;
p(1:iDegree+1) = element(1):h:element(2);
pointCount = iDegree + 1;
for k=2:nElement
    element = Mesh.node(Mesh.element(:,k));
    h = (element(2) - element(1))/iDegree;
    p(pointCount+1:pointCount+iDegree) = element(1)+h:h:element(2);
    pointCount = pointCount+iDegree;
end
Fem = struct('point', p, 'T', t, 'degree', iDegree);
end
