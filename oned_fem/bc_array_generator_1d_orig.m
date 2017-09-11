function [uDof, nodeType] = bc_array_generator_1d_orig(Fem, domain, iBoundaryCondition)
nodeType = zeros(size(Fem.point));
%--------------------------------------------------------------------------
%boundary_condition_structure_generator_1d:
%   generates ``uDof'' and ``nodeType'' data structures for use in solving
%   boundary value problems.
%
%   Parameters  :   domain - a two-element array containing the
%                            coordinates of the endpoints.
%                   Fem - the Fem data structure which contains
%                         x-coordinates of of the finite element space, the
%                         connectivity array, and the degree of the basis
%                         polynomials.
%                   iBoundaryCondition - two-element array specifying the
%                                        boundary conditions. The entry
%                                        equals 0 if that boundary point is
%                                        essential/dirichlet, or equals -1
%                                        if it is natural/neumann.
%
%   Return      :   nodeType - array where a column-index corresponds to a
%                              a finite element index point (Fem.point). The
%                              entry equals 0 if that point is on the
%                              essential/dirichlet boundary condition. The
%                              entry equals -1 if that point is on the
%                              natural/neumann boundary condition. The
%                              entry equals -200 if that point is on the
%                              interior. 
%                   uDof - array where column-index corresponds to a finite
%                          element index point (Fem.point), except the
%                          essential/dirchlet boundary condition is
%                          skipped. The j-th point in N_{h,u}, where N_{h,u}
%                          is the set of finite element points excluding
%                          those on the essential/dirichlet boundary, is
%                          the uDof(j)-th point N_{h,g}, which is all of
%                          the finite element points (Fem.point).
%   
%   Example     :   
%       Consider the BVP
%           -(au')' + cu = f, x \in [0,1]
%           u(0) = g_0; a(1)u'(1) = g_1;
%       we can generate our data structures for this BVP as follows
%
%       domain = [0,1]; nElements = 3;
%       Mesh = mesh_generator_1d(domain, nElements);
%       iDegree = 2; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
%       iBoundaryCondition = [0,-1];
%       [uDof, nodeType] = bc_array_generator_1d(Fem, ...
%           domain, iBoundaryCondition);
%
%       This produces
%           uDof = [2, 3, 4, 5, 6, 7]
%           nodeType = [0, -200, -200, -200, -200, -200, -1]
%
%       The 1st point in the N_{h,u} is the uDof(1)th (2nd) point in Fem.point.
%       The only point not in uDof is the point not in N_{h,u}.
%--------------------------------------------------------------------------
tolerance = eps; % use Matlab’s tolerance
for i = 1:length(Fem.point)
    if (abs(Fem.point(i) - domain(1)) <= tolerance)
        nodeType(i) = iBoundaryCondition(1);
    elseif (abs(Fem.point(i) - domain(2)) <= tolerance)
        nodeType(i) = iBoundaryCondition(2);
    else
        nodeType(i) = -200;
    end
end
uDof = find(nodeType < 0);
end

