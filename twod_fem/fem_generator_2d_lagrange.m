function Fem = fem_generator_2d_lagrange(Mesh, iDegree)

% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu

mesh.p = Mesh.node; mesh.t = Mesh.element; 
mesh.e = Mesh.edge; M_ne = Mesh.M_ne;
if iDegree == 1
    p = Mesh.node; t = Mesh.element;
    Fem = struct('point',p,'T',t,'degree',iDegree);
elseif iDegree == 2
    p = zeros(2, length(mesh.p) + length(mesh.e));
    p(:, 1:length(mesh.p)) = mesh.p;
    t = zeros(6, length(mesh.t)); t(1:3, :) = mesh.t;
    p_count = length(mesh.p); e_ind = zeros(1, length(mesh.e));
    for k = 1:length(mesh.t)
        vert = mesh.p(:, mesh.t(:, k));
        i = 1; j = 2;
        if e_ind(M_ne(mesh.t(i, k), mesh.t(j, k))) == 0 % not treated
            p_count = p_count + 1;
            p(:, p_count) = (vert(:, i) + vert(:, j))/2;
            e_ind(M_ne(mesh.t(i, k), mesh.t(j, k))) = p_count;
            t(4, k) = p_count;
        else
            t(4, k) = e_ind(M_ne(mesh.t(i, k), mesh.t(j, k)));
        end
        i = 2; j = 3;
        if e_ind(M_ne(mesh.t(i, k), mesh.t(j, k))) == 0 % not treated
            p_count = p_count + 1;
            p(:, p_count) = (vert(:, i) + vert(:, j))/2;
            e_ind(M_ne(mesh.t(i, k), mesh.t(j, k))) = p_count;
            t(5, k) = p_count;
        else
            t(5, k) = e_ind(M_ne(mesh.t(i, k), mesh.t(j, k)));
        end
        i = 3; j = 1;
        if e_ind(M_ne(mesh.t(i, k), mesh.t(j, k))) == 0 % not treated
            p_count = p_count + 1;
            p(:, p_count) = (vert(:, i) + vert(:, j))/2;
            e_ind(M_ne(mesh.t(i, k), mesh.t(j, k))) = p_count;
            t(6, k) = p_count;
        else
            t(6, k) = e_ind(M_ne(mesh.t(i, k), mesh.t(j, k)));
        end
    end
    Fem = struct('point',p,'T',t,'degree',iDegree);
end
end