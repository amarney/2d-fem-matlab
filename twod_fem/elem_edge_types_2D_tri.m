function [Kt, et] = elem_edge_types_2D_tri(Mesh, bdName)
Kt = 2*ones(1, length(Mesh.element)); et = -200*ones(1, length(Mesh.edge));
for k = 1:length(Mesh.element)
    et_loc = [Mesh.M_e(Mesh.element(1, k), Mesh.element(2, k)), ...
        Mesh.M_e(Mesh.element(2, k), Mesh.element(3, k)), ...
        Mesh.M_e(Mesh.element(3, k), Mesh.element(1, k))];
    if min(et_loc) == 1
        Kt(k) = 1;          % boundary element
        elem = Mesh.node(:, Mesh.element(:, k));
        for i = 1:3 % check all 3 edges
            if et_loc(i) == 1
                X1 = elem(:, i); I = Mesh.element(i, k);
                if i <= 2
                    X2 = elem(:, i+1); J = Mesh.element(i+1, k);
                else
                    X2 = elem(:, 1); J = Mesh.element(1, k);
                end
                if (bdName(X1(1), X1(2)) == 0) ...
                        && (bdName(X2(1), X2(2)) == 0)
                    et(Mesh.M_ne(I, J)) = 0; % Dirichlet edge
                else
                    et(Mesh.M_ne(I, J)) = -1; % Neumann edge
                end
            end
        end
    end
end
% Kt = 2 -> interior element
% Kt = 1 -> boundary element
% et = 0 -> dirichlet edge
% et = -1 -> neumann edge
% et = -200 -> interior edge
end