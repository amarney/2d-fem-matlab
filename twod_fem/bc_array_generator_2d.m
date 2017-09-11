function [dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName)
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
tol = 5*eps; % eps is the machine epsilon
B_val = feval(bName, Fem.point(1, :), Fem.point(2, :));
I_b = find(abs(B_val) <= tol); % which nodes are on the boundary
B_D_val = feval(bdName, Fem.point(1, I_b), Fem.point(2, I_b));
I_D = find(abs(B_D_val) <= tol); I_D = I_b(I_D);
I_N = find(abs(B_D_val) > tol); I_N = I_b(I_N);
nt = -200*ones(1, length(Fem.point));
nt(I_D) = zeros(1, length(I_D));
nt(I_N) = -ones(1, length(I_N));
dof_u = find(nt < 0);
end
