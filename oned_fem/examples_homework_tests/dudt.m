function dudt = dudt(u, q_0, N_00, M_22, bc_M22, bc_N00)
len_u = length(u)/2;
dudt(1:len_u) = u(1:len_u); 
dudt((len_u+1):(2*len_u)) = N_00\(q_0 + bc_M22 + bc_N00 - M_22*u(1:len_u));
end
