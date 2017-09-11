function F = beam_semi_FE_1d(t, u, Mesh, Fem, qName, bc_eM22, bc_eN00, N_00, M_22)
qh_0 = load_vector_assembler_1d_lagrange_global_t(qName, Mesh, Fem, 0, 5, t);
q_0 = qh_0(Fem.dof_u);
F = dudt(u, q_0, N_00, M_22, bc_eM22, bc_eN00);
F = F';
end