load('mesh_rec_for_FE_mesh_viewer_2D.mat')
Mesh.node = mesh.p; Mesh.element = mesh.t;
elementLabel = 1; nodeLabel = 1; edgeLabel = 1;
figure(1); clf;
mesh_viewer_2d(Mesh, Mesh.node, elementLabel, nodeLabel, edgeLabel)
axis([-0.1, 1.1, -0.1, 1.1]);