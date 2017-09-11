function vecGlob = linformglobherm1d_fsi(pressure, beamPosition, MeshFluid, MeshBeam, FemPressure, FemBeam, iDerivative, nGaussPoint)
% assemble load vector with hermite c1 elements
beamDegree = FemBeam.degree; pressureDegree = FemPressure.degree;
vecGlob = zeros(size(FemBeam.point, 2), 1);
pressureVals = zeros(1,nGaussPoint); % pressure values at gauss points
for k = 1:size(MeshBeam.element,2) %loop over all elements
    elementBeam = MeshBeam.node(MeshBeam.element(:,k));
    localBeam = beamPosition(FemBeam.T(:,k));
    gaussPoint = gauss_node_generator_1d_local(elementBeam,nGaussPoint); %x-vals
    gaussWeight = gauss_weight_generator_1d_local(elementBeam,nGaussPoint);
    beamPoint = evalfeherm1d(gaussPoint, localBeam, elementBeam, beamDegree, 0); % y-vals
    pt.x = gaussPoint; pt.y = beamPoint;
    for kk = 1:size(MeshFluid.element,2) 
        elementFluid = MeshFluid.node(:, MeshFluid.element(:,kk));
        point_is_in = ptInTriangle(pt, elementFluid);
        if (all(point_is_in(:))) %% it is important that MeshFluid and MeshBeam correspond at interface, so that ALL points are in
            localPressure = pressure(FemPressure.T(:,kk));
            pressureVals = evaluate_fe_function_2d_lagrange_tri(gaussPoint, beamPoint, localPressure, ...
                elementFluid, pressureDegree, 0, 0); % 0 for x and y derivatives
            break
        end
    end
    vecLoc = linformlocherm1d_fsi(pressureVals, elementBeam, beamDegree, iDerivative, gaussPoint, gaussWeight);
    vecGlob(FemBeam.T(:,k)) = vecGlob(FemBeam.T(:,k)) + vecLoc;
end
return;