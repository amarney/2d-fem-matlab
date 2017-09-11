function essbc2Name_time = essbc2Name_time(x, y, beamVelocity, MeshBeam, wFem)
% x: points to evaluate at
% y: the bottom of the FluidMesh (old beam position)
% beamVelocity: the w_t(x)
% MeshBeam: the mesh of the beam over domain [0,1]
% wFem: the C1 cubic-hermite FE space for beam
%%length(x) = length(y)

% reso = 20;
% % wInterpolate = zeros(1,(reso+1)*size(MeshBeam.element,2));
% xInterpolate  = zeros(1,(reso+1)*size(MeshBeam.element,2));
% % for k = 1:size(MeshBeam.element,2)
% %     elementBeam = MeshBeam.node(MeshBeam.element(:,k));
% %     wLocal = beamDisplacement(wFem.T(:,k));
% %     xl = linspace(elementBeam(1), elementBeam(2), reso+1);
% %     xInterpolate(((k-1)*(reso+1)+1):(k*(reso+1))) = xl;
% %     wInterpolate(((k-1)*(reso+1)+1):(k*(reso+1))) = evalfeherm1d(xl, wLocal, ...
% %         elementBeam, 1, 0); % wInterp = u1*N1 + u2*N2 + u1'*N3 + u2'*N4
% % end
%
%
%
%
% wtInterpolate = zeros(1,(reso+1)*size(MeshBeam.element,2)); % bottom boundary
% for k = 1:size(MeshBeam.element,2)
%     elementBeam = MeshBeam.node(MeshBeam.element(:,k));
%     wLocal = beamVelocity(wFem.T(:,k));
%
%     %xl = linspace(elementBeam(1), elementBeam(2), reso+1);
%     wtInterpolate(((k-1)*(reso+1)+1):(k*(reso+1))) = evalfeherm1d(xl, wLocal, ...
%         elementBeam, 1, 0); % wInterp = u1*N1 + u2*N2 + u1'*N3 + u2'*N4
% end

% sort wInterpolate, and swap xInterpolate and wtInterpolate to match
%[SwInterp, IS] = sort(wInterpolate);
%SxInterp = xInterpolate(IS);
%SwtInterp = wtInterpolate(IS);

% eliminate repition
%[wInterpolate,IU,~] = unique(wInterpolate);
%xInterpolate = xInterpolate(IU);
%wtInterpolate = wtInterpolate(IU);

%[SCwInterpolate, SIC] = sort(CInterpolate);
%SCxInterpolate = CxInterpolate(SIC);
%SCwtInterpolate = CwtInterpolate(SIC);
for entry = 1:length(x)
    if ((x(entry) >= 0) & (x(entry) <= 1) & (y(entry) == 1)) % top
        essbc2Name_time(entry) = 0;
    elseif ((x(entry) == 0) & (y(entry) >= 0) & (y(entry) <= 1)) % left
        essbc2Name_time(entry) = 0;
    elseif ((x(entry) == 1) & (y(entry) >= 0) & (y(entry) <= 1)) % right
        essbc2Name_time(entry) = 0;
    else
        %% We are passing in x and y, but in reality the y value
        %  is not evaluated at (beam only moves in 1D).
        for k = 1:size(MeshBeam.element,2)
            elementBeam = MeshBeam.node(MeshBeam.element(:,k));
            if (x(entry) - elementBeam(1))*(x(entry) - elementBeam(2)) <= 0
                wtLocal = beamVelocity(wFem.T(:,k)); % we have found the correct element
                essbc2Name_time(entry) = evalfeherm1d(x(entry), wtLocal, elementBeam, 0, 0);
            end
        end
    end
    
    
    
    % else, on bottom
    %if kt == 1
    % essbc2Name_time = interp1(xInterpolate, wtInterpolate, x, 'pchip');
    %else
    % use griddatta instead of interp2
    %F = scatteredInterpolant(xInterpolate',wInterpolate',wtInterpolate','natural');
    % essbc2Name_time = F(x,y);
    %end
end
%% griddata returns a lot of NaN values, because x,y lies outside of convex region specified by
% xInterpolate, wInterpolate


end