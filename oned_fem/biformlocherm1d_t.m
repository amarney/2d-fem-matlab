function matLoc = biformlocherm1d_t(aName, element, iDegree1, iDerivative1, ...
    iDegree2, iDerivative2, nGaussPoint, t)
gaussPoint = gauss_node_generator_1d_local(element,nGaussPoint);
gaussWeight = gauss_weight_generator_1d_local(element,nGaussPoint);
matLoc = zeros(iDegree1 + 3, iDegree2 + 3); %change +1 to +3
jaco = (element(2) - element(1))/2;
for iShape1 = 1:(iDegree1+3) % +1 to +3
    for iShape2 = 1:(iDegree2+3) % +1 to +3
        if ((iShape1 == 1) || (iShape1 == 2)) && ((iShape2 == 1) || (iShape2 == 2))
        integrand = feval(aName,gaussPoint,t).*shapeloc1dherm(gaussPoint,element,iDegree1,iDerivative1,iShape1).*...
            shapeloc1dherm(gaussPoint,element,iDegree2,iDerivative2,iShape2);
        matLoc(iShape1,iShape2) = sum(gaussWeight.*integrand);
        end
        if ((iShape1 == 3) || (iShape1 == 4)) && ((iShape2 == 1) || (iShape2 == 2))
        integrand = jaco*feval(aName,gaussPoint,t).*shapeloc1dherm(gaussPoint,element,iDegree1,iDerivative1,iShape1).*...
            shapeloc1dherm(gaussPoint,element,iDegree2,iDerivative2,iShape2);
        matLoc(iShape1,iShape2) = sum(gaussWeight.*integrand);            
        end
        if ((iShape1 == 1) || (iShape1 == 2)) && ((iShape2 == 3) || (iShape2 == 4))
        integrand = jaco*feval(aName,gaussPoint,t).*shapeloc1dherm(gaussPoint,element,iDegree1,iDerivative1,iShape1).*...
            shapeloc1dherm(gaussPoint,element,iDegree2,iDerivative2,iShape2);
        matLoc(iShape1,iShape2) = sum(gaussWeight.*integrand);            
        end      
        if ((iShape1 == 3) || (iShape1 == 4)) && ((iShape2 == 3) || (iShape2 == 4))
        integrand = jaco*jaco*feval(aName,gaussPoint,t).*shapeloc1dherm(gaussPoint,element,iDegree1,iDerivative1,iShape1).*...
            shapeloc1dherm(gaussPoint,element,iDegree2,iDerivative2,iShape2);
        matLoc(iShape1,iShape2) = sum(gaussWeight.*integrand);            
        end        
    end
end
end


