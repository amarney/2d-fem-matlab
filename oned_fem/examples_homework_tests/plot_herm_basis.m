t = -1:0.01:1; iDegree = 1; iDerivative = 0;
f1 = shaperef1dherm(t, iDegree, iDerivative, 1)
f2 = shaperef1dherm(t, iDegree, iDerivative, 2)
f3 = shaperef1dherm(t, iDegree, iDerivative, 3)
f4 = shaperef1dherm(t, iDegree, iDerivative, 4)
%plot(t,f1,t,f2,t,f3,t,f4)

iDegree = 3
      f1 = shape_function_generator_1d_lagrange_reference(t, ... 
      iDegree, iDerivative, 1);
      f2 = shape_function_generator_1d_lagrange_reference(t, ... 
      iDegree, iDerivative, 2);
      f3 = shape_function_generator_1d_lagrange_reference(t, ... 
      iDegree, iDerivative, 3);
      f4 = shape_function_generator_1d_lagrange_reference(t, ... 
      iDegree, iDerivative, 4);
      plot(t,f1,t,f2,t,f3,t,f4)