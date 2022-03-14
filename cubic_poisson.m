function [ y ] = cubic_poisson( Nx, zinit, w0, w1 )
   %[ solution ] = poisson_cubic( cardinality of the partition of (0,1),
   % initialization for fixed point algorithm (vector),
   % boundary value0, boundary value1 )
   
   % This function solve the problem
   % -\Delta y+y^3=h   (0,1)
   % y(0)=w0, y(1)=w1
   %For any z, consider y_z
   %-\Delta y_z+z^2y_z=h (0,1)
   %y_z(0)=w0, y_z(1)=w1.
   % Look for a fixed point of the map z to y_z

   deltax = 1/(Nx-1);
   spacepartition = linspace(0, 1, Nx);
   tol = 10^(-3);
   error = 10;
   Nmax = 10000;
   z = zinit;
   iter = 1;

   while (error > tol && iter < Nmax)
      iter = iter + 1;
      A = 2*eye(Nx-2,Nx-2);
      for i=1:Nx-3
         A(i+1,i) = -1;
         A(i,i+1) = -1;
      end
      A = A/((deltax)^2);
    
      tildeA = A+diag(z(2:Nx-1).^2);
      y=zeros(Nx,1);
      y(1) = w0;
      y(Nx) = w1;
    
      rhs = zeros(Nx-2,1);
      for i=1:Nx-2
         rhs(i) = 0; 
      end
      rhs(1) = rhs(1) + w0/((deltax)^2);
      rhs(Nx-2) = rhs(Nx-2) + w1/((deltax)^2);
      y(2:Nx-1)=(tildeA^(-1))*(rhs);
      z=(0.5)*z+(0.5)*y;
      
      ynormsquared=integral(@(x) interp1(spacepartition(1:Nx), y, x).^2, 0, 1);
      integrand = (A*y(2:Nx-1)+y(2:Nx-1).^3-rhs).^2;
      tildeerrornormsquared = sum(integrand*(deltax));
      
      if (ynormsquared<=10^(-1))
         error = sqrt(tildeerrornormsquared); 
      else
         error = sqrt(tildeerrornormsquared/ynormsquared);
      end
   end
end

