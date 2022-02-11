function [ y ] = poisson_cubic_nohomoboundary( Nx, zinit, w0, w1 )
%[ solution ] = poisson_cubic( cardinality of the partition of (0,1),
%initialization for fixed point algorithm (vector),
%boundary value0, boundary value1 )
%This function solve the problem
%-\Delta y+y^3=h       \mbox{in}  (0,1)
%y(0)=w0, y(1)=w1.
%We follow the apoproach of
%page 47 of \cite{boyer2013penalised}.
%By applying the fixed point method explained below.
%For any $z \in C^0([0,1];\mathbb{R})$,
%we define $y_z\in H^2(0,1)$ solution to
%-\Delta y_z+z^2y_z=h       \mbox{in}  (0,1)
%y_z(0)=w0, y_z(1)=w1.
%We look for a fixed point of the map
%\Gamma:C^0([0,1];\mathbb{R})\longrightarrow C^0([0,1];\mathbb{R})
%z\longmapsto y_z.

%Discretization in space.
deltax = 1/(Nx-1);

%Space partition.
spacepartition = linspace(0, 1, Nx);

%Tolerance.
tol=10^(-3);

%Initial error.
error = 10;

%maximal number of iterations.
Nmax=10000;

%STEP 1. Initialization.
%WARNING: "zinit"\in \mathbb{R}^{Nx}!!!
z=zinit;

%Iteration counter.
iter = 1;

while (error > tol && iter < Nmax)
    % Update iteration counter.
    iter = iter + 1;
    
    %STEP 3. Approximate solution to:
    %-\Delta y_z+z^2y_z=h       \mbox{in}  (0,1)
    %y_z(0)=y_z(1)=0.
    
    %We introduce the three diagonal Toeplitz matrix.
    A=2*eye(Nx-2,Nx-2);
    for i=1:Nx-3
       A(i+1,i)=-1;
       A(i,i+1)=-1;
    end
    A=A/((deltax)^2);
    
    %We introduce the matrix, which is
    %the discrete version of -\Delta +z^2 I_{L^2(0,1)}.
    tildeA=A+diag(z(2:Nx-1).^2);
    
    %We compute the approximate solution.
    y=zeros(Nx,1);
    y(1) = w0;
    y(Nx) = w1;
    
    rhs=zeros(Nx-2,1);
    for i=1:Nx-2
       rhs(i)=0; %h(spacepartition(i+1));
    end
    %WARNING. Differently from \cite[section 2.4, pages 15-16]{leveque2007finite},
    %in this code, the linear operator is -\Delta.
    rhs(1) = rhs(1) + w0/((deltax)^2);
    rhs(Nx-2) = rhs(Nx-2) + w1/((deltax)^2);
    
    y(2:Nx-1)=(tildeA^(-1))*(rhs);
    
    %STEP 4. We update z.
    z=(0.5)*z+(0.5)*y;
    
    %STEP 5. Computation of the error.
    
    %L^2 norm squared of "y".
    ynormsquared=integral(@(x) interp1(spacepartition(1:Nx), y, x).^2, 0, 1);
    
    %L^2 norm squared of "-\Delta y+y^3-h".
    integrand=(A*y(2:Nx-1)+y(2:Nx-1).^3-rhs).^2;
    
    tildeerrornormsquared=sum(integrand*(deltax));
    %tildeerrornormsquared=integral(@(x) interp1(spacepartition(2:Nx-1),
    %integrand, x), 0, 1);
    
    %STEP 8. Computation of the error.
    
    if (ynormsquared<=10^(-1))
       error = sqrt(tildeerrornormsquared); 
    else
       error = sqrt(tildeerrornormsquared/ynormsquared);
    end
    
    %STEP 9. We print on the "MATLAB's Command Window" the status of the algorithm.
    %fprintf("Iteration %i - Error %g\n", iter, error);
    
end


    %figure(1)
    %clf(1);
    %plot(spacepartition, y, 'r', 'LineWidth', 2);
    







end

