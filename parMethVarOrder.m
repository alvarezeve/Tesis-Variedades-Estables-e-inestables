clear
format long

%==========================================================================
%Progaram Controlls:
%==========================================================================

%Course Map Computation Parameters:
    epsilon = 1.35;
    num_x = 25;               
    num_theta = 25;
    num_iterates = 200;

    %Note: num_x = 61;
    %      num_theta = 3;
    %      num_iterates = 175;
    %Parameters still fast...
   
%Parameterization of Unstable manifold:
    numManifoldPoints = 5000;
    order = 60;
    endPoint = 11.0;
    step = (endPoint - 0)/numManifoldPoints;
    t = 0;
    parmIterates = 1;    
    
%Unstable Manifold Computation Parameters:
    %berturbe from fixed point along manifild parameter: 
    %(ur variables are for the upper right hand branch of the manifolds)
    mu = 0.0000001
    num_manPoints = 100            %1200: many folds. 25: begin manifold.  
    num_manIterates = 15          %15 is a good value when epsilon is high.
    ur_corner = [2*pi; 2*pi];

    %NOTE:  When epsilon is low takt 'num_manPoints' low (say 100) and 
    %        'num_manIterates' high (say 50).
    %       When epsilon is high take 'num_manPoints' high (maybe 600-1000)
    %       but 'num_manIterates' low (say 10-15).
    

%Stable Manifold Computation Parameters:
    %perturbe from fixed point along manifild parameter: 
    %(ur variables are for the upper right hand branch of the manifolds)
    delta = mu
    num_smanPoints = num_manPoints
    num_smanIterates = num_manIterates
    lr_corner = [2*pi; 0];
    ul_corner = [0; 2*pi];


%Eigenvalues and Eigenvectors of Differential at Origin:
  %Numerical Values (test)  
  D0 = [1, epsilon; 1, 1+epsilon];
  [V, D] = eigs(D0);

  %Analytic Expressions:
    %eigenvalues:
    lambda1 = 1 + epsilon/2 + sqrt(epsilon*(4+epsilon))/2
    lambda2 = 1 + epsilon/2 - sqrt(epsilon*(4+epsilon))/2

    %Eigenvectors
    Xi1 = [-epsilon/2 + sqrt(epsilon*(epsilon+4))/2; 1];
    Xi2 = [-epsilon/2 - sqrt(epsilon*(epsilon + 4))/2; 1];

    %normalized eigenvectors:
    u1 = Xi1/norm(Xi1)
    u2 = Xi2/norm(Xi2)

%==========================================================================
%Computations:
%==========================================================================
    
    
%Compute a course pic of the maps behavior:
this_x = 0.0;

for i = 1:num_x
    i
    this_x = this_x + (i-1)*((2*pi)/(num_x - 1)); 
    this_theta = 0.0;
    thisThetaSection = [0, 0];
    
    for j = 1:num_theta
        this_theta = this_theta + (j-1)*((2*pi)/(num_theta - 1));
        x = this_x;
        theta = this_theta;
        thisIteration = [0, 0];
        
        for k = 1:num_iterates
            thisImage = standardMap(x, theta, epsilon);
            thisIteration(k, 1) = thisImage(1);
            thisIteration(k, 2) = thisImage(2);
            x = thisImage(1);
            theta = thisImage(2);
        end
        thisThetaSection((1+(j-1)*num_iterates):(j*num_iterates), 1:2) = thisIteration;
    end
    phasePlot((1 + num_iterates*num_theta*(i-1)):(num_iterates*num_theta*i), 1:2) = thisThetaSection;
end


%Unstable Manifold Brute Force Computation:
for i = 1:num_manPoints
    this_u = mu*(i)*u1;
    ur_u = ur_corner - this_u;
    x = this_u(1);
    theta = this_u(2);
    ur_x = ur_u(1);
    ur_theta = ur_u(2);
    for j = 1:num_manIterates
        thisImage = standardMap(x, theta, epsilon); 
        urImage = standardMap(ur_x, ur_theta, epsilon);
        USM(j, 1) = thisImage(1);
        USM(j, 2) = thisImage(2);
        urUSM(j, 1) = urImage(1);
        urUSM(j, 2) = urImage(2);
        x = thisImage(1);
        theta = thisImage(2);
        ur_x = urImage(1);
        ur_theta = urImage(2);
    end
    unstableManifold((1+(i-1)*num_manIterates):(i*num_manIterates), 1:2) = USM;
    ur_unstableManifold((1+(i-1)*num_manIterates):(i*num_manIterates), 1:2) = urUSM;
end

%Unstable manifold parameterization computation:
lambda = 1 + epsilon/2 + sqrt(epsilon*(4+epsilon))/2
alpha = 1

a1 = alpha
b1 = a1*(1/2 + sqrt(epsilon*(4+epsilon))/(2*epsilon))

alpha0 = 1
beta0 = 0
alpha1 = 0
beta1 = b1


A(2) = a1;
B(2) = b1;
Alpha(1) = 1;
Beta(1) = 0;
Alpha(2) = 0;
Beta(2) = b1;

%Compute the coefficents from the recursion formula:
for n = 1:(order-1)
   coef_a = -(epsilon/(n+1))*((1-lambda^(n+1))/((1-lambda^(n+1))*(1+epsilon-lambda^(n+1)) - epsilon));
   coef_b = (epsilon/(n+1))*((lambda^(n+1))/((1-lambda^(n+1))*(1+epsilon-lambda^(n+1)) - epsilon));
   
   thisTerm = 0;
   for k = 0:(n-1)
       thisTerm = thisTerm + (k+1)*Alpha(n-k+1)*B(k+2);
   end
   A(n+2) = coef_a*thisTerm;
   B(n+2) = coef_b*thisTerm;

   thisAlphaTerm = 0;
   thisBetaTerm = 0;
   
   for k = 0:n
       thisAlphaTerm = thisAlphaTerm + (k+1)*Beta(n-k+1)*B(k+2);
       thisBetaTerm = thisBetaTerm + (k+1)*Alpha(n-k+1)*B(k+2);
   end
   Alpha(n+2) = (-1/(n+1))*thisAlphaTerm;
   Beta(n+2) = (1/(n+1))*thisBetaTerm;
end

%Alpha'
%Beta'
%A'
%B'



%a3 = ((1-lambda^3)/((1-lambda^3)*(1+epsilon-lambda^3)-epsilon))*(epsilon/factorial(3))*b1^3
%b3 = ((-lambda^3)/((1-lambda^3)*(1+epsilon-lambda^3)-epsilon))*(epsilon/factorial(3))*b1^3

%a5 = ((1-lambda^5)/((1-lambda^5)*(1+epsilon-lambda^5)-epsilon))*epsilon*((3*b1^2*b3)/(factorial(3))-(b1^5)/(factorial(5)))
%b5 = ((-lambda^5)/((1-lambda^5)*(1+epsilon-lambda^5)-epsilon))*epsilon*((3*b1^2*b3)/(factorial(3))-(b1^5)/(factorial(5)))

%for i = 1:numManifoldPoints
%    t = t + step;
%    Px(i) = a1*t + a3*t^3 + a5*t^5;
%    Ptheta(i) = b1*t + b3*t^3 + b5*t^5;
%end    

%Using the coefficents of the parameterization, compute the 
%unstable manifold:
t = 0;
t_minus = 0;

for i = 1:numManifoldPoints
    t = t + step;
    t_minus = t_minus - step;
    %x = a1*t + a3*t^3 + a5*t^5;
    %theta = b1*t + b3*t^3 + b5*t^5;
    x = 0;
    theta = 0;
    x_minus = 0;
    theta_minus = 0;
    
    for k = 1:(order+1)
      x = x + A(k)*t^(k-1);
      theta = theta + B(k)*t^(k-1);
  
      x_minus = x_minus + A(k)*(t_minus)^(k-1);
      theta_minus = theta_minus + B(k)*(t_minus)^(k-1);
      
    end
    P_USM(1, 1) = x;
    P_USM(1, 2) = theta;
    
    P_minusUSM(1, 1) = x_minus + 2*pi;
    P_minusUSM(1, 2) = theta_minus + 2*pi;
    
    for j = 2:parmIterates
        thisImage = standardMap(x, theta, epsilon); 
        P_USM(j, 1) = thisImage(1);
        P_USM(j, 2) = thisImage(2);        
        x = thisImage(1);
        theta = thisImage(2);
    
        thisImage_minus = standardMap(P_minusUSM(j-1,1), P_minusUSM(j-1,2), epsilon); 
        P_minusUSM(j, 1) = thisImage_minus(1);
        P_minusUSM(j, 2) = thisImage_minus(2);        
        x_minus = thisImage_minus(1);
        theta_minus = thisImage_minus(2);
    end
    P((1+(i-1)*parmIterates):(i*parmIterates),1) = P_USM(:,1);
    P((1+(i-1)*parmIterates):(i*parmIterates),2) = P_USM(:,2);
    PM((1+(i-1)*parmIterates):(i*parmIterates),1) = P_minusUSM(:,1);
    PM((1+(i-1)*parmIterates):(i*parmIterates),2) = P_minusUSM(:,2);
end

%plot(Px(:), Ptheta(:), 'g.');


%==========================================================================
%Plotting:
%==========================================================================


figure 
hold on
plot(0,0, 'bo', 2*pi, 0, 'bo', 0, 2*pi, 'bo', 2*pi, 2*pi, 'bo');
plot(phasePlot(:, 1), phasePlot(:, 2), 'k.')

%plot(unstableManifold(:, 1), unstableManifold(:, 2), 'r.')
%plot(ur_unstableManifold(:, 1), ur_unstableManifold(:, 2), 'r.')

plot(P(:,1), P(:,2), 'g.');
plot(PM(:,1), PM(:,2), 'b.');

axis([-0.4 2*pi*1.1 -0.4 2*pi*1.1])