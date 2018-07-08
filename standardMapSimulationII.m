clear
format long

%==========================================================================
%Progaram Controlls:
%==========================================================================

%Course Map Computation Parameters:
    epsilon = 1.3;
    num_x = 25;               
    num_theta = 25;
    num_iterates = 150;

    %Note: num_x = 61;
    %      num_theta = 3;
    %      num_iterates = 175;
    %Parameters still fast...
    
%Unstable Manifold Computation Parameters:
    %berturbe from fixed point along manifild parameter: 
    %(ur variables are for the upper right hand branch of the manifolds)
    mu = 0.0000001
    num_manPoints = 1000            %1200: many folds. 25: begin manifold.  
    num_manIterates = 300          %15 is a good value when epsilon is high.
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


    
%Stable Manifold Brute Force Computation:
for i = 1:num_smanPoints
    this_s = delta*(i)*u2;
    lr_s = lr_corner + this_s;
    ul_s = ul_corner - this_s;
    ul_x = ul_s(1);
    ul_theta = ul_s(2);
    lr_x = lr_s(1);
    lr_theta = lr_s(2);
    for j = 1:num_smanIterates
        ulImage = inverse_standardMap(ul_x, ul_theta, epsilon); 
        lrImage = inverse_standardMap(lr_x, lr_theta, epsilon);
        ulSM(j, 1) = ulImage(1);
        ulSM(j, 2) = ulImage(2);
        lrSM(j, 1) = lrImage(1);
        lrSM(j, 2) = lrImage(2);
        ul_x = ulImage(1);
        ul_theta = ulImage(2);
        lr_x = lrImage(1);
        lr_theta = lrImage(2);
    end
    ul_stableManifold((1+(i-1)*num_manIterates):(i*num_manIterates), 1:2) = ulSM;
    lr_stableManifold((1+(i-1)*num_manIterates):(i*num_manIterates), 1:2) = lrSM;
end

sizeSM = size(lr_stableManifold);
pointCount = sizeSM(1,1);

x_0 = ul_stableManifold(pointCount,1);
theta_0 = ul_stableManifold(pointCount,2);
x = x_0;
theta = theta_0;
testOrbit(1, 1:2) = [x, theta];

for i = 2:450
      testOrbit(i, 1:2) = (standardMap(x, theta, epsilon))';
      x = testOrbit(i, 1);
      theta = testOrbit(i, 2);
      distOrig(i-1) = norm(testOrbit(i, 1:2) - [0, 2*pi]);
end    

plot(distOrig, '*');

%==========================================================================
%Plotting:
%==========================================================================


figure 
hold on
plot(0,0, 'bo', 2*pi, 0, 'bo', 0, 2*pi, 'bo', 2*pi, 2*pi, 'bo');
plot(phasePlot(:, 1), phasePlot(:, 2), 'k.')

plot(testOrbit(:,1), testOrbit(:,2), '*r');
plot(testOrbit(1,1), testOrbit(1,2), 'or');

%plot(unstableManifold(:, 1), unstableManifold(:, 2), 'r.')
%plot(ur_unstableManifold(:, 1), ur_unstableManifold(:, 2), 'r.')

%plot(ul_stableManifold(:, 1), ul_stableManifold(:, 2), 'g.')
%plot(lr_stableManifold(:, 1), lr_stableManifold(:, 2), 'g.')

axis([-0.4 2*pi*1.1 -0.4 2*pi*1.1])
