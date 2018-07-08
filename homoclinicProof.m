clear
format long

%==========================================================================
%Progaram Controlls:
%==========================================================================

%Course Map Computation Parameters:
    epsilon = 1.3;
    num_x = 3;               
    num_theta = 15;
    num_iterates = 100;

    %Note: num_x = 61;
    %      num_theta = 3;
    %      num_iterates = 175;
    %Parameters still fast...
   
%Parameterization of Unstable manifold:
    numManifoldPoints = 4000;
    order = 60;
    endPoint = 11.0;
    step = (endPoint - 0)/numManifoldPoints;
    t = 0;
    parmIterates = 1;    
    

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

%Stable and Unstable manifold parameterization computation:

%un barred are unstable and barred are stable:
lambda = 1 + epsilon/2 + sqrt(epsilon*(4+epsilon))/2
alpha = 1
lambda_bar = 1 + epsilon/2 - sqrt(epsilon*(4+epsilon))/2
alpha_bar = alpha

a1 = alpha
b1 = a1*(1/2 + sqrt(epsilon*(4+epsilon))/(2*epsilon))

a_bar1 = alpha_bar
b_bar1 = a_bar1*(1/2 - sqrt(epsilon*(4+epsilon))/(2*epsilon))

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

A_bar(2) = a_bar1;
B_bar(2) = b_bar1;
Alpha_bar(1) = 1;
Beta_bar(1) = 0;
Alpha_bar(2) = 0;
Beta_bar(2) = b_bar1;



%Compute the coefficents from the recursion formula:
for n = 1:(order-1)
   coef_a = -(epsilon/(n+1))*((1-lambda^(n+1))/((1-lambda^(n+1))*(1+epsilon-lambda^(n+1)) - epsilon));
   coef_b = (epsilon/(n+1))*((lambda^(n+1))/((1-lambda^(n+1))*(1+epsilon-lambda^(n+1)) - epsilon));
   
   coef_aBar = -(epsilon/(n+1))*((1-lambda_bar^(n+1))/((1-lambda_bar^(n+1))*(1+epsilon-lambda_bar^(n+1)) - epsilon));
   coef_bBar = (epsilon/(n+1))*((lambda_bar^(n+1))/((1-lambda_bar^(n+1))*(1+epsilon-lambda_bar^(n+1)) - epsilon));
   
   thisTerm = 0;
   thisTermBar = 0;
   for k = 0:(n-1)
       thisTerm = thisTerm + (k+1)*Alpha(n-k+1)*B(k+2);
       thisTermBar = thisTermBar + (k+1)*Alpha_bar(n-k+1)*B_bar(k+2);
   end
   A(n+2) = coef_a*thisTerm;
   B(n+2) = coef_b*thisTerm;

   A_bar(n+2) = coef_aBar*thisTermBar;
   B_bar(n+2) = coef_bBar*thisTermBar;
   
   thisAlphaTerm = 0;
   thisBetaTerm = 0;
   
   thisAlphaTermBar= 0;
   thisBetaTermBar = 0;
   
   for k = 0:n
       thisAlphaTerm = thisAlphaTerm + (k+1)*Beta(n-k+1)*B(k+2);
       thisBetaTerm = thisBetaTerm + (k+1)*Alpha(n-k+1)*B(k+2);
   
       thisAlphaTermBar = thisAlphaTermBar + (k+1)*Beta_bar(n-k+1)*B_bar(k+2);
       thisBetaTermBar = thisBetaTermBar + (k+1)*Alpha_bar(n-k+1)*B_bar(k+2);
   end
   Alpha(n+2) = (-1/(n+1))*thisAlphaTerm;
   Beta(n+2) = (1/(n+1))*thisBetaTerm;

   Alpha_bar(n+2) = (-1/(n+1))*thisAlphaTermBar;
   Beta_bar(n+2) = (1/(n+1))*thisBetaTermBar;
end

%Newton Method
t_0 = 5.0
t_n = t_0;

t_0check = parmEval(t_0, A, B, order)
t_0bar_check = parmEval(t_0, A_bar, B_bar, order)

for j = 1:6
   F = parmEval(t_n, A, B, order) - parmEval(t_n, A_bar, B_bar, order);
   Fprime = derivEval(t_n, A, B, order) - derivEval(t_n, A_bar, B_bar, order);
   F2 = F(2)
   Fprime2 = Fprime(2)   
   t_n_plus_1 = t_n - F2/Fprime2
   t_n = t_n_plus_1;
end    

t_n
int_check = parmEval(t_n, A, B, order) - parmEval(t_n, A_bar, B_bar, order) 

deriv_check1 = derivEval(t_n, A, B, order)
deriv_check2 = derivEval(t_n, A_bar, B_bar, order)
angle = acos(dot(deriv_check1, deriv_check2))

q_orbit(1, 1:2) = parmEval(t_n, A, B, order)
forwardNorms(1) = abs(norm(q_orbit(1, 1:2))-2*pi);


orbit_num = 15
for n = 1:orbit_num
   x = q_orbit(n,1);
   theta = q_orbit(n,2);
   q_orbit(n+1, 1:2) = standardMap(x, theta, epsilon);
   forwardNorms(n+1) = abs(norm(q_orbit(n+1, 1:2))-2*pi);
end

bq_orbit(1, 1:2) = parmEval(t_n, A, B, order)
backwardNorms(1) = norm(bq_orbit(1, 1:2));

for n = 1:orbit_num
   x = bq_orbit(n,1);
   theta = bq_orbit(n,2);
   bq_orbit(n+1, 1:2) = inverse_standardMap(x, theta, epsilon);
   backwardNorms(n+1) = norm(bq_orbit(n+1, 1:2));
end    



%==========================================================================
%Plotting:
%==========================================================================

figure 
hold on

%plot homoclinic orbit
plot(bq_orbit(:, 1), bq_orbit(:,2), 'or');
plot(q_orbit(:, 1), q_orbit(:,2), 'or');
plot(bq_orbit(:, 1), bq_orbit(:,2), '*r');
plot(q_orbit(:, 1), q_orbit(:,2), '*r');
plot(q_orbit(1, 1), q_orbit(1, 2), '*k')
plot(q_orbit(1, 1), q_orbit(1, 2), 'ok')

%Plot Global Dynamics Sample:
plot(0,0, 'bo', 2*pi, 0, 'bo', 0, 2*pi, 'bo', 2*pi, 2*pi, 'bo');
plot(phasePlot(:, 1), phasePlot(:, 2), 'k.')


axis([-0.4 2*pi*1.1 -0.4 2*pi*1.1])


figure
hold on

plot(forwardNorms(:), 'r')
plot(backwardNorms(:), 'b')

%plot(abs(q_orbit(:,1)), 'r')
%plot(abs(q_orbit(:,2)-2*pi), 'b')
%plot(abs(bq_orbit(:,1)), 'g')
%plot(abs(bq_orbit(:,2)), 'y')

