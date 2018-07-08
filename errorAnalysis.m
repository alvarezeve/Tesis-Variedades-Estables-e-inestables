clear
format long eng

%==========================================================================
%Progaram Controlls:
%==========================================================================

%Course Map Computation Parameters:
    epsilon = 1.5;
    num_x = 25;               
    num_theta = 25;
    num_iterates = 100;

    %Note: num_x = 61;
    %      num_theta = 3;
    %      num_iterates = 175;
    %Parameters still fast...
   
%Parameterization of the Stable and Unstable manifolds:
    numManifoldPoints = 150;
    order = 1;
    endPoint = 2.5*10^(-5);
    step = (endPoint - 0)/numManifoldPoints;
    t = 0;
     
    

    
%==========================================================================
%Computations:
%==========================================================================
    
    
%Compute a course pic of the maps behavior:
this_x = 0.0;

for i = 1:num_x
    %i
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
lambda = 1 + epsilon/2 + sqrt(epsilon*(4+epsilon))/2;
alpha = 1;
lambda_bar = 1 + epsilon/2 - sqrt(epsilon*(4+epsilon))/2;
alpha_bar = alpha;

a1 = alpha;
b1 = a1*(1/2 + sqrt(epsilon*(4+epsilon))/(2*epsilon));

a_bar1 = alpha_bar;
b_bar1 = a_bar1*(1/2 - sqrt(epsilon*(4+epsilon))/(2*epsilon));

alpha0 = 1;
beta0 = 0;
alpha1 = 0;
beta1 = b1;


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


%Using the coefficents of the parameterization, compute the 
%stable and unstable manifolds:

figure
hold on 

colorMatrix(1:6, 1:3)=[1, 0, 0;
                       0, 0, 1;
                       0, 1, 0;
                       0, 0, 0;
                       1, 1, 0;
                       1, 0, 1];
colorMatrix(7:order, 1:3)= hot((order-6));   

parameterSpace = linspace(step, endPoint, numManifoldPoints);

for p = order:order
    p
    t = 0;
    for i = 1:numManifoldPoints
        t = t + (step);
        
        P_U(i, 1:2) = parmEval(t, A, B, p);
        P_S(i, 1:2) = parmEval(t, A_bar, B_bar, p);

        P_of_t_over_lambda = parmEval(t/lambda, A, B, p);
        f_composed_P = standardMap(P_of_t_over_lambda(1), P_of_t_over_lambda(2), epsilon);
        Error_U(i, 1:2) =  f_composed_P' - P_U(i,:);
        absError_U(i) = norm(Error_U(i,:));
    end    
    
   if (p <=5) 
       plot(parameterSpace(:), absError_U(:), '.','Color', colorMatrix(p,:));
   end
   
   if ((mod(p, 10) == 1) && (p >5)) 
       plot(parameterSpace(:), absError_U(:), '.','Color', colorMatrix(p,:));
   end  
   stable_manifold = P_S;
   unstable_manifold = P_U;
   P_U = 0;
   P_S = 0;

end

actualStop = endPoint

plot(parameterSpace(:), absError_U(:), '*r');

%==========================================================================
%Plotting:
%==========================================================================

figure 
hold on

%Plot Global Dynamics Sample:
plot(0,0, 'bo', 2*pi, 0, 'bo', 0, 2*pi, 'bo', 2*pi, 2*pi, 'bo');
plot(phasePlot(:, 1), phasePlot(:, 2), 'k.')

%Plot Parameterization of Unstable Manifold: 
plot(unstable_manifold(:,1), unstable_manifold(:,2), 'g.');
%plot(stable_manifold(:,1), stable_manifold(:,2), 'b.');
%plot(P_U2(:,1), P_U2(:,2), 'r.');

%Plot Parameterization of Stable Manifold:
%plot(PS(:,1), PS(:,2), 'r.');
%plot(PSM(:,1), PSM(:,2), 'y.');

axis([-0.4 2*pi*1.1 -0.4 2*pi*1.1])