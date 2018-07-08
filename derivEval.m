function P_prime = derivEval(t, A, B, order)

x = 0;
theta = 0;

for k = 2:(order+1)
      x = x + (k-1)*A(k)*t^(k-2);
      theta = theta + (k-1)*B(k)*t^(k-2);
end

P_prime = [x, theta];