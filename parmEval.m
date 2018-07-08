function P_of_t = parmEval(t, A, B, order)

x = 0;
theta = 0;

for k = 1:(order+1)
      x = x + A(k)*t^(k-1);
      theta = theta + B(k)*t^(k-1);
end

P_of_t = [mod(x, 2*pi), mod(theta, 2*pi)];