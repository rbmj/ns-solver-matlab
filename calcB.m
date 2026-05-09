function B = calcB(P, e, gamma, m, n, rho)
  x0 = rho.^(-2);
  x1 = n.*x0;
  x2 = m.*x1;
  x3 = 1./rho;
  x4 = n.*x3;
  x5 = m.*x3;
  x6 = n.^2;
  x7 = gamma - 1;
  x8 = -e.*x7;
  x9 = -x7;
  B = [0 0 1 0; -x2 x4 x5 0; x0.*(rho.*(-P - x8) - x6) x5.*x9 x4.*(3 - gamma) x7; x1.*(-2*P - e - x8) x2.*x9 x0.*(rho.*(P + e) - x6.*x7) gamma.*x4];
end
