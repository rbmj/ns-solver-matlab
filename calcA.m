function A = calcA(P, e, gamma, m, n, rho)
  x0 = rho.^(-2);
  x1 = m.^2;
  x2 = gamma - 1;
  x3 = -e.*x2;
  x4 = 1./rho;
  x5 = m.*x4;
  x6 = -x2;
  x7 = n.*x4;
  x8 = m.*x0;
  x9 = n.*x8;
  A = [0 1 0 0; x0.*(rho.*(-P - x3) - x1) x5.*(3 - gamma) x6.*x7 x2; -x9 x7 x5 0; x8.*(-2*P - e - x3) x0.*(rho.*(P + e) - x1.*x2) x6.*x9 gamma.*x5];
end
