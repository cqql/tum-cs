function [theta,theta_xi,theta_tau]=BaseFuncNodalTP_legendre(xi,tau)

global N;

size = (N+1)^2;

size_psi = N+1;

psi = zeros(size_psi,1);
psi_t = zeros(size_psi,1);
psi_xi = zeros(size_psi,1);
psi_tau = zeros(size_psi,1);

theta = zeros(size,1);
theta_xi = zeros(size,1);
theta_tau = zeros(size,1);

for i=1:size_psi
   psi(i,1) = 0;
   psi_t(i,1) = 0;
   psi_xi(i,1) = 0;
   psi_tau(i,1) = 0;
   for j=1:i
       psi(i,1) = psi(i,1) + (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-xi)^(j-1);
       psi_t(i,1) = psi_t(i,1) + (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-tau)^(j-1);
   end
   for j=2:i
       psi_xi(i,1) = psi_xi(i,1) - (j-1) * (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-xi)^(j-2);
       psi_tau(i,1) = psi_tau(i,1) - (j-1) * (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-tau)^(j-2);
   end
end
index = 1;

for i=1:size_psi
    for j=1:size_psi
        theta(index,1) = psi(i,1) * psi_t(j,1);
        theta_xi(index,1) = psi_xi(i,1) * psi_t(j,1);
        theta_tau(index,1) = psi(i,1) * psi_tau(j,1);
        index = index +1;
    end
end




end