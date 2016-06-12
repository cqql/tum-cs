function [theta,theta_xi,theta_eta,theta_tau]=BaseFuncNodalTP2D_legendre(xi,eta,tau)

global N;

size = (N+1)^3;

size_1 = N+1;

psi_tmp_xi = zeros(size_1,1);
psi_tmp_eta  = zeros(size_1,1);
psi_tmp_tau  = zeros(size_1,1);
psi_tmp_xixi  = zeros(size_1,1);
psi_tmp_etaeta  = zeros(size_1,1);
psi_tmp_tautau  = zeros(size_1,1);


theta = zeros(size,1);
theta_xi = zeros(size,1);
theta_eta = zeros(size,1);
theta_tau = zeros(size,1);

for i=1:size_1
   psi_tmp_xi(i,1) = 0;
   psi_tmp_xixi(i,1) = 0;
   for j=1:i
       psi_tmp_xi(i,1) = psi_tmp_xi(i,1) + (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-xi)^(j-1);
   end
   for j=2:i
       psi_tmp_xixi(i,1) = psi_tmp_xixi(i,1) - (j-1) * (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-xi)^(j-2);
   end
end

for i=1:size_1
   psi_tmp_eta(i,1) = 0;
   psi_tmp_etaeta(i,1) = 0;
   for j=1:i
       psi_tmp_eta(i,1) = psi_tmp_eta(i,1) + (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-eta)^(j-1);
   end
   for j=2:i
       psi_tmp_etaeta(i,1) = psi_tmp_etaeta(i,1) - (j-1) * (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-eta)^(j-2);
   end
end

for i=1:size_1
   psi_tmp_tau(i,1) = 0;
   psi_tmp_tautau(i,1) = 0;
   for j=1:i
       psi_tmp_tau(i,1) = psi_tmp_tau(i,1) + (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-tau)^(j-1);
   end
   for j=2:i
       psi_tmp_tautau(i,1) = psi_tmp_tautau(i,1) - (j-1) * (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-tau)^(j-2);
   end
end

index = 1;

for i=1:size_1
    for j=1:size_1
         for k=1:size_1
            theta(index,1) = psi_tmp_xi(i,1) * psi_tmp_eta(j,1) * psi_tmp_tau(k,1);
            theta_xi(index,1) = psi_tmp_xixi(i,1) * psi_tmp_eta(j,1) * psi_tmp_tau(k,1);
            theta_eta(index,1) = psi_tmp_xi(i,1) * psi_tmp_etaeta(j,1) * psi_tmp_tau(k,1);
            theta_tau(index,1) = psi_tmp_xi(i,1) * psi_tmp_eta(j,1) * psi_tmp_tautau(k,1);
            index = index +1;
         end
    end
end

end
















