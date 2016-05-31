function [psi,psi_xi]=SpaceBaseFunc_2(xi)

global N;

size = N+1;
psi = zeros(size,1);
psi_xi = zeros(size,1);
for i=1:size
   psi(i,1) = 0;
   psi_xi(i,1) = 0;
   for j=1:i
       psi(i,1) = psi(i,1) + (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-xi)^(j-1);
   end
   for j=2:i
       psi_xi(i,1) = psi_xi(i,1) - (j-1) * (-1)^(i-1) * nchoosek(i-1,j-1) * nchoosek(i+j-2,j-1) * (-xi)^(j-2);
   end
end

end