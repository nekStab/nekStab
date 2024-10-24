function [V,chk] = bcA(A,cc, Ny)
global gamma alpha D2 yvecs

B = ones(Ny,1)*(0+0i);
B(end)= 0+0i;
B(1) = 1+1i*0;    % v = 1 at the top
B(2)= cc;         % v = cc at the line under the top

[L,U] = lu(A);    % solve the system
yy = L\B;
V = U\yy;

% evaluate the error on the unbounded BC (Jacobs & Durbin 1998)
chk = (D2(2,:)*V+gamma^2*V(2))-exp(alpha*(yvecs(3)-yvecs(2)))*(D2(3,:)*V+gamma^2*V(3));
end