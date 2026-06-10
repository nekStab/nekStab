function [D1,D2,D3,D4,y]=chebCL(N,X1,Xint)

N=N-1;
%%%%% Spatial Discretization
ksi=zeros(N,1);
for j=1:N+1
    ksi(j)=cos(pi*(j-1)/N);
end
a=Xint*X1/(X1-2*Xint);
b=1+2*a/X1;
y=a*(1+ksi)./(b-ksi);


c=zeros(1,N+1);
for i=2:N
c(i)=1;
end
c(N+1)=2;
c(1)=2;

%%%%%%%% Spectral numerical method

d=zeros(N+1,N+1);
for i=1:N+1
   for j=1:N+1
      if i~=j
      d(i,j)=((c(i)/c(j))*((-1)^(i+j)))/(ksi(i)-ksi(j));
      else
      d(i,j)=-ksi(i)/(2*(1-ksi(i)^2));
      end
   end
end
d(1,1)=(2*N.^2+1)/6;
d(N+1,N+1)=-(2*N.^2+1)/6;

%%%%%%% S dy/d(ksi) Matrix

S1=a*(b+1)./(y+a).^2;
S=diag(S1);

%%%%%%% physical derivatives

D1=S*d;
D2=D1*D1;
D3=D1*D2;
D4=D2*D2;
