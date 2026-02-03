function [Blas,Blas_p,Blas_pp]=blasius_profil(H,yvecs)
alt=1;
err=1e-10;
a=0;
b=1.1;
W1=2;
W2=0.5;
W=1.480620611513586;

% while alt>=1E-12;
%    W=(W1+W2)/2;
%    arret=1;
%    a=0;
%    b=1.1;
%    while arret >=err;
%        c=(a+b)/2;
        c=0.571371930837631;
        f0=[0;0;c];
        [y,f]=ode45(@(y,f) blasiusa(y,f,W),[0 H],f0);
%         plot(f(:,2),y,'linewidth',1)
%        set(gca,'linewidth',1,'fontsize',18),% pause(.5), drawnow
%        if f(end,2)>1; b=c; else a=c; end
%        arret=abs(f(end,2)-1);
%    end
%    INT=simps(f(:,2),y);
%   if INT>1; W2=W; else W1=W; end
%   alt=abs(INT-1);
% end

f(:,4)=-W*f(:,1).*f(:,3);
INT=simps(f(:,2),y);
Blas=f(:,2);
Blas_p=f(:,3);
Blas_pp=f(:,4);

Blas=interp1(y,Blas,yvecs,'spline');
Blas_p=interp1(y,Blas_p,yvecs,'spline');
Blas_pp=interp1(y,Blas_pp,yvecs,'spline');


% disp('delta star*')
% disp(INT)
% disp('U_{inf}')
% disp(f(end,2))
format long
