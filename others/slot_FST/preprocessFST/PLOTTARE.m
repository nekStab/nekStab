function PLOTTARE(A,col)

for i = 1:length(A(1,:))
    plot3(A(1,i),A(2,i),A(3,i),'ro','linewidth',3)
    txt = [' ' num2str(i)];
    text(A(1,i),A(2,i),A(3,i),txt);
    hold on
end

B=[A(:,13),A(:,1)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,13),A(:,14)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,2),A(:,14)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,2),A(:,17)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,1),A(:,17)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,20),A(:,17)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,20),A(:,5)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,9),A(:,5)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,9),A(:,1)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,9),A(:,10)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,3),A(:,10)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,3),A(:,13)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,5),A(:,15)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,7),A(:,15)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,7),A(:,10)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,7),A(:,19)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,18),A(:,19)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,18),A(:,3)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,18),A(:,4)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,14),A(:,4)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,4),A(:,12)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,11),A(:,12)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,11),A(:,2)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,15),A(:,16)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,8),A(:,16)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,8),A(:,19)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,16),A(:,6)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,8),A(:,12)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,11),A(:,6)];
plot3(B(1,:),B(2,:),B(3,:),col)
B=[A(:,20),A(:,6)];
plot3(B(1,:),B(2,:),B(3,:),col)

axis equal
% xlabel('$\omega$','interpreter','Latex','fontsize',20)
% ylabel('$\gamma$','interpreter','Latex','fontsize',20)
% zlabel('$\beta$','interpreter','Latex','fontsize',20)
























