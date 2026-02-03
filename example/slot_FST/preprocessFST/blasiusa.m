function df=blasiusa(y,f,W)
df(1,1)=f(2);
df(2,1)=f(3);
df(3,1)=-W*f(1)*f(3);
