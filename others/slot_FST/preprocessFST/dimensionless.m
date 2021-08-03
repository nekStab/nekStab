function dimensionless(numk,nmodes,delta,H)

mkdir RESULTATS2
for nfile = 1:numk*nmodes
    if nfile<10
       numero = [num2str(0) num2str(0) num2str(nfile)];
    elseif nfile<100
       numero = [num2str(0) num2str(nfile)];
    else
       numero = num2str(nfile);
    end
    name = ['./RESULTATS/wavenumber' numero '.dat'];
    name2 = ['./RESULTATS2/wavenumber' numero '.dat'];
    fid = fopen(name,'r');
    omega=fscanf(fid,'%f',1);  % pulsation
    gamma=fscanf(fid,'%f',1);  % wall-normal wavenumber 
    beta=fscanf(fid,'%f',1);   % spanwise wavenumber
    fclose(fid);
    %Re-dimensioneless!
    omega = omega*H/delta;
    gamma = gamma*H/delta;
    beta = beta*H/delta;
  
    fid = fopen(name2,'w');
    fprintf(fid,'%16f\n',omega);
    fprintf(fid,'%16f\n',gamma);
    fprintf(fid,'%16f\n',beta);
    fclose(fid);
end

for nfile = 1:numk*nmodes
    if nfile<10
       numero = [num2str(0) num2str(0) num2str(nfile)];
    elseif nfile<100
       numero = [num2str(0) num2str(nfile)];
    else
       numero = num2str(nfile);
    end
    name = ['./RESULTATS/velocity' numero '.dat'];
    name2 = ['./RESULTATS2/velocity' numero '.dat'];
    fid = fopen(name,'r');
    npoint=fscanf(fid,'%f',1);  
    VEL=fscanf(fid,'%f%f%f%f%f%f%f',[7 npoint]); 
    fclose(fid);
    VEL = VEL';
    VEL(:,1)= VEL(:,1)*delta/H;
    fid = fopen(name2,'w');
    fprintf(fid,'          %i\n',npoint);
    for ii = 1:npoint
        fprintf(fid,'%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e\n',VEL(ii,1),VEL(ii,2),VEL(ii,3),VEL(ii,4),VEL(ii,5),VEL(ii,6),VEL(ii,7));
    end
    fclose(fid);
end