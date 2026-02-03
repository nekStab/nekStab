function FST_modes(Re, Ny, Ly,plotres,numk,nmode)

global gamma kk D2 yvecs alpha omega

Lreal = Ly;
y = linspace(0,Lreal,Ny);

for nfile = 1:numk*nmode
    fprintf('FILE NUMBER = %i \n',nfile);
    if nfile<10
       numero = [num2str(0) num2str(0) num2str(nfile)];
    elseif nfile<100
       numero = [num2str(0) num2str(nfile)];
    else
       numero = num2str(nfile);
    end
    name = ['./FST_data/wavenumber' numero '.dat'];
    fprintf('Opening file: %s', name)
    fid = fopen(name,'r');
    if fid == -1
        error('Author:Function:OpenFile', 'Cannot open file: %s', name);
    end
    omega=fscanf(fid,'%f',1);  % pulsation
    fclose(fid);
    
    fid = fopen(name,'r');
    if fid == -1
        error('Author:Function:OpenFile', 'Cannot open file: %s', name);
    end
    gamma=fscanf(fid,'%f',1);  % wall-normal wavenumber 
    fclose(fid);
    
    fid = fopen(name,'r');
    if fid == -1
        error('Author:Function:OpenFile', 'Cannot open file: %s', name);
    end
    beta=fscanf(fid,'%f',1);   % spanwise wavenumber
    fclose(fid);
    
    CASE = 'TEMPORAL';  % if SPATIAL enter real omega and find complex alpha
                        % if TEMPORAL enter real alpha and find complex omega
    
    switch CASE
        case 'SPATIAL'
            alpha = 1i/2*(sqrt(Re^2+4*(beta^2+gamma^2)-4i*omega*Re)-Re); %(Jacobs 2001)
            c = omega/alpha;
            kk = sqrt(alpha^2+beta^2);
        case 'TEMPORAL'
            alpha = omega;       % Taylor approximation ==> alpha = omega*U_{inf} 
            kk = sqrt(alpha^2+beta^2);
            c = 1-1i*(1+gamma^2/kk^2)*kk^2/(alpha*Re); %(Brandt 2002)
    end
    
    %DERIVATION MATRICES
    [D1,D2,~,D4,yvecs] = chebCL(Ny,Ly,2*Ly/5);
    I = eye(Ny);
    
    %GENERATE BLASIUS PROFILE dimensioneless by delta*
    [U,~,U_PP]=blasius_profil(Ly,yvecs);
    %figure(1);
    %plot(U,yvecs,'b.',U_P,yvecs,'r.',U_PP,yvecs,'g.')
    %legend('U','dU','ddU');
    %title('Blasius profile and its first and second derivative')
    
    
    %OS-OPERATOR
    A = (D4-2*D2*kk^2+kk^4*I)-1i*Re*alpha*((diag(U)-I*c)*(D2-kk^2*I)-diag(U_PP));
    A(end,:) = zeros(1,Ny);
    A(end,end) = 1;            % Dirichlet at the wall
    A(end-1,:) = D1(end,:);    % Neuman at the wall
    A(1,:) = zeros(1,Ny);
    A(1,1) = 1+1i*0;                % Dirichlet desomogeneization at the top
    A(2,:) = zeros(1,Ny);   
    A(2,2) = 1;                % Dirichlet desomogeneization v = ? to find 
                               % with Newthon iteration
    
    %Newton iteration on the "v" value on a point in the freestream
    cc1 = 1;      % first iteration value
    chk1 = 10;    % condition that have to be respected:
    
                  % (D2*v+gamma^2*v)_(y1)
                  % ---------------------  =  exp(kk(y2-y1)) (Jacobs 1998)
                  % (D2*v+gamma^2*v)_(y2)
    ite = 1;              
    while abs(chk1)>1e-9
        [V,chk1] = bcA(A,cc1,Ny);      % impose BC and solve the OS system
    
        cc2 = cc1+(rand(1,1)+rand(1,1)*1i)*1e-7;
        [~,chk2] = bcA(A,cc2,Ny);
        
        dF = (chk2-chk1)/(cc2-cc1); % first order derivative evaluation
        cc1 = cc1-chk1/dF;          % Newthon update
            
        fprintf('  %i Newthon residues = %d \n',ite, abs(chk1));
        ite = ite+1;
    end
        
    % find eta value
    S = 1i*Re*alpha*(diag(U)-c*I)-(D2-kk^2*I);
    B = V*0;    %B = -1i*Re*beta*U_P.*V; %REMOVE COUPLING TERM
    
    S(1,:)     = zeros(1,Ny);
    S(1,1)     = 1;            % Dirichlet on the Bottom
    S(end,:)   = zeros(1,Ny);
    S(end,end) = 1;            % Dirichlet on the bottom
    
    B(1)       = 1+1i*0;       % eta = 1 on the top
    B(end)     = 0;            % eta = 0 on the bottom
    
    [L,Q] = lu(S);             
    yy = L\B;
    E = Q\yy;                  % solve SQ system
    
    ydm  = Lreal- 0.2*Lreal;
    ymax = Lreal;
    Ss = V*0;
    for i = 1:Ny
        ys = 1- (yvecs(i)-ydm)/(ymax - ydm);
        if     ys<=0
            Ss(i) = 0;
        elseif ys>0 && ys<1
            Ss(i) = 1/(1+exp(1/(ys-1)+1/ys));
            if Ss(i)<1e-50
               Ss(i)=0;
            end
        elseif ys>=1
            Ss(i) = 1;
        end
    end
    
    
    V  = V.*Ss;
    dV = D1*V.*Ss+V.*(D1*Ss);
    E = E.*Ss;
    
    if strcmp(plotres,'true') == 1
       figure(3)
       subplot(1,2,1)
       plot(real(V),yvecs,'--r',imag(V),yvecs,'--b',abs(V),yvecs,'-g')
       legend('v_r','v_i','|v|');
       title('v')
       subplot(1,2,2)
       plot(real(E),yvecs,'--r',imag(E),yvecs,'--b',abs(E),yvecs,'-g')
       legend('eta_r','eta_i','|eta|');
       title('eta')
    end
    
    
    % (v, eta) -------> (u, v, w)
    Uos = 1i*alpha/kk^2*dV;
    Wos = 1i*beta/kk^2*dV;
    
    Usq = -1i*beta/kk^2*E;
    Wsq = 1i*alpha*E/kk^2;
    
    %interpolate on constant grid
    phi = rand(1,1)*2*pi;
    Uos = (interp1(yvecs,Uos,y,'spline'))*exp(1i*phi);
    Vos = (interp1(yvecs,V,y,'spline'))*exp(1i*phi);
    Wos = (interp1(yvecs,Wos,y,'spline'))*exp(1i*phi);
    
    phi = rand(1,1)*2*pi;
    Usq = (interp1(yvecs,Usq,y,'spline'))*exp(1i*phi);
    Wsq = (interp1(yvecs,Wsq,y,'spline'))*exp(1i*phi);
    
    phi = rand(1,1)*2*pi;
    U = cos(phi)*Uos+sin(phi)*Usq;
    V = cos(phi)*Vos;
    W = cos(phi)*Wos+sin(phi)*Wsq;
    
    % the normalization energy take into account just the domain out to 
    % the BL ----> is sure to have the imposed tu
    y1 = linspace(5,ydm,Ny);    
    U1 = interp1(y,U,y1,'spline');
    V1 = interp1(y,V,y1,'spline');
    W1 = interp1(y,W,y1,'spline');

    % Normalization to have unit energy + gives random phase
    ENERGY = U1.*conj(U1) + V1.*conj(V1) + W1.*conj(W1);
    ENERGIA = 0.5*abs(trapz(y1,ENERGY))/abs(y1(1)-y1(end));
    
    U = U/sqrt(ENERGIA);
    V = V/sqrt(ENERGIA);
    W = W/sqrt(ENERGIA);
    
    %check
    ENERGY = U.*conj(U) + V.*conj(V) + W.*conj(W);
    
    if strcmp(plotres,'true') == 1
       figure(4)
       subplot(1,4,1)
       plot(real(U),y,'--r',imag(U),y,'--b',abs(U),y,'-g')
       xlim([-2 2])
       %legend('u_r','u_i','|u|');
       title('u')
       subplot(1,4,2)
       plot(real(V),y,'--r',imag(V),y,'--b',abs(V),y,'-g')
       xlim([-2 2])
       %legend('v_r','v_i','|v|');
       title('v')
       subplot(1,4,3)
       plot(real(W),y,'--r',imag(W),y,'--b',abs(W),y,'-g')
       xlim([-2 2])
       %legend('w_r','w_i','|w|');
       title('w')
       subplot(1,4,4)
       plot(ENERGY,y,'-m')
       title('energy')
    end

    U=flipud(U);
    V=flipud(V);
    W=flipud(W);
    y=flipud(y);

    
    numnp = length(y);
    name = ['./FST_data/velocity' numero '.dat'];
    fid = fopen(name,'w');
    fprintf(fid,'          %i\n',numnp);
    for ii = 1:numnp
        fprintf(fid,'%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e\n',y(ii),real(U(ii)),imag(U(ii)),real(V(ii)),imag(V(ii)),real(W(ii)),imag(W(ii)));
    end
    fclose(fid);
end
