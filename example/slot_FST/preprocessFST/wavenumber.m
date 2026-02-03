function [count]=wavenumber(numk,kkini,kkfin,plotres)
%numk = Number of wavenumber k = mod (alpha, beta, gamma)

%kkini   Minimal frequency kmin
%kkfin   Maximal frecuency k max
dkk = (kkfin-kkini)/(numk-1);

fprintf('kmin %d, kmax %d, deltak %d \n',kkini,kkfin,dkk)

%     ============        DODECAEDRON DEFINITION        ================

      phi = (1.+sqrt(5.))/2.;
      radius0 = sqrt(3.);

      vec_c(1,1)=1.;   % Omega
      vec_c(2,1)=1.;   % Gamma
      vec_c(3,1)=1.;   % Beta

      vec_c(1,2)=1.;
      vec_c(2,2)=1.;
      vec_c(3,2)=-1.;

      vec_c(1,3)=1.;
      vec_c(2,3)=-1.;
      vec_c(3,3)=1.;

      vec_c(1,4)=1.;
      vec_c(2,4)=-1.;
      vec_c(3,4)=-1.;

      vec_c(1,5)=-1.;
      vec_c(2,5)=1.;
      vec_c(3,5)=1.;

      vec_c(1,6)=-1.;
      vec_c(2,6)=1.;
      vec_c(3,6)=-1.;

      vec_c(1,7)=-1.;
      vec_c(2,7)=-1.;
      vec_c(3,7)=1.;

      vec_c(1,8)=-1.;
      vec_c(2,8)=-1.;
      vec_c(3,8)=-1.;


      vec_c(1,9)=0.;
      vec_c(2,9)=1./phi;
      vec_c(3,9)=phi;

      vec_c(1,10)=0.;
      vec_c(2,10)=-1./phi;
      vec_c(3,10)=phi;

      vec_c(1,11)=0.;
      vec_c(2,11)=1./phi;
      vec_c(3,11)=-phi;

      vec_c(1,12)=0.;
      vec_c(2,12)=-1./phi;
      vec_c(3,12)=-phi;

      vec_c(1,13)=phi;
      vec_c(2,13)=0.;
      vec_c(3,13)=1./phi;

      vec_c(1,14)=phi;
      vec_c(2,14)=0.;
      vec_c(3,14)=-1./phi;

      vec_c(1,15)=-phi;
      vec_c(2,15)=0.;
      vec_c(3,15)=1./phi;


      vec_c(1,16)=-phi;
      vec_c(2,16)=0.;
      vec_c(3,16)=-1./phi;

      vec_c(1,17)=1./phi;
      vec_c(2,17)=phi;
      vec_c(3,17)=-0.;

      vec_c(1,18)=1./phi;
      vec_c(2,18)=-phi;
      vec_c(3,18)=0.;

      vec_c(1,19)=-1./phi;
      vec_c(2,19)=-phi;
      vec_c(3,19)=0.;
      
      vec_c(1,20)=-1./phi;
      vec_c(2,20)=phi;
      vec_c(3,20)=0.;
%=================================================================
if strcmp(plotres,'true') == 1
   figure(1)
   PLOTTARE(vec_c,'-b')
end

kk = linspace(kkini,kkfin,numk);
nfile = 0;

mkdir FST_data

for i = 1:numk    
    fprintf('wavenumber sphere discretization = %.3f \n',kk(i))
    ok = 0;
    while ok == 0
        random = rand(1,2);
    
        %from cartesian to spherical --> reshape rho --> from spherical to cartesian 
        [vec_p(1,:),vec_p(2,:),vec_p(3,:)] = cart2sph(vec_c(1,:),vec_c(2,:),vec_c(3,:));
        vec_cn = vec_p;
    
        vec_pn(3,:) = vec_p(3,:)/radius0*kk(i);
        vec_pn(1,:) = vec_p(1,:);
        vec_pn(2,:) = vec_p(2,:);
        [vec_cn(1,:),vec_cn(2,:),vec_cn(3,:)] = sph2cart(vec_pn(1,:),vec_pn(2,:),vec_pn(3,:));
        
        %from cartesian to polar--> z dir rotation --> from  polar to cartesian 
        [vec_pn(1,:),vec_pn(2,:)] = cart2pol(vec_cn(1,:),vec_cn(2,:));
        vec_pn(1,:) = vec_pn(1,:)+2*pi*random(1);
        [vec_cn(1,:),vec_cn(2,:)]= pol2cart(vec_pn(1,:),vec_pn(2,:));
    
        %from cartesian to polar--> x dir rotation --> from  polar to cartesian 
        [vec_pn(2,:),vec_pn(3,:)] = cart2pol(vec_cn(2,:),vec_cn(3,:));
        vec_pn(2,:) = vec_pn(2,:)+pi*random(2);
        [vec_cn(2,:),vec_cn(3,:)]= pol2cart(vec_pn(2,:),vec_pn(3,:));
        
        if strcmp(plotres,'true') == 1
           figure(2)
           hold off
           PLOTTARE(vec_cn,'-b')
        end
        mArrow3([0 0 0],[kkfin*(1+0.5) 0 0],'color','k','stemWidth',0.04);
        mArrow3([0 0 0],[0 kkfin*(1+0.5) 0],'color','k','stemWidth',0.04);
        mArrow3([0 0 0],[0 0 kkfin*(1+0.5)],'color','k','stemWidth',0.04);
        set(gca,'fontsize',16)
        view([135 30])
        ax= 3.5;
        axis([-ax ax -ax ax -ax ax])
        titolo = ['k = ' num2str(kk(i))];
        title (titolo)
        text(kkfin*(1+0.8),0,0,['   ' ...
        '$\omega$'],'HorizontalAlignment','left','FontSize',20,'interpreter','Latex');
        text(0,kkfin*(1+0.5),0,['   ' ...
        '$\gamma$'],'HorizontalAlignment','left','FontSize',20,'interpreter','Latex');
        text(0,0,kkfin*(1+0.5),['   ' ...
        '$\beta$'],'HorizontalAlignment','left','FontSize',20,'interpreter','Latex');
        
        % takes the symmetic dodechaedron respect to gamma
        for j = 1:length(vec_cn(1,:))
            vec_cn(1,length(vec_c(1,:))+j)= vec_cn(1,j);
            vec_cn(2,length(vec_c(1,:))+j)= -vec_cn(2,j);
            vec_cn(3,length(vec_c(1,:))+j)= vec_cn(3,j);
        end
    
%         if strcmp(plotres,'true') == 1
%            PLOTTARE(vec_cn(:,length(vec_c)+1:end),'-g')
%            titolo = ['wavenumber = ' num2str(kk)];
%            title (titolo)
%         end
    
        count = 0;
        for j = 1:length(vec_cn)
            if vec_cn(1,j)>0 && vec_cn(2,j)>0  %take positive omega and gamma
                 STORE(:,1+count) = vec_cn(:,j);
                 count = count+1; 
            end
        end
        
        if count == length(vec_cn)/4   %take positive omega and gamma
        ok = 1;
            fprintf('---> rotation = %5.3f, elevation= %5.3f \n',180*2*random(1),180*random(2))
            for m = 1:count
                nfile = nfile +1;
                if nfile<10
                    numero = [num2str(0) num2str(0) num2str(nfile)];
                elseif nfile<100
                    numero = [num2str(0) num2str(nfile)];
                else
                    numero = num2str(nfile);
                end
                ola = sqrt(STORE(1,m)^2+STORE(2,m)^2+STORE(3,m)^2);
                fprintf('---> omega = %5.3f, gamma = %5.3f, beta = %5.3f, check = %5.3f \n',STORE(1,m),STORE(2,m),STORE(3,m),ola)
                name = ['./FST_data/wavenumber' numero '.dat'];
                fid = fopen(name,'w');
                fprintf(fid,'%16f\n',STORE(1,m));
                fprintf(fid,'%16f\n',STORE(2,m));
                fprintf(fid,'%16f\n',STORE(3,m));
                fclose(fid);
            end
        end  
    end    
    pause(.2)
end
 
      
      
      
      
