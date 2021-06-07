function [dfdy,dfdalr,dmdy,dmdalr]=diff_ring_liner_contact_vg(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yr,alr,PR,sigmap,zk,Pk,omega)

Npe=length(yr);

%
z0_1=(rbn-(a11+alr)/2/a12);
z0_1=min(z0_1,rb2+rbn);
z0_1=max(z0_1,-rb1+rbn);
z0_2=(rbn-(a21+alr)/2/a22);
z0_2=min(z0_2,rb2+rbn);
z0_2=max(z0_2,-rb1+rbn);
checktwopoints=(z0_1<rbn).*(z0_2>rbn);
if nnz(checktwopoints)
    warning('bad ring profile');
    Nbx=100;
    
    dx=(rb1+rb2)/Nbx;
    x1=-rb1:dx:0;
    xx1=repmat(x1',1,Npe);
    x2=dx:dx:rb2;
    xx2=repmat(x2',1,Npe);
    xx=[xx1; xx2]+rbn;
    
    delybr=repmat(Db/2+yb-Rr-yr-arm,size(xx,1),1);
    alr2=repmat(alr,size(xx,1),1);
    ringface=[(a12*xx1.^2+a11*xx1+a10); (a22*xx2.^2+a21*xx2+a20)];
    
    h=delybr+ringface+alr2.*xx;
    
    xxn=alr2.*(arm-ringface)+xx;
    
    Pcl=PR*Pk*(omega-h/sigmap).^zk.*(h/sigmap<=omega);
    
    dPcldy=PR*Pk/sigmap*zk*(omega-h/sigmap).^(zk-1).*(h/sigmap<=omega);
    dPcldalr=PR*Pk*zk*(omega-h/sigmap).^(zk-1).*(-xx/sigmap).*(h/sigmap<=omega);
    % dPcldalr=PR*Pgt*z*(omega-h/sigmap).^(z-1).*(-xx/sigmap).*(h/sigmap<=omega)*0; %for testing
    
    dmcldy=dPcldy.*xxn;
    dmcldalr=dPcldalr.*xxn+Pcl.*(arm-ringface);
    % dmcldalr=dPcldalr.*xxn+Pcl.*(yc-ringface)*0;
    
    dfdy=-sum((dPcldy(1:end-1,:)+dPcldy(2:end,:))/2*dx);
    dfdalr=-sum((dPcldalr(1:end-1,:)+dPcldalr(2:end,:))/2*dx);
    
    dmdy=sum((dmcldy(1:end-1,:)+dmcldy(2:end,:))/2*dx);
    dmdalr=sum((dmcldalr(1:end-1,:)+dmcldalr(2:end,:))/2*dx);
else   
    Nbz=50;
    
    dfdy=zeros(1,Npe);
    dfdalr=zeros(1,Npe);
    dmdy=zeros(1,Npe);
    dmdalr=zeros(1,Npe);
    
    z0=z0_1.*(z0_1<rbn).*(z0_2<rbn)+z0_2.*(z0_2>rbn).*(z0_1>rbn)+rbn*(z0_1>=rbn).*(z0_2<=rbn);
    
    hmin=(a12*(z0-rbn).^2+a11*(z0-rbn)+a10+Db/2+yb-Rr-yr-arm+alr.*z0).*(z0<=rbn)...
        +(a22*(z0-rbn).^2+a21*(z0-rbn)+a20+Db/2+yb-Rr-yr-arm+alr.*z0).*(z0>rbn);
    
    contact=hmin<(omega*sigmap);
    %contact region
    ncontact=nnz(contact);
    if ncontact
        alrc=alr(contact);
        ybc=yb(contact);
        yrc=yr(contact);
        z0c=z0(contact);
        za1=(-(a11+alrc)-sqrt((a11+alrc).^2-4*a12*(a10+alrc*rbn+Db/2+ybc-Rr-yrc-arm-omega*sigmap)))/2/a12+rbn;
        za2=(-(a21+alrc)-sqrt((a21+alrc).^2-4*a22*(a20+alrc*rbn+Db/2+ybc-Rr-yrc-arm-omega*sigmap)))/2/a22+rbn;
        zb1=(-(a11+alrc)+sqrt((a11+alrc).^2-4*a12*(a10+alrc*rbn+Db/2+ybc-Rr-yrc-arm-omega*sigmap)))/2/a12+rbn;
        zb2=(-(a21+alrc)+sqrt((a21+alrc).^2-4*a22*(a20+alrc*rbn+Db/2+ybc-Rr-yrc-arm-omega*sigmap)))/2/a22+rbn;
        za=za1.*(z0c<=rbn)+(za2.*(za2>rbn)+za1.*(za2<rbn)).*(z0c>rbn);
        zb=zb2.*(z0c>=rbn)+(zb1.*(zb1<rbn)+zb2.*(zb1>rbn)).*(z0c<rbn);
        za=max(za,-rb1+rbn);
        za=min(za,rb2+rbn);
        zb=max(zb,-rb1+rbn);
        zb=min(zb,rb2+rbn);
        dz=(zb-za+1e-12)/Nbz;
        %dz2=repmat(dz,Nbz,1);
        zz=repmat(za,Nbz+1,1)+repmat([0:Nbz]',1,ncontact).*repmat(dz,Nbz+1,1);
        
        delybr=repmat(Db/2+ybc-Rr-yrc-arm,Nbz+1,1);
        alr2=repmat(alrc,Nbz+1,1);
        ringface=(a12*(zz-rbn).^2+a11*(zz-rbn)+a10).*(zz<=rbn)+(a22*(zz-rbn).^2+a21*(zz-rbn)+a20).*(zz>rbn);
        h=delybr+ringface+alr2.*zz;
        zn=alr2.*(arm-ringface)+zz;
        checkcontact=h<(omega*sigmap);
        Pcl=PR*Pk*(omega-h/sigmap).^zk.*checkcontact;
        
        dPcldy=PR*Pk/sigmap*zk*(omega-h/sigmap).^(zk-1).*checkcontact;
        %dPcldalr=PR*Pk*zk*(omega-h/sigmap).^(zk-1).*(-zz/sigmap).*(h<(omega*sigmap));
        dPcldalr=dPcldy.*(-zz);
        dmcldy=dPcldy.*zn;
        dmcldalr=dPcldalr.*zn+Pcl.*(arm-ringface);
        dfdy(contact)=-sum((dPcldy(1:end-1,:)+dPcldy(2:end,:))/2).*dz;
        dfdalr(contact)=-sum((dPcldalr(1:end-1,:)+dPcldalr(2:end,:))/2).*dz;
        dmdy(contact)=sum((dmcldy(1:end-1,:)+dmcldy(2:end,:))/2).*dz;
        dmdalr(contact)=sum((dmcldalr(1:end-1,:)+dmcldalr(2:end,:))/2).*dz;
    end
end