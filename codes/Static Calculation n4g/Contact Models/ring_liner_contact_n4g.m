function [fcl,mcl,hmin,z0]=ring_liner_contact_n4g(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yr,alr,alp,PR,sigmap,zk,Pk,omega)
%alr=alr-alp;
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
    z0=rbn*ones(nnz(checktwopoints),1);
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
    
    hmin=min(h);
    %hmin=hmin.*(hmin>0)+0.01*sigmap*(hmin<=0);
    xxn=alr2.*(arm-ringface)+xx;
    
    Pcl=PR*Pk*(omega-h/sigmap).^zk.*(h/sigmap<=omega);
    dmcl=Pcl.*xxn;
    
    fcl=-sum((Pcl(1:end-1,:)+Pcl(2:end,:))/2*dx);
    mcl=sum((dmcl(1:end-1,:)+dmcl(2:end,:))/2*dx);

else
    Nbz=100;
    
    fcl=zeros(1,Npe);
    mcl=zeros(1,Npe);
    
    z0=z0_1.*(z0_1<rbn).*(z0_2<rbn)+z0_2.*(z0_2>rbn).*(z0_1>rbn)+rbn*(z0_1>=rbn).*(z0_2<=rbn);
    
    hmin=(a12*(z0-rbn).^2+a11*(z0-rbn)+a10+Db/2+yb-Rr-yr-arm+alr.*z0).*(z0<=rbn)...
        +(a22*(z0-rbn).^2+a21*(z0-rbn)+a20+Db/2+yb-Rr-yr-arm+alr.*z0).*(z0>rbn);
    %hmin=hmin.*(hmin>0)+0.01*sigmap*(hmin<=0);
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
        Pcl=PR*Pk*(omega-h/sigmap).^zk.*(h<(omega*sigmap));
        dmcl=Pcl.*zn;
        fcl(contact)=-sum((Pcl(1:end-1,:)+Pcl(2:end,:))/2).*dz;
        mcl(contact)=sum((dmcl(1:end-1,:)+dmcl(2:end,:))/2).*dz;
    end
end





