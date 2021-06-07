function [dfhldy,dffldy,dmhldy,dfhldalr,dffldalr,dmhldalr]=diff_ring_liner_hydro_vg(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yr,alr,hOCR,KOCR,Kp,Kf,Ph_OCR,ap,F0,b,c,d,e,rw,muUo,templ,isfuele,Vp,oil,sigmap,oiltreshhold,viscosityfactor)

%Nbz=50;
dfhldy=zeros(1,length(yr));
dfhldalr=zeros(1,length(yr));
dffldy=zeros(1,length(yr));
dffldalr=zeros(1,length(yr));
dmhldy=zeros(1,length(yr));
dmhldalr=zeros(1,length(yr));
%Oil data
zk=oil.zk;
temp1=oil.temp1;
temp2=oil.temp2;
rho_oil=oil.rho_oil;
hlratio=oil.hlratio;
bta1=oil.bta1;
bta2=oil.bta2;
zm0=oil.zm0;

z0_1=(rbn-(a11+alr)/2/a12);
z0_1=min(z0_1,rb2+rbn);
z0_1=max(z0_1,-rb1+rbn);
z0_2=(rbn-(a21+alr)/2/a22);
z0_2=min(z0_2,rb2+rbn);
z0_2=max(z0_2,-rb1+rbn);
checktwopoints=(z0_1<rbn).*(z0_2>rbn);
if nnz(checktwopoints)
    warning('bad ring profile');
    z0=rbn;
    hmin=(a12*(z0-rbn).^2+a11*(z0-rbn)+a10+Db/2+yb-Rr-yr-arm+alr.*z0);
    %hmin=hmin.*(hmin>0)+0.01*sigmap*(hmin<=0);
    dz0da=0;
else
    z0=z0_1.*(z0_1<rbn).*(z0_2<rbn)+z0_2.*(z0_2>rbn).*(z0_1>rbn)+rbn*(z0_1>=rbn).*(z0_2<=rbn);
    dz0da=-1/2/a12.*(z0_1<rbn).*(z0_2<rbn).*(z0_1>(-rb1+rbn))-1/2/a22.*(z0_2>rbn).*(z0_1>rbn).*(z0_2<(rb2+rbn));
    
    hmin=(a12*(z0-rbn).^2+a11*(z0-rbn)+a10+Db/2+yb-Rr-yr-arm+alr.*z0).*(z0<=rbn)...
        +(a22*(z0-rbn).^2+a21*(z0-rbn)+a20+Db/2+yb-Rr-yr-arm+alr.*z0).*(z0>rbn);
    %hmin=hmin.*(hmin>0)+0.01*sigmap*(hmin<=0);
end

%viscosity calculation
mu_0=rho_oil*zk*exp(temp1./(temp2+templ))*1e-6;
beta_mu=10.^(bta1+bta2*templ);
gama_mu=abs(Vp./hmin);
mu=mu_0.*(1+hlratio*(gama_mu./beta_mu))./(1+gama_mu./beta_mu);
mu(isfuele)=mu(isfuele)/viscosityfactor;

%checkbg=hOCR>(a12*rb1^2-rb1*a11+a10+hmin);
checkbg=hOCR>(oiltreshhold*sigmap);
checknbg=~checkbg;
nnbg=nnz(checknbg);
nbg=nnz(checkbg);

if nnbg
    hOCRnbg=hOCR(checknbg);
    hminnbg=hmin(checknbg);
    Kpnbg=Kp(checknbg);
    Kfnbg=Kf(checknbg);
    munbg=mu(checknbg);
    z0nbg=z0(checknbg);
    alrnbg=alr(checknbg);
    dz0danbg=dz0da(checknbg);
    a=KOCR+Kpnbg;
    Pom=(hOCRnbg/sigmap).^(Kpnbg).*munbg*abs(Vp)/muUo*ap*Ph_OCR;
    Fom=F0*munbg*Vp/sigmap.*(hOCRnbg/sigmap).^Kfnbg;
    hSigp=hminnbg>sigmap;
    h_hd=max(sigmap,hminnbg);
    
    Phydro=(Pom.*(h_hd/sigmap).^(-a)).*(hSigp)+Pom.*(1+a.*(1-hminnbg/sigmap)).*(~hSigp);
    
    dPhydrody=Pom.*(h_hd/sigmap).^(-a-1).*a/sigmap.*(hSigp)...
        +Pom.*a/sigmap.*(~hSigp);
    dPhydroda=-dPhydrody.*(z0nbg+alrnbg.*dz0danbg);
    
    dfsheardy=F0*(munbg*Vp)./(h_hd.^2).*((hOCRnbg./h_hd).^Kfnbg).*(1+Kfnbg).*(hSigp)...
        +Fom.*(1+Kfnbg)/sigmap.*(~hSigp);
    dfshearda=-dfsheardy.*(z0nbg+alrnbg.*dz0danbg);
    
    dfhldy(checknbg)=-dPhydrody*rw;
    dfhldalr(checknbg)=-dPhydroda*rw;
    dffldy(checknbg)=-dfsheardy*rw;
    dffldalr(checknbg)=-dfshearda*rw;
    dmhldy(checknbg)=rw*(dPhydrody.*z0nbg-dfsheardy*arm);
    dmhldalr(checknbg)=rw*(dPhydroda.*z0nbg+Phydro.*dz0danbg-dfshearda*arm);
end

if nbg
    hminbg=hmin(checkbg);
    hminbg=hminbg.*(hminbg>(0.01*sigmap))+0.01*sigmap*(hminbg<=(0.01*sigmap));
    mubg=mu(checkbg);
    z0bg=z0(checkbg);
    alrbg=alr(checkbg);
    dz0dabg=dz0da(checkbg);
    
    an=a12*rb1^2./hminbg;
    fhydrofl=2*c*6*mubg*abs(Vp)./hminbg/a12.*((1+tanh((log10(an)+d)/e))/2).^b;
    dfhydrofldh=-6*mubg*abs(Vp)./(a12*hminbg.^2)*2*c.*((1+tanh((log10(an)+d)/e))/2).^b...
        -6*mubg*abs(Vp)./(a12*hminbg)*2*c*b.*((1+tanh((log10(an)+d)/e))/2).^(b-1)*(1/2).*sech((log10(an)+d)/e).^2/e./(hminbg*log(10));
    
    th1=atan(-sqrt(an));
    th2=fullyfloodedtrailingedge(th1);
%     z2=tan(th2).*sqrt(hminbg/a12);
%     h2=a12*z2.^2+a11*z2+a10+hminbg;
%     dz=(z2+rb1)/Nbz;
%     dz2=repmat(dz,Nbz,1);
%     zz=repmat(-rb1,Nbz+1,nbg)+repmat([0:Nbz]',1,nbg).*repmat(dz,Nbz+1,1);
%     h=a12*zz.^2+a11*zz+a10+repmat(hminbg,Nbz+1,1);
    dth2dth1=(12+16*cos(2*th1)+4*cos(4*th1)-16*cos(th2).^2-16*cos(th2).^2.*cos(2*th1))...
        ./(12+16*cos(2*th2)+4*cos(4*th2)-16*cos(th2).^2-16*cos(th2).^2.*cos(2*th2)+8*sin(2*th2).*(2*th2-2*th1+sin(2*th2)-sin(2*th1)));
    
    %dsheardh=-repmat(mubg,Nbz+1,1)*Vp./(h.^2);
    %dz2dh=sqrt(hminbg/a12).*sec(th2).^2.*dth2dth1.*(sqrt(a12*rb1^2)/2*hminbg.^(-3/2))./(1+an)+1/2./sqrt(a12*hminbg).*tan(th2);
    %dfsheardh=sum((dsheardh(1:end-1,:)+dsheardh(2:end,:))/2).*dz+mubg*Vp./h2.*dz2dh;
    dfsheardh=-mubg*Vp/2/sqrt(a12).*(hminbg.^(-3/2)).*(th2-th1)...
        +mubg*Vp./sqrt(a12*hminbg).*(dth2dth1.*sqrt(a12*rb1^2)/2.*hminbg.^(-3/2)./(1+an)-sqrt(a12*rb1^2)/2*hminbg.^(-3/2)./(1+an));
    
    dfhydrofldy=dfhydrofldh*(-1);
    dfhydroflda=dfhydrofldh.*(z0bg+alrbg.*dz0dabg);
    dfshearfldy=dfsheardh*(-1);
    dfshearflda=dfsheardh.*(z0bg+alrbg.*dz0dabg);
    
    dfhldy(checkbg)=-dfhydrofldy;
    dfhldalr(checkbg)=-dfhydroflda;
    dffldy(checkbg)=-dfshearfldy;
    dffldalr(checkbg)=-dfshearflda;
    dmhldy(checkbg)=dfhydrofldy.*z0bg-dfshearfldy*arm;
    dmhldalr(checkbg)=dfhydroflda.*z0bg+fhydrofl.*dz0dabg-dfshearflda*arm;
    
end