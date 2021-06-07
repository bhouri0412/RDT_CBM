function [hmin,hot,fhl,ffl,mhl,z0]=ring_liner_hydro_vg(arm,rb1,rb2,rbn,a10,a11,a12,a20,a21,a22,Db,yb,Rr,yr,alr,hOCR,KOCR,Kp,Kf,Ph_OCR,ap,F0,b,c,d,e,rw,muUo,templ,isfuele,Vp,oil,sigmap,oiltreshhold,viscosityfactor)

%Nbz=50;
fhl=zeros(1,length(yr));
ffl=zeros(1,length(yr));
mhl=zeros(1,length(yr));
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
    z0=rbn*ones(nnz(checktwopoints),1);
    hmin=(a12*(z0-rbn).^2+a11*(z0-rbn)+a10+Db/2+yb-Rr-yr-arm+alr.*z0);
    %hmin=hmin.*(hmin>0)+0.01*sigmap*(hmin<=0);
else
    z0=z0_1.*(z0_1<rbn).*(z0_2<rbn)+z0_2.*(z0_2>rbn).*(z0_1>rbn)+rbn*(z0_1>=rbn).*(z0_2<=rbn);
    
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
    Pom=(hOCRnbg/sigmap).^(Kpnbg).*munbg*abs(Vp)/muUo*ap*Ph_OCR;
    Fom=F0*munbg*Vp/sigmap.*(hOCRnbg/sigmap).^Kfnbg;
    a=KOCR+Kpnbg;
    hSigp=hminnbg>sigmap;
    h_hd=max(sigmap,hminnbg);
    %haziz=hminnbg<hOCR;
    
    Phydro=(Pom.*(h_hd/sigmap).^(-a)).*(hSigp)+Pom.*(1+a.*(1-hminnbg/sigmap)).*(~hSigp);
    %Phydro=Phydro.*(haziz)+0.*(~haziz);
    shear=(F0*munbg*Vp./h_hd.*(hOCRnbg./h_hd).^Kfnbg).*hSigp+Fom.*(1+(1+Kfnbg).*(1-hminnbg/sigmap)).*(~hSigp);  
    %shear=shear.*(haziz)+0.*(~haziz);
    fhl(checknbg)=-Phydro*rw;
    ffl(checknbg)=-shear*rw;
    mhl(checknbg)=Phydro*rw.*z0(checknbg)-shear*rw*arm;
    hot(checknbg)=hmin(checknbg)/2;
end

if nbg
    hminbg=hmin(checkbg);
    hminbg=hminbg.*(hminbg>(0.01*sigmap))+0.01*sigmap*(hminbg<=(0.01*sigmap));
    mubg=mu(checkbg);
    
    an=a12*rb1^2./hminbg;
    fhydrofl=2*c*6*mubg*abs(Vp)./hminbg/a12.*((1+tanh((log10(an)+d)/e))/2).^b;
    th1=atan(-sqrt(an));
    th2=fullyfloodedtrailingedge(th1);
    z2=tan(th2).*sqrt(hminbg/a12);
    h2=a12*z2.^2+a11*z2+a10+hminbg;
    %dz=(z2+rb1)/Nbz;
    %dz2=repmat(dz,Nbz,1);
    %zz=repmat(-rb1,Nbz+1,nbg)+repmat([0:Nbz]',1,nbg).*repmat(dz,Nbz+1,1);
    %h=a12*zz.^2+a11*zz+a10+repmat(hminbg,Nbz+1,1);
    %dshear=repmat(mubg,Nbz+1,1)*Vp./h;
    %fshear=sum((dshear(1:end-1,:)+dshear(2:end,:))/2).*dz;
    fshear=mubg*Vp./sqrt(a12*hminbg).*(th2-th1);
    fhl(checkbg)=-fhydrofl;
    ffl(checkbg)=-fshear;
    mhl(checkbg)=fhydrofl.*z0(checkbg)-fshear*arm;
    hot(checkbg)=h2/2;
end





