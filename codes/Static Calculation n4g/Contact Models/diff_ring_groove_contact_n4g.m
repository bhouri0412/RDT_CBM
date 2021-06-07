function [dfdz,dfdalr,dmdz,dmdalr]=...
    diff_ring_groove_contact_n4g(gcl,dz_g,zg,zr,alr,alp,beta,beta_th,y1,y2,sigmag,z,Pk,omega,isu,thrt,thrb)
%alr=alr-alp;
Npe=length(zr);
Nby=50;
%thr=(isu==1)*thrt+(isu==-1)*thrb;

dfdz=zeros(1,Npe);
dfdalr=zeros(1,Npe);
dmdz=zeros(1,Npe);
dmdalr=zeros(1,Npe);

h0=gcl/2+isu*(zg+dz_g-zr);
al=(alr-beta-beta_th);
al=al.*(abs(al)>1e-8)+1e-8*(abs(al)<1e-8);


minh=(h0-isu*y2*al).*(isu*al>0) ...
    +(h0-isu*y1*al).*(isu*al<0);
contact=minh<(4*sigmag);
%contact region
ncontact=nnz(contact);
if ncontact
    yc=(4*sigmag-h0(contact))./(-isu*al(contact));
    yc=min(yc,y2);
    yc=max(yc,y1);
    ya=y1*(isu*al(contact)<0)+yc.*(isu*al(contact)>0);
    ya=((y2-ya)>1e-8).*ya+(y2-1e-8).*((y2-ya)<1e-8);  %ya cannot excess y2
    yb=y2*(isu*al(contact)>0)+yc.*(isu*al(contact)<0);
    yb=((yb-y1)>1e-8).*yb+(y1+1e-8).*((yb-y1)<1e-8);    %yb cannot smaller than y1
    
    dy=(yb-ya)/Nby;
    dy2=repmat(dy,Nby,1);
    yy=repmat(ya,Nby+1,1)+repmat([0:Nby]',1,ncontact).*repmat(dy,Nby+1,1);
    h0_2=repmat(h0(contact),Nby+1,1);
    al_2=repmat(al(contact),Nby+1,1);
    h=h0_2-isu*al_2.*yy;
    %h=h.*(h>1e-8)+1e-8*(h<1e-8);
    
    dPcgdz=Pk*z*(omega-h/sigmag).^(z-1)*(isu/sigmag).*(h/sigmag<=omega);
    dPcgdalr=Pk*z*(omega-h/sigmag).^(z-1)*(isu/sigmag).*yy.*(h/sigmag<=omega);
    dmcgdz=dPcgdz.*yy;
    dmcgdalr=dPcgdalr.*yy;
    
    dfdz(contact)=sum(-isu*(dPcgdz(1:end-1,:)+dPcgdz(2:end,:))/2.*dy2);
    dfdalr(contact)=sum(-isu*(dPcgdalr(1:end-1,:)+dPcgdalr(2:end,:))/2.*dy2);
    dmdz(contact)=sum(-isu*(dmcgdz(1:end-1,:)+dmcgdz(2:end,:))/2.*dy2);
    dmdalr(contact)=sum(-isu*(dmcgdalr(1:end-1,:)+dmcgdalr(2:end,:))/2.*dy2);
end

end
