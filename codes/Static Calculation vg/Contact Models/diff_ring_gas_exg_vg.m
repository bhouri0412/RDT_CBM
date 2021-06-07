function [dfrgasdalr,dmgasdalr]=diff_ring_gas_exg_vg(Pu,Pd,rbn,rb1,rb2,a11,a12,a21,a22,alr)

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
    dz0da=0;
else
    z0=z0_1.*(z0_1<rbn).*(z0_2<rbn)+z0_2.*(z0_2>rbn).*(z0_1>rbn)+rbn*(z0_1>=rbn).*(z0_2<=rbn);
    dz0da=-1/2/a12.*(z0_1<rbn).*(z0_2<rbn).*(z0_1>(-rb1+rbn))-1/2/a22.*(z0_2>rbn).*(z0_1>rbn).*(z0_2<(rb2+rbn));
end

dfrgasdalr=Pu.*dz0da-Pd.*dz0da;

dmgasdalr=Pu.*(-z0.*dz0da)+Pd.*z0.*dz0da;

end