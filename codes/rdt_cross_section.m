function [Ac,y_c,z_c,auo,huo,aui,hui,ali,hli,alo,hlo,thrt,thrb,ylep,zlep,arm,rbn,yuep,zuep,rb1,rb2,a20,a21,a22,a10,a11,a12,Iyr,Izr,Iyz,Izz,Iyy,Jt,alp,Ip]=rdt_cross_section()

prompt='Type 1 to give points coordinate, 0 to define manually the parameters: ';
IsSerieofpoint=input(prompt);
if IsSerieofpoint
    prompt='Number of points defining the cross section, provide the points in clockwise order with OD on the left: ';
    Npoints=input(prompt);
    z=zeros(Npoints,1);
    y=zeros(Npoints,1);
    for i=1:Npoints %ask for points cooridnates in the clockwise order
        prompt = 'Radial cooridnate (mm): ';
        y(i)=input(prompt)*1e-3;
        prompt = 'Axial cooridnate (mm): ';
        z(i)=input(prompt)*1e-3;
    end
    
    y=-y;
    
    prompt='id for upper OD point; ';
    induo=input(prompt);
    prompt='id for upper ID point: ';
    indui=input(prompt);
    prompt='id for lower ID point: ';
    indli=input(prompt);
    prompt='id for lower OD point: ';
    indlo=input(prompt);
    prompt='id for lower end point: ';
    indle=input(prompt);
    prompt='id for minimum point: ';
    indmp=input(prompt);
    prompt='id for upper end point: ';
    indue=input(prompt);
    
    Ac=0;
    for i=1:Npoints-1
        Ac=Ac+y(i)*z(i+1)-y(i+1)*z(i);
    end
    Ac=Ac+y(Npoints)*z(1)-y(1)*z(Npoints);
    Ac=Ac/2;   
    
    y_c=0;
    z_c=0;
    
    for i=1:Npoints-1
        y_c=y_c+(y(i)+y(i+1))*(y(i)*z(i+1)-y(i+1)*z(i));
        z_c=z_c+(z(i)+z(i+1))*(y(i)*z(i+1)-y(i+1)*z(i));
    end
    y_c=y_c+(y(Npoints)+y(1))*(y(Npoints)*z(1)-y(1)*z(Npoints));
    z_c=z_c+(z(Npoints)+z(1))*(y(Npoints)*z(1)-y(1)*z(Npoints));
    y_c=y_c/(6*Ac);
    z_c=z_c/(6*Ac);
        
    auo=y(induo)-y_c;
    huo=z(induo)-z_c; 
    aui=y_c-y(indui);
    hui=z(indui)-z_c;
    ali=y_c-y(indli); 
    hli=z_c-z(indli);
    alo=y(indlo)-y_c; 
    hlo=z_c-z(indlo);  
    ylep=y(indle)-y_c; 
    zlep=z(indle)-z_c; 
    arm=y(indmp)-y_c; 
    rbn=z(indmp)-z_c; 
    yuep=y(indue)-y_c; 
    zuep=z(indue)-z_c;
    thrt=atan((huo-hui)/(auo+aui));
    thrb=atan((hlo-hli)/(alo+ali));
    rb1=rbn-zlep;
    rb2=zuep-rbn; 
    a20=0;
    prompt='Linear cofficient for upper edge shape factor: ';
    a21=input(prompt);
    a22=(arm-yuep-a21*rb2)/rb2^2;
    a10=0;
    prompt='Linear cofficient for lower edge shape factor: ';
    a11=input(prompt);
    a12=(arm-ylep-a11*rb1)/rb1^2;
    
	yc=y-y_c;
    zc=z-z_c;
    
    y_c=-y_c;
   
    Izr=0;
    for i=1:Npoints-1
        Izr=Izr+(yc(i)*zc(i+1)-yc(i+1)*zc(i))*(yc(i)^2+yc(i)*yc(i+1)+yc(i+1)^2);
    end
    Izr=Izr+(yc(Npoints)*zc(1)-yc(1)*zc(Npoints))*(yc(Npoints)^2+yc(Npoints)*yc(1)+yc(1)^2);
    Izr=Izr/12;
    
    Iyr=0;
    for i=1:Npoints-1
        Iyr=Iyr+(yc(i)*zc(i+1)-yc(i+1)*zc(i))*(zc(i)^2+zc(i)*zc(i+1)+zc(i+1)^2);
    end
    Iyr=Iyr+(yc(Npoints)*zc(1)-yc(1)*zc(Npoints))*(zc(Npoints)^2+zc(Npoints)*zc(1)+zc(1)^2);
    Iyr=Iyr/12;
    
    Iyz=0; 
    for i=1:Npoints-1
        Iyz=Iyz+(yc(i)*zc(i+1)-yc(i+1)*zc(i))*(2*yc(i)*zc(i)+yc(i)*zc(i+1)+yc(i+1)*zc(i)+2*yc(i+1)*zc(i+1));
    end
    Iyz=Iyz+(yc(Npoints)*zc(1)-yc(1)*zc(Npoints))*(2*yc(Npoints)*zc(Npoints)+yc(Npoints)*zc(1)+yc(1)*zc(Npoints)+2*yc(1)*zc(1));
    Iyz=Iyz/24;
    
    alp=1/2*atan(2*Iyz/(Iyr-Izr)); 
    
    Izz=Iyr*(sin(alp))^2+Izr*(cos(alp))^2-Iyz*sin(2*alp); 
    Iyy=Iyr*(cos(alp))^2+Izr*(sin(alp))^2+Iyz*sin(2*alp); 
    
    alp=-alp;
    
    Ip=Izz+Iyy;  
    
    b=(12^2*Izz^3/Iyy)^(1/8);
    h=(12*Iyy/b)^(1/3);
    L=max(b,h);
    l=min(b,h);
    Jt=L*l^3*(1/3-0.21*l/L*(1-l^4/12/L^4));
    
    display('center of mass radial coordinate (mm): ');
    display(y_c*1e3);
    display('center of mass axial coordinate (mm): ');
    display(z_c*1e3);
    display('upper OD width (mm): ');
    display(auo*1e3);
    display('upper OD height (mm): ');
    display(huo*1e3);
    display('upper ID width (mm): ');
    display(aui*1e3);
    display('upper ID height (mm): ');
    display(hui*1e3);
    display('lower OD width (mm): ');
    display(alo*1e3);
    display('lower OD height (mm): ');
    display(hlo*1e3);
    display('lower ID width (mm): ');
    display(ali*1e3);
    display('lower ID height (mm): ');
    display(hli*1e3);
    display('upper end point width (mm): ');
    display(yuep*1e3);
    display('upper end point axial location (mm): ');
    display(zuep*1e3);
    display('lower end point width (mm): ');
    display(ylep*1e3);
    display('lower end point axial location (mm): ');
    display(zlep*1e3);
    display('minimum point width (mm): ');
    display(arm*1e3);
    display('minimum point axial location (mm): ');
    display(rbn*1e3);
    display('ring upper flank angle (deg): ');
    display(thrt*180/pi);
    display('ring lower flank angle (deg): ');
    display(thrb*180/pi);
    display('lower edge width (mm): ');
    display(rb1*1e3);
    display('upper edge width (mm): ');
    display(rb2*1e3);
    display('quadratic coefficient in upper edge shape factor (1/m): ');
    display(a22);
    display('quadratic coefficient in upper edge shape factor (1/m): ');
    display(a12);
    display('cross section area (mm^2): ');
    display(Ac*1e6);
    display('moment of inertia in plane Iz (mm^4): ');
    display(Izr*1e12);
    display('moment of inertia out of plane Iy (mm^4): ');
    display(Iyr*1e12);
    display('product of intertia Iyz (mm^4): ');
    display(Iyz*1e12);
    display('principal moment of inertial in plane Izp (mm^4): ');
    display(Izz*1e12);
    display('principal moment of inertial out of plane Iyp (mm^4): ');
    display(Iyy*1e12);
    display('principal angle (deg): ');
    display(alp*180/pi);
    display('polar moment of inertia (mm^4): ');
    display(Ip*1e12);
    display('torsional factor (mm^4): ');
    display(Jt*1e12);
    
    y=-y;
    
    figure(1);
    for i=1:Npoints-1
        if y(i)==y(i+1)
            hold on;
            m=min(z(i),z(i+1));
            M=max(z(i),z(i+1));
            ind=(M-m)/1000;
            vect_z=m:ind:M;
            vect_y=y(i)*ones(1,length(vect_z));
            plot(vect_y,vect_z,'b');
        else
            hold on;
            coeff_a=(z(i+1)-z(i))/(y(i+1)-y(i));
            coeff_b=z(i+1)-coeff_a*y(i+1);
            m=min(y(i),y(i+1));
            M=max(y(i),y(i+1));
            ind=(M-m)/1000;
            vect_y=m:ind:M;
            vect_z=coeff_a*vect_y+coeff_b;
            plot(vect_y,vect_z,'b');
        end
    end
    if y(Npoints)==y(1)
        hold on;
        m=min(z(Npoints),z(1));
        M=max(z(Npoints),z(1));
        ind=(M-m)/1000;
        vect_z=m:ind:M;
        vect_y=y(Npoints)*ones(1,length(vect_z));
        plot(vect_y,vect_z,'b');
    else
        hold on;
        coeff_a=(z(1)-z(Npoints))/(y(1)-y(Npoints));
        coeff_b=z(1)-coeff_a*y(1);
        m=min(y(Npoints),y(1));
        M=max(y(Npoints),y(1));
        ind=(M-m)/1000;
        vect_y=m:ind:M;
        vect_z=coeff_a*vect_y+coeff_b;
        plot(vect_y,vect_z,'b');
    end
    
else   
    
    prompt='Upper OD width (mm): ';
    auo=input(prompt)*1e-3;
    prompt='Upper OD height (mm): ';
    huo=input(prompt)*1e-3;
    prompt='Upper ID width (mm): ';
    aui=input(prompt)*1e-3;
    prompt='Upper ID height (mm): ';
    hui=input(prompt)*1e-3;
    prompt='Lower ID width (mm): ';
    ali=input(prompt)*1e-3;
    prompt='Lower ID height (mm): ';
    hli=input(prompt)*1e-3;
    prompt='Lower OD width (mm): ';
    alo=input(prompt)*1e-3;
    prompt='Lower OD height (mm): ';
    hlo=input(prompt)*1e-3;
    prompt='Lower end point width (mm): ';
    ylep=input(prompt)*1e-3;
    prompt='Lower end point axial location (mm): ';
    zlep=input(prompt)*1e-3;
    prompt='Minimum point width (mm): ';
    arm=input(prompt)*1e-3;
    prompt='Minimum point axial location (mm): ';
    rbn=input(prompt)*1e-3;
    prompt='Upper end point width (mm): ';
    yuep=input(prompt)*1e-3;
    prompt='Upper end point axial location (mm): ';
    zuep=input(prompt)*1e-3;
    prompt='Linear cofficient for upper edge shape factor: ';
    a21=input(prompt);
    prompt='Linear cofficient for lower edge shape factor: ';
    a11=input(prompt);
    
    y_c=0;
    z_c=0;
    
    induo=1;
    indui=2;
    indli=3;
    indlo=4;
    indle=5;
    indmp=6;
    indue=7;
    
    Npoints=7;
    
    z=zeros(Npoints,1);
    y=zeros(Npoints,1);
    
    y(induo)=auo; 
    z(induo)=huo; 
    y(indui)=-aui; 
    z(indui)=hui; 
    y(indli)=-ali;   
    z(indli)=-hli; 
    y(indlo)=alo;  
    z(indlo)=-hlo;  
    y(indle)=ylep; 
    z(indle)=zlep; 
    y(indmp)=arm; 
    z(indmp)=rbn;
    y(indue)=yuep; 
    z(indue)=zuep;
    
    thrt=atan((huo-hui)/(auo+aui)); 
    thrb=atan((hlo-hli)/(alo+ali)); 
    rb1=rbn-zlep; 
    rb2=zuep-rbn; 
    
    a20=0;
    a22=(arm-yuep-a21*rb2)/rb2^2;
    a10=0;
    a12=(arm-ylep-a11*rb1)/rb1^2;
    
    Ac=0;
    for i=1:Npoints-1
        Ac=Ac+y(i)*z(i+1)-y(i+1)*z(i);
    end
    Ac=Ac+y(Npoints)*z(1)-y(1)*z(Npoints);
    Ac=Ac/2; 
    
    yc=y-y_c;
    zc=z-z_c;
   
    Izr=0; 
    for i=1:Npoints-1
        Izr=Izr+(yc(i)*zc(i+1)-yc(i+1)*zc(i))*(yc(i)^2+yc(i)*yc(i+1)+yc(i+1)^2);
    end
    Izr=Izr+(yc(Npoints)*zc(1)-yc(1)*zc(Npoints))*(yc(Npoints)^2+yc(Npoints)*yc(1)+yc(1)^2);
    Izr=Izr/12;
    
    Iyr=0;
    for i=1:Npoints-1
        Iyr=Iyr+(yc(i)*zc(i+1)-yc(i+1)*zc(i))*(zc(i)^2+zc(i)*zc(i+1)+zc(i+1)^2);
    end
    Iyr=Iyr+(yc(Npoints)*zc(1)-yc(1)*zc(Npoints))*(zc(Npoints)^2+zc(Npoints)*zc(1)+zc(1)^2);
    Iyr=Iyr/12;
    
    Iyz=0; 
    for i=1:Npoints-1
        Iyz=Iyz+(yc(i)*zc(i+1)-yc(i+1)*zc(i))*(2*yc(i)*zc(i)+yc(i)*zc(i+1)+yc(i+1)*zc(i)+2*yc(i+1)*zc(i+1));
    end
    Iyz=Iyz+(yc(Npoints)*zc(1)-yc(1)*zc(Npoints))*(2*yc(Npoints)*zc(Npoints)+yc(Npoints)*zc(1)+yc(1)*zc(Npoints)+2*yc(1)*zc(1));
    Iyz=Iyz/24;
    
    alp=1/2*atan(2*Iyz/(Iyr-Izr));
    
    Izz=Iyr*(sin(alp))^2+Izr*(cos(alp))^2-Iyz*sin(2*alp); 
    Iyy=Iyr*(cos(alp))^2+Izr*(sin(alp))^2+Iyz*sin(2*alp); 
    
    alp=-alp;
    
    Ip=Izz+Iyy;   
    
    b=(12^2*Izz^3/Iyy)^(1/8);
    h=(12*Iyy/b)^(1/3);
    L=max(b,h);
    l=min(b,h);
    Jt=L*l^3*(1/3-0.21*l/L*(1-l^4/12/L^4));
    
    display('ring upper flank angle (deg): ');
    display(thrt*180/pi);
    display('ring lower flank angle (deg): ');
    display(thrb*180/pi);
    display('lower edge width (mm): ');
    display(rb1*1e3);
    display('upper edge width (mm): ');
    display(rb2*1e3);
    display('quadratic coefficient in upper edge shape factor (1/m): ');
    display(a22);
    display('quadratic coefficient in lower edge shape factor (1/m): ');
    display(a12);
    display('cross section area (mm^2): ');
    display(Ac*1e6);
    display('moment of inertia in plane Iz (mm^4): ');
    display(Izr*1e12);
    display('moment of inertia out of plane Iy (mm^4): ');
    display(Iyr*1e12);
    display('product of intertia Iyz (mm^4): ');
    display(Iyz*1e12);
    display('principal moment of inertial in plane Izp (mm^4): ');
    display(Izz*1e12);
    display('principal moment of inertial out of plane Iyp (mm^4): ');
    display(Iyy*1e12);
    display('principal angle (deg): ');
    display(alp*180/pi);
    display('polar moment of inertia (mm^4): ');
    display(Ip*1e12);
    display('torsional factor Jt (mm^4): ');
    display(Jt*1e12);
    
    y=-y;
    
    figure(1);
    for i=1:Npoints-1
        if y(i)==y(i+1)
            hold on;
            m=min(z(i),z(i+1));
            M=max(z(i),z(i+1));
            ind=(M-m)/1000;
            vect_z=m:ind:M;
            vect_y=y(i)*ones(1,length(vect_z));
            plot(vect_y,vect_z,'b');
        else
            hold on;
            coeff_a=(z(i+1)-z(i))/(y(i+1)-y(i));
            coeff_b=z(i+1)-coeff_a*y(i+1);
            m=min(y(i),y(i+1));
            M=max(y(i),y(i+1));
            ind=(M-m)/1000;
            vect_y=m:ind:M;
            vect_z=coeff_a*vect_y+coeff_b;
            plot(vect_y,vect_z,'b');
        end
    end
    if y(Npoints)==y(1)
        hold on;
        m=min(z(Npoints),z(1));
        M=max(z(Npoints),z(1));
        ind=(M-m)/1000;
        vect_z=m:ind:M;
        vect_y=y(Npoints)*ones(1,length(vect_z));
        plot(vect_y,vect_z,'b');
    else
        hold on;
        coeff_a=(z(1)-z(Npoints))/(y(1)-y(Npoints));
        coeff_b=z(1)-coeff_a*y(1);
        m=min(y(Npoints),y(1));
        M=max(y(Npoints),y(1));
        ind=(M-m)/1000;
        vect_y=m:ind:M;
        vect_z=coeff_a*vect_y+coeff_b;
        plot(vect_y,vect_z,'b');
    end
    
end


end