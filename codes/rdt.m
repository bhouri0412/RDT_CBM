function rdt()

[Ac,y_c,z_c,auo,huo,aui,hui,ali,hli,alo,hlo,thrt,thrb,ylep,zlep,arm,rbn,yuep,zuep,rb1,rb2,a20,a21,a22,a10,a11,a12,Iyr,Izr,Iyz,Izz,Iyy,Jt,alp,Ip]=rdt_cross_section();

inputs.ring.auo=auo;
inputs.ring.aui=aui;
inputs.ring.alo=alo;
inputs.ring.ali=ali;
inputs.ring.huo=huo;
inputs.ring.hui=hui;
inputs.ring.hlo=hlo;
inputs.ring.hli=hli;
inputs.ring.thrt=thrt;           
inputs.ring.thrb=thrb;           
inputs.ring.rb1=rb1;
inputs.ring.rb2=rb2;
inputs.ring.rbn=rbn;
inputs.ring.a10=a10;
inputs.ring.a11=a11;
inputs.ring.a12=a12;
inputs.ring.a20=a20;
inputs.ring.a21=a21;
inputs.ring.a22=a22;
inputs.ring.arm=arm;
inputs.ring.Ac=Ac;
inputs.ring.Izz=Izz;
inputs.ring.Iyy=Iyy;
inputs.ring.Izr=Izr;
inputs.ring.Ip=Ip;
inputs.ring.Jt=Jt;
inputs.ring.alp=alp;
inputs.ring.y_c=y_c;
inputs.ring.z_c=z_c;
inputs.ring.ylep=ylep;
inputs.ring.zlep=zlep;
inputs.ring.yuep=yuep;
inputs.ring.zuep=zuep;
inputs.ring.Iyr=Iyr;
inputs.ring.Iyz=Iyz;

prompt='Type 1 to compute ovality and pressure distribution for circular shape from free shape, 0 otherwise: ';
Isfromfs=input(prompt);
while Isfromfs
    [rforce,rov]=rdt_from_fs_to_pr_and_ov(inputs);
    %save the results
    %this code computes the radial linear force: 1st column contains angles in degrees and second one force par 
    %circumference length in N/m
    %and the ovality in usual representation: 1st column contains angles in degrees and second one coordintes in microns  
    %given the freeshape in millimeters with angles in degrees
    dlmwrite('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\rforce_fromfs.txt',rforce,'delimiter','\t','precision','%12.8f','newline','pc');
    dlmwrite('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\rov_fromfs.txt',rov,'delimiter','\t','precision','%12.8f','newline','pc');
    prompt='Type 1 to compute ovality and pressure distribution for circular shape from free shape again, 0 otherwise: ';
    Isfromfs=input(prompt);
end

prompt='Type 1 to compute ovality and free shape from pressure distribution for circular shape again, 0 otherwise: ';
Isfrompr=input(prompt);
while Isfrompr
    [fs,rov]=rdt_from_pr_to_fs_and_ov(inputs);
    %save the results
    %this code computes the free shape in millimeters and ovality in microns in usual representation 
    %from radial linear force as input: 
    %1st column should contain angles in degrees and second one force per circumference length in N/m
    dlmwrite('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\fs_frompr.txt',fs,'delimiter','\t','precision','%12.8f','newline','pc');
    dlmwrite('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\rov_frompr.txt',rov,'delimiter','\t','precision','%12.8f','newline','pc');
    prompt='Type 1 to compute ovality and free shape from pressure distribution for circular shape, 0 otherwise: ';
    Isfrompr=input(prompt);
end

prompt='Type 1 to compute free shape and pressure distribution for circular shape from ovality again, 0 otherwise: ';
Isfromov=input(prompt);
while Isfromov
    [fs,rforce]=rdt_from_ov_to_fs_and_pr(inputs);
    %save the results
    %this code computes the radial linear force: 1st column contains angles in degrees and second one force par 
    %circumference length in N/m
    %and the freeshape in usual representation;1st column contains angles in degrees and second one coordintes in millimeters  
    %given the ovality in microns with angles in degrees
    dlmwrite('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\rforce_fromov.txt',rforce,'delimiter','\t','precision','%12.8f','newline','pc');
    dlmwrite('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\fs_fromov.txt',fs,'delimiter','\t','precision','%12.8f','newline','pc');
    prompt='Type 1 to compute free shape and pressure distribution for circular shape from ovality again, 0 otherwise: ';
    Isfromov=input(prompt);
end

prompt='Type 1 to launch the conformability study for one specific ring gap position, 0 otherwise: ';
Isconf1g=input(prompt);
while Isconf1g
    [res_theer,res_hmin,res_zr,res_yr,res_alr,res_z0,res_uoc,res_uic,res_loc,res_lic,...
        res_thee,res_fy,res_fl,res_fz,res_m,res_Mfinal,res_S_uo,res_S_ui,res_S_lo,res_S_li,Nplot,Nr,ybr,Db,Tg,IsBDist,gap]=rdt_conformability(inputs);
    mat_result=[res_theer,res_hmin,res_zr,res_yr,res_alr,res_z0,res_uoc,res_uic,res_loc,res_lic,...
        res_thee,res_fy,res_fl,res_fz,res_m,res_Mfinal,res_S_uo,res_S_ui,res_S_lo,res_S_li];
    %save the results
    dlmwrite('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\res_conf_one_gap.txt',mat_result,'delimiter','\t','precision','%12.8f','newline','pc');   
    prompt='Type 1 to create circular plot of distorted bore and ring with radial force representation, 0 otherwise: ';
    Iscircle=input(prompt);
    while Iscircle
        rdt_conformability_radial_plot(Nplot,Nr,ybr,Db,Tg,res_hmin*1e-6,res_fy,res_fl,IsBDist,gap)
        prompt='Type 1 to create a new circular plot of distorted bore and ring with radial force representation, 0 otherwise: ';
        Iscircle=input(prompt);
    end
    prompt='Type 1 to launch the conformability study again for one specific ring gap position, 0 otherwise: ';
    Isconf1g=input(prompt);
end

prompt='Type 1 to launch the conformability study for different ring gap positions equally spaced between 0 and 359 degrees, 0 otherwise: ';
Isconfdiffg=input(prompt);
while Isconfdiffg
    [res_Tg,res_theer,res_hmin,res_zr,res_yr,res_alr,res_z0,res_uoc,res_uic,res_loc,res_lic,res_ybr,...
        res_thee,res_fy,res_fl,res_fz,res_m,res_Mfinal,res_S_uo,res_S_ui,res_S_lo,res_S_li,util_data]=rdt_conformability_video(inputs);
    %save the results
    
    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_theer.txt',res_theer,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_hmin.txt',res_hmin,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_zr.txt',res_zr,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_yr.txt',res_yr,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_alr.txt',res_alr,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_z0.txt',res_z0,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_uoc.txt',res_uoc,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_uic.txt',res_uic,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_loc.txt',res_loc,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_lic.txt',res_lic,'delimiter','\t','precision','%12.8f','newline','pc');

 
    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_thee.txt',res_thee,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_fy.txt',res_fy,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_fl.txt',res_fl,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_fz.txt',res_fz,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_m.txt',res_m,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_Mfinal.txt',res_Mfinal,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_S_uo.txt',res_S_uo,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_S_lo.txt',res_S_lo,'delimiter','\t','precision','%12.8f','newline','pc');

    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_S_ui.txt',res_S_ui,'delimiter','\t','precision','%12.8f','newline','pc');
 
    dlmwrite('/Users/mohamedazizbhouri/Desktop/RDT_CBM/codes/texts/res_S_li.txt',res_S_li,'delimiter','\t','precision','%12.8f','newline','pc');

    
    prompt='Type 1 to create circular plot of distorted bore and ring with radial force representation, 0 otherwise: ';
    Iscircle=input(prompt);
    while Iscircle
        rdt_conformability_video_animation(res_Tg,res_hmin,res_ybr,util_data,res_fy,res_fl)
        prompt='Type 1 to create a new circular plot of distorted bore and ring with radial force representation, 0 otherwise: ';
        Iscircle=input(prompt);
    end
    prompt='Type 1 to launch the conformability study again for different ring gap positions equally spaced between 0 and 359 degrees, 0 otherwise: ';
    Isconfdiffg=input(prompt);
end

prompt='Type 1 to launch the static twist study, 0 otherwise: ';
Isstattwist=input(prompt);
while Isstattwist
    [res_thee,res_yr,res_hmin,res_alr,res_loc,res_lic,res_fy,res_fz,res_m]=rdt_static_twist(inputs);
    mat_result=[res_thee,res_yr,res_hmin,res_alr,res_loc,res_lic,res_fy,res_fz,res_m];
    %save the results
    dlmwrite('C:\Users\YangLiu\Desktop\Aziz Perso\papers mit\RA\RDT\res_stat_twist.txt',mat_result,'delimiter','\t','precision','%12.8f','newline','pc');   
    prompt='Type 1 to launch the static twist study again, 0 otherwise: ';
    Isstattwist=input(prompt);
end

end