%parameters
pi=3.1415926;
h=6.626e-34/(2*pi);
e=-1.6e-19;
gama=28e9*2*pi;
Hex=-40;  % the antiferromagnetic
Hkxy=-1.43; % the xy direction
Hkz=-4.02; % the z direction
dt=2e-15;
miu0=1;
alpha=0.005;
n=50000;
g11=1;
g22=1;
g12=1;
g21=1;
for fain=1:1:1
fai=fain*5*3.14/180;
for fn=1:1:50
f=450+fn*2;
w=f*2*pi*1e9;
 %m0
m1(1,1)=-0.964723;
m1(1,2)=0.26326;
m1(1,3)=0;
m1(1,4)=1;
m2(1,1)=0.95555;
m2(1,2)=-0.294803;
m2(1,3)=0;
m2(1,4)=1; 
%H0
H0x=0.0001;
H0y=0.0001;
H0z=0.0001;
%Hk
Hk1(1,1)=-Hkxy*m1(1,1);  
Hk1(1,2)=Hkxy*m1(1,2);  
Hk1(1,3)=Hkz*m1(1,3);  
Hk2(1,1)=Hkxy*m2(1,1);  
Hk2(1,2)=-Hkxy*m2(1,2);  
Hk2(1,3)=Hkz*m2(1,3);
%Hex
Hex1(1,1)=Hex*m2(1,1);  
Hex1(1,2)=Hex*m2(1,2);  
Hex1(1,3)=Hex*m2(1,3);
Hex2(1,1)=Hex*m1(1,1);  
Hex2(1,2)=Hex*m1(1,2);  
Hex2(1,3)=Hex*m1(1,3);
%Heff
Heff1(1,1)=H0x+Hex1(1,1)+Hk1(1,1);
Heff1(1,2)=H0y+Hex1(1,2)+Hk1(1,2);
Heff1(1,3)=H0z+Hex1(1,3)+Hk1(1,3);
Heff2(1,1)=H0x+Hex2(1,1)+Hk2(1,1);
Heff2(1,2)=H0y+Hex2(1,2)+Hk2(1,2);
Heff2(1,3)=H0z+Hex2(1,3)+Hk2(1,3);
t(1,1)=0;
for i=1:1:n    
    t(i+1,1)=i*dt;
%LLG
    rt1(i+1,1)=-1*gama*miu0*(m1(i,2)*Heff1(i,3)-m1(i,3)*Heff1(i,2))-miu0*alpha*gama*(m1(i,1)*m1(i,2)*Heff1(i,2)-m1(i,2)*m1(i,2)*Heff1(i,1)-m1(i,3)*m1(i,3)*Heff1(i,1)+m1(i,1)*m1(i,3)*Heff1(i,3));
    rt1(i+1,2)=-1*gama*miu0*(m1(i,3)*Heff1(i,1)-m1(i,1)*Heff1(i,3))-miu0*alpha*gama*(m1(i,3)*m1(i,2)*Heff1(i,3)-m1(i,3)*m1(i,3)*Heff1(i,2)-m1(i,1)*m1(i,1)*Heff1(i,2)+m1(i,1)*m1(i,2)*Heff1(i,1));
    rt1(i+1,3)=-1*gama*miu0*(m1(i,1)*Heff1(i,2)-m1(i,2)*Heff1(i,1))-miu0*alpha*gama*(m1(i,1)*m1(i,3)*Heff1(i,1)-m1(i,1)*m1(i,1)*Heff1(i,3)-m1(i,2)*m1(i,2)*Heff1(i,3)+m1(i,2)*m1(i,3)*Heff1(i,2));
    m1(i+1,1)=m1(i,1)+dt*rt1(i+1,1);
    m1(i+1,2)=m1(i,2)+dt*rt1(i+1,2);
    m1(i+1,3)=m1(i,3)+dt*rt1(i+1,3);
    m1(i+1,4)=sqrt(m1(i+1,1)^2+m1(i+1,2)^2+m1(i+1,3)^2);
    m1(i+1,1)=m1(i+1,1)/m1(i+1,4);
    m1(i+1,2)=m1(i+1,2)/m1(i+1,4);
    m1(i+1,3)=m1(i+1,3)/m1(i+1,4);
    
    rt2(i+1,1)=-1*gama*miu0*(m2(i,2)*Heff2(i,3)-m2(i,3)*Heff2(i,2))-miu0*alpha*gama*(m2(i,1)*m2(i,2)*Heff2(i,2)-m2(i,2)*m2(i,2)*Heff2(i,1)-m2(i,3)*m2(i,3)*Heff2(i,1)+m2(i,1)*m2(i,3)*Heff2(i,3));
    rt2(i+1,2)=-1*gama*miu0*(m2(i,3)*Heff2(i,1)-m2(i,1)*Heff2(i,3))-miu0*alpha*gama*(m2(i,3)*m2(i,2)*Heff2(i,3)-m2(i,3)*m2(i,3)*Heff2(i,2)-m2(i,1)*m2(i,1)*Heff2(i,2)+m2(i,1)*m2(i,2)*Heff2(i,1));
    rt2(i+1,3)=-1*gama*miu0*(m2(i,1)*Heff2(i,2)-m2(i,2)*Heff2(i,1))-miu0*alpha*gama*(m2(i,1)*m2(i,3)*Heff2(i,1)-m2(i,1)*m2(i,1)*Heff2(i,3)-m2(i,2)*m2(i,2)*Heff2(i,3)+m2(i,2)*m2(i,3)*Heff2(i,2));
    m2(i+1,1)=m2(i,1)+dt*rt2(i+1,1);
    m2(i+1,2)=m2(i,2)+dt*rt2(i+1,2);
    m2(i+1,3)=m2(i,3)+dt*rt2(i+1,3);
    m2(i+1,4)=sqrt(m2(i+1,1)^2+m2(i+1,2)^2+m2(i+1,3)^2);
    m2(i+1,1)=m2(i+1,1)/m2(i+1,4);
    m2(i+1,2)=m2(i+1,2)/m2(i+1,4);
    m2(i+1,3)=m2(i+1,3)/m2(i+1,4);
%Hk
    Hk1(i+1,1)=-Hkxy*m1(i+1,1);  
    Hk1(i+1,2)=Hkxy*m1(i+1,2);  
    Hk1(i+1,3)=Hkz*m1(i+1,3);  
    Hk2(i+1,1)=Hkxy*m2(i+1,1);  
    Hk2(i+1,2)=-Hkxy*m2(i+1,2);  
    Hk2(i+1,3)=Hkz*m2(i+1,3);
%Hex
    Hex1(i+1,1)=Hex*m2(i+1,1);  
    Hex1(i+1,2)=Hex*m2(i+1,2);  
    Hex1(i+1,3)=Hex*m2(i+1,3);
    Hex2(i+1,1)=Hex*m1(i+1,1);  
    Hex2(i+1,2)=Hex*m1(i+1,2);  
    Hex2(i+1,3)=Hex*m1(i+1,3);
%Hp
    Hp0=0.05;
    Hp(i+1,3)=0;
    Hp(i+1,2)=Hp0*cos(w*t(i+1,1))*sin(fai);
    Hp(i+1,1)=Hp0*cos(w*t(i+1,1))*cos(fai);
%Heff
    Heff1(i+1,1)=H0x+Hp(i+1,1)+Hex1(i+1,1)+Hk1(i+1,1);
    Heff1(i+1,2)=H0y+Hp(i+1,2)+Hex1(i+1,2)+Hk1(i+1,2);
    Heff1(i+1,3)=H0z+Hp(i+1,3)+Hex1(i+1,3)+Hk1(i+1,3);
    Heff2(i+1,1)=H0x+Hp(i+1,1)+Hex2(i+1,1)+Hk2(i+1,1);
    Heff2(i+1,2)=H0y+Hp(i+1,2)+Hex2(i+1,2)+Hk2(i+1,2);
    Heff2(i+1,3)=H0z+Hp(i+1,3)+Hex2(i+1,3)+Hk2(i+1,3);  
 %spin current   
    Is1(i+1,1)=m1(i+1,2)*rt1(i+1,3)-m1(i+1,3)*rt1(i+1,2);
    Is1(i+1,2)=m1(i+1,3)*rt1(i+1,1)-m1(i+1,1)*rt1(i+1,3);
    Is1(i+1,3)=m1(i+1,1)*rt1(i+1,2)-m1(i+1,2)*rt1(i+1,1);
    Is2(i+1,1)=m2(i+1,2)*rt2(i+1,3)-m2(i+1,3)*rt2(i+1,2);
    Is2(i+1,2)=m2(i+1,3)*rt2(i+1,1)-m2(i+1,1)*rt2(i+1,3);
    Is2(i+1,3)=m2(i+1,1)*rt2(i+1,2)-m2(i+1,2)*rt2(i+1,1);
    Is3(i+1,1)=m1(i+1,2)*rt2(i+1,3)-m1(i+1,3)*rt2(i+1,2);
    Is3(i+1,2)=m1(i+1,3)*rt2(i+1,1)-m1(i+1,1)*rt2(i+1,3);
    Is3(i+1,3)=m1(i+1,1)*rt2(i+1,2)-m1(i+1,2)*rt2(i+1,1);
    Is4(i+1,1)=m2(i+1,2)*rt1(i+1,3)-m2(i+1,3)*rt1(i+1,2);
    Is4(i+1,2)=m2(i+1,3)*rt1(i+1,1)-m2(i+1,1)*rt1(i+1,3);
    Is4(i+1,3)=m2(i+1,1)*rt1(i+1,2)-m2(i+1,2)*rt1(i+1,1);
    Is=(h/e)*(g11*Is1+g22*Is2+g12*Is3+g21*Is4);    
end
%every 10 points
tn=n/10;
for j=1:1:tn
    mt(j,1)=m1(j*10-9,1);
    mt(j,2)=m1(j*10-9,2);
    mt(j,3)=m1(j*10-9,3);
    mt(j,4)=m2(j*10-9,1);
    mt(j,5)=m2(j*10-9,2);
    mt(j,6)=m2(j*10-9,3); 
    mt(j,7)=mt(j,1)+mt(j,4);
    mt(j,8)=mt(j,2)+mt(j,5);
    mt(j,9)=mt(j,3)+mt(j,6); 
    mt(j,10)=t(j*10-9,1);
    Ist(j,1)=Is(j*10-9,1);
    Ist(j,2)=Is(j*10-9,2);
    Ist(j,3)=Is(j*10-9,3);  
    Is1t(j,1)=Is1(j*10-9,1);
    Is1t(j,2)=Is1(j*10-9,2);
    Is1t(j,3)=Is1(j*10-9,3); 
    Is2t(j,1)=Is2(j*10-9,1);
    Is2t(j,2)=Is2(j*10-9,2);
    Is2t(j,3)=Is2(j*10-9,3); 
    Is3t(j,1)=Is3(j*10-9,1);
    Is3t(j,2)=Is3(j*10-9,2);
    Is3t(j,3)=Is3(j*10-9,3); 
    Is4t(j,1)=Is4(j*10-9,1);
    Is4t(j,2)=Is4(j*10-9,2);
    Is4t(j,3)=Is4(j*10-9,3); 
end
T=1/f;
Tn=T*1e-10/dt;
NI=floor(10*Tn);
Iselectz(:,1)=(h/e)*Is1t(tn-NI:tn,3);
Iselectz(:,2)=(h/e)*Is2t(tn-NI:tn,3);
Iselectz(:,3)=(h/e)*Is3t(tn-NI:tn,3);
Iselectz(:,4)=(h/e)*Is4t(tn-NI:tn,3);
Iselectz(:,5)=Ist(tn-NI:tn,3);
Isztot(fain,1)=sum(Iselectz(:,1))/10;
Isztot(fain,2)=sum(Iselectz(:,2))/10;
Isztot(fain,3)=sum(Iselectz(:,3))/10;
Isztot(fain,4)=sum(Iselectz(:,4))/10;
Isztot(fain,5)=sum(Iselectz(:,5))/10;
Iselectx(:,1)=(h/e)*Is1t(tn-NI:tn,1);
Iselectx(:,2)=(h/e)*Is2t(tn-NI:tn,1);
Iselectx(:,3)=(h/e)*Is3t(tn-NI:tn,1);
Iselectx(:,4)=(h/e)*Is4t(tn-NI:tn,1);
Iselectx(:,5)=Ist(tn-NI:tn,1);
Isxtot(fain,1)=sum(Iselectx(:,1))/10;
Isxtot(fain,2)=sum(Iselectx(:,2))/10;
Isxtot(fain,3)=sum(Iselectx(:,3))/10;
Isxtot(fain,4)=sum(Iselectx(:,4))/10;
Isxtot(fain,5)=sum(Iselectx(:,5))/10;
Iselecty(:,1)=(h/e)*Is1t(tn-NI:tn,2);
Iselecty(:,2)=(h/e)*Is2t(tn-NI:tn,2);
Iselecty(:,3)=(h/e)*Is3t(tn-NI:tn,2);
Iselecty(:,4)=(h/e)*Is4t(tn-NI:tn,2);
Iselecty(:,5)=Ist(tn-NI:tn,2);
Isytot(fain,1)=sum(Iselecty(:,1))/10;
Isytot(fain,2)=sum(Iselecty(:,2))/10;
Isytot(fain,3)=sum(Iselecty(:,3))/10;
Isytot(fain,4)=sum(Iselecty(:,4))/10;
Isytot(fain,5)=sum(Iselecty(:,5))/10;
faid=fain*5;
strm1 = ['E:\AFM\mt-fai',int2str(faid),'-f',int2str(fn),'.xlsx'];
xlswrite (strm1,mt);
clear m1;
clear m2;
clear rt1;
clear rt2;
clear Is1;
clear Is2;
clear Is3;
clear Is4;
clear Iselectx;
clear Iselecty;
clear Iselectz;
clear mt
end
strm4 = ['E:\AFM\Isx-fai',int2str(faid),'.xls'];
xlswrite (strm4,Isxtot);
strm5 = ['E:\AFM\Isy-fai',int2str(faid),'.xls'];
xlswrite (strm5,Isytot);
strm6 = ['E:\AFM\Isz-fai',int2str(faid),'.xls'];
xlswrite (strm6,Isztot);
end