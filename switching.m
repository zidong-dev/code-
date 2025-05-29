%parameters
h=6.626e-34/(2*pi);
e=1.6e-19;
pi=3.1415926;
alpha=0.03;
gama=28e9*2*pi;
Ms=2.79e4;
Hexaf=-40;  % the antiferromagnetic
Hkxy=-1.43; % the xy direction
Hkz=-4.02; % the z direction
eathetax=0*pi/180; % the x easy axis
eathetay=0*pi/180; % the y easy axis 
dt=5e-15;
miu0=1;
n=30000;
for hn=1:1:100
    H0=(hn)*0.1; 
%H0
fai=0*3.14/180;
H0x=H0*cos(fai);
H0y=H0*sin(fai);
H0z=0;
%m0
m1(1,1)=0.964723;
m1(1,2)=0.26326;
m1(1,3)=0;
m1(1,4)=1;
m2(1,1)=-0.95555;
m2(1,2)=-0.294803;
m2(1,3)=0;
m2(1,4)=1; 
%Hk
Hk1(1,1)=-Hkxy*m1(1,1)*cos(eathetax);  
Hk1(1,2)=Hkxy*m1(1,2)*cos(eathetay);  
Hk1(1,3)=Hkz*m1(1,3);  
Hk2(1,1)=Hkxy*m2(1,1)*cos(eathetax);  
Hk2(1,2)=-Hkxy*m2(1,2)*cos(eathetay);  
Hk2(1,3)=Hkz*m2(1,3);
%Hex
Hex1(1,1)=Hexaf*m2(1,1);  
Hex1(1,2)=Hexaf*m2(1,2);  
Hex1(1,3)=Hexaf*m2(1,3);
Hex2(1,1)=Hexaf*m1(1,1);  
Hex2(1,2)=Hexaf*m1(1,2);  
Hex2(1,3)=Hexaf*m1(1,3);
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
Hk1(i+1,1)=-Hkxy*m1(i+1,1)*cos(eathetax);  
Hk1(i+1,2)=Hkxy*m1(i+1,2)*cos(eathetay);  
Hk1(i+1,3)=Hkz*m1(i+1,3);  
Hk2(i+1,1)=Hkxy*m2(i+1,1)*cos(eathetax);  
Hk2(i+1,2)=-Hkxy*m2(i+1,2)*cos(eathetay);  
Hk2(i+1,3)=Hkz*m2(i+1,3);
%Hex
Hex1(i+1,1)=Hexaf*m2(i+1,1);  
Hex1(i+1,2)=Hexaf*m2(i+1,2);  
Hex1(i+1,3)=Hexaf*m2(i+1,3);
Hex2(i+1,1)=Hexaf*m1(i+1,1);  
Hex2(i+1,2)=Hexaf*m1(i+1,2);  
Hex2(i+1,3)=Hexaf*m1(i+1,3);
%Htemp

%Heff
Heff1(i+1,1)=H0x+Hex1(i+1,1)+Hk1(i+1,1);
Heff1(i+1,2)=H0y+Hex1(i+1,2)+Hk1(i+1,2);
Heff1(i+1,3)=H0z+Hex1(i+1,3)+Hk1(i+1,3);
Heff2(i+1,1)=H0x+Hex2(i+1,1)+Hk2(i+1,1);
Heff2(i+1,2)=H0y+Hex2(i+1,2)+Hk2(i+1,2);
Heff2(i+1,3)=H0z+Hex2(i+1,3)+Hk2(i+1,3);
end
%every 10 points
tn=n/10;
clear mt;
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
end
strm1 = ['E:\SOT-trilayers\m-hx-h',int2str(hn),'.xlsx'];
xlswrite (strm1,mt);
end