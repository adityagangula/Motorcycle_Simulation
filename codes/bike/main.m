%generalised coords in order 
syms x(t) y(t) z(t) yaw(t) roll(t) pitch(t) delta(t) extr(t) extf(t) thetar(t) thetaf(t)

%other important values (some may be unused)
syms swingarmang(t) xi(t) frontWcamber(t) kinsteeringang(t)
syms sideslipangr(t) sideslipanglef(t) pneumatictrailr(t) pneumatictrailf(t)
%xi - angle between steering assembly plane (F) and front wheel plane (G)

syms q [11 1]
q(1)=x;
q(2)=y;
q(3)=z;
q(4)=yaw;
q(5)=roll;
q(6)=pitch;
q(7)=delta;
q(8)=extr;
q(9)=extf;
q(10)=thetar;
q(11)=thetaf;
dq = diff(q,t);
syms Lag [11 1]
syms eqns [11 1]
syms Y %for ODE solver

epsilon=0.43; %steering head angle 
kr=150000;kf=5000/0.12; % spring constants for rear and front suspension
drear=7375;dfront=1062.5; % damping coefficients for suspensions
daero=0.28; %aero damping coeff
ksr=92000;ksf=84400; % lateral stiffnesses of tires
lr=0.3;lf=0.18; % rest length of rear and front suspension
l1=0.25;l2=0.5; % length from rear frame CoM to swingarm joint and steering column joint respectively
l3=0.265;l4=0.265; % length from swingarm joint to swingarm CoM and length from CoM to rear wheel respectively
l5=0.098;l6=0.268; % length from steering joint to steering CoM and length from CoM to where suspension begins respectively
l7=0.122; % length from where suspension begins to CoM of iverted fork
l8=0.1; % length from swingarm joint to where rear suspension meets rear frame
Rr=0.317;Rf=0.292; % rolling radius rear and front
tr=0.097;tf=0.062; % radius of cross section rear and front
Ir=0.66;If=0.47; % Moment of inertia tires
g=9.8;
Ma=223; % rear frame mass including rider
Mb=10.0; % swingarm mass 
Mc=7.0; % front unsprung mass
Md=16.2; % rear wheel mass
Mg=12.0; % front wheel mass
Mf=8.75; % steering column mass

%important values definitions
xi(t) = epsilon + pitch(t);
frontWcamber(t) = asin(cos(delta)*sin(roll) - cos(roll)*sin(delta)*sin(xi));
kinsteeringang(t) = atan((sin(delta)*cos(xi))/(cos(roll)*cos(delta) + sin(delta)*sin(roll)*sin(xi)));
swingarmang(t) = acos(((l3+l4)^2 + l8^2 - (lr + extr)^2)/(2*(l3+l4)*l8));

%transformation matrices
syms mn dm ad ba ea fe gf hg mh [3 3]

mn(1,1)=cos(yaw);mn(1,2)=0;mn(1,3)=-sin(yaw);
mn(2,1)=0;mn(2,2)=1;mn(2,3)=0;
mn(3,1)=sin(yaw);mn(3,2)=0;mn(3,3)=cos(yaw);

dm(1,1)=1;dm(1,2)=0;dm(1,3)=0;
dm(2,1)=0;dm(2,2)=cos(roll);dm(2,3)=sin(roll);
dm(3,1)=0;dm(3,2)=-sin(roll);dm(3,3)=cos(roll);

ad(1,1)=cos(pitch);ad(1,2)=sin(pitch);ad(1,3)=0;
ad(2,1)=-sin(pitch);ad(2,2)=cos(pitch);ad(2,3)=0;
ad(3,1)=0;ad(3,2)=0;ad(3,3)=1;

ba(1,1)=cos(swingarmang);ba(1,2)=sin(swingarmang);ba(1,3)=0;
ba(2,1)=-sin(swingarmang);ba(2,2)=cos(swingarmang);ba(2,3)=0;
ba(3,1)=0;ba(3,2)=0;ba(3,3)=1;

ea(1,1)=cos(epsilon);ea(1,2)=sin(epsilon);ea(1,3)=0;
ea(2,1)=-sin(epsilon);ea(2,2)=cos(epsilon);ea(2,3)=0;
ea(3,1)=0;ea(3,2)=0;ea(3,3)=1;

fe(1,1)=cos(delta);fe(1,2)=0;fe(1,3)=-sin(delta);
fe(2,1)=0;fe(2,2)=1;fe(2,3)=0;
fe(3,1)=sin(delta);fe(3,2)=0;fe(3,3)=cos(delta);

gf(1,1)=cos(xi);gf(1,2)=-sin(xi);gf(1,3)=0;
gf(2,1)=sin(xi);gf(2,2)=cos(xi);gf(2,3)=0;
gf(3,1)=0;gf(3,2)=0;gf(3,3)=1;

hg(1,1)=1;hg(1,2)=0;hg(1,3)=0;
hg(2,1)=0;hg(2,2)=cos(frontWcamber);hg(2,3)=-sin(frontWcamber);
hg(3,1)=0;hg(3,2)=sin(frontWcamber);hg(3,3)=cos(frontWcamber);

mh(1,1)=cos(kinsteeringang);mh(1,2)=0;mh(1,3)=sin(kinsteeringang);
mh(2,1)=0;mh(2,2)=1;mh(2,3)=0;
mh(3,1)=-sin(kinsteeringang);mh(3,2)=0;mh(3,3)=cos(kinsteeringang);

nm=inv(mn);md=inv(dm);da=inv(ad);ab=inv(ba);ae=inv(ea);ef=inv(fe);fg=inv(gf);gh=inv(hg);hm=inv(mh);

% coords of 6 bodies CoM in N(ground/universal) frame
% A - rear frame
% B - swimngarm
% C - front unsprung mass
% D - rear wheel
% G - front wheel
% F - steering column
syms A B D F C G [1 3] %#ok<*NASGU> 

A(1)=x;A(2)=y;A(3)=z;
B = A - [l1 0 0]*nm*md*da - [l3 0 0]*nm*md*da*ab; %#ok<*MINV> 
D = B - [l4 0 0]*nm*md*da*ab;
F = A + [l2 0 0]*nm*md*da - [0 l5 0]*nm*md*da*ae*ef;
C = F - [0 l6+l7+extf 0]*nm*md*da*ae*ef;
G = F - [0 l6+lf+extf 0]*nm*md*da*ae*ef;

Ax=A(1);Ay=A(2);Az=A(3);
Bx=B(1);By=B(2);Bz=B(3);
Dx=D(1);Dy=D(2);Dz=D(3);
Fx=F(1);Fy=F(2);Fz=F(3);
Cx=x(t) - (cos(delta(t))*((l5*cos(pitch(t))*cos(roll(t))*sin(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(roll(t))*sin(pitch(t))*cos(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) - (cos(delta(t))*((cos(pitch(t))*cos(roll(t))*sin(epsilon)*(l6 + l7 + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (cos(roll(t))*sin(pitch(t))*cos(epsilon)*(l6 + l7 + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) - (sin(delta(t))*sin(roll(t))*(l6 + l7 + extf(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) - (l5*sin(delta(t))*sin(roll(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l2*cos(pitch(t))*cos(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (l2*sin(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2));
Cy=y(t) - (l2*cos(yaw(t))*sin(pitch(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (cos(pitch(t))*cos(roll(t))*cos(epsilon)*(l6 + l7 + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l2*cos(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (cos(roll(t))*sin(pitch(t))*sin(epsilon)*(l6 + l7 + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) - (l5*cos(pitch(t))*cos(roll(t))*cos(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(roll(t))*sin(pitch(t))*sin(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2));
Cz=z(t) - (sin(delta(t))*((l5*cos(pitch(t))*cos(roll(t))*sin(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(roll(t))*sin(pitch(t))*cos(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) - (sin(delta(t))*((cos(pitch(t))*cos(roll(t))*sin(epsilon)*(l6 + l7 + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (cos(roll(t))*sin(pitch(t))*cos(epsilon)*(l6 + l7 + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) + (cos(delta(t))*sin(roll(t))*(l6 + l7 + extf(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(delta(t))*sin(roll(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l2*cos(roll(t))*sin(yaw(t)))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2));
Gx=x(t) - (cos(delta(t))*((l5*cos(pitch(t))*cos(roll(t))*sin(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(roll(t))*sin(pitch(t))*cos(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) - (cos(delta(t))*((cos(pitch(t))*cos(roll(t))*sin(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (cos(roll(t))*sin(pitch(t))*cos(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) - (sin(delta(t))*sin(roll(t))*(l6 + lf + extf(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) - (l5*sin(delta(t))*sin(roll(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l2*cos(pitch(t))*cos(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (l2*sin(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2));
Gy=y(t) - (l2*cos(yaw(t))*sin(pitch(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (cos(pitch(t))*cos(roll(t))*cos(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l2*cos(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (cos(roll(t))*sin(pitch(t))*sin(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) - (l5*cos(pitch(t))*cos(roll(t))*cos(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(roll(t))*sin(pitch(t))*sin(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2));
Gz=z(t) - (sin(delta(t))*((l5*cos(pitch(t))*cos(roll(t))*sin(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(roll(t))*sin(pitch(t))*cos(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) - (sin(delta(t))*((cos(pitch(t))*cos(roll(t))*sin(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (cos(roll(t))*sin(pitch(t))*cos(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) + (cos(delta(t))*sin(roll(t))*(l6 + lf + extf(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(delta(t))*sin(roll(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l2*cos(roll(t))*sin(yaw(t)))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2));

% sideslip angles
sideslipangr(t) = acos(([diff(Dx,t) diff(Dy,t) diff(Dz,t)]*nm*md*[1;1;1])/((sqrt((diff(Dx,t))^2 + (diff(Dy,t))^2 + (diff(Dz,t))^2))));
sideslipangf(t) = acos(([diff(Gx,t) diff(Gy,t) diff(Gz,t)]*nm*md*[1;1;1])/((sqrt((diff(Gx,t))^2 + (diff(Gy,t))^2 + (diff(Gz,t))^2))));

% contact points rear and front
syms Cr Cf [1 3]

Cr = D + [0 -Rr*cos(roll) Rr*sin(roll)]*nm*md;
Cf = G + [0 -Rf*cos(frontWcamber) Rf*sin(frontWcamber)]*nm*md*da*ae*ef*fg;
Crx = x(t) + ((l3*cos(yaw(t))*sin(pitch(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (l3*cos(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)))*sqrt((- l3^4 - 4*l3^3*l4 - 6*l3^2*l4^2 + 2*l3^2*l8^2 + 2*l3^2*lr^2 + 4*l3^2*lr*extr(t) + 2*l3^2*extr(t)^2 - 4*l3*l4^3 + 4*l3*l4*l8^2 + 4*l3*l4*lr^2 + 8*l3*l4*lr*extr(t) + 4*l3*l4*extr(t)^2 - l4^4 + 2*l4^2*l8^2 + 2*l4^2*lr^2 + 4*l4^2*lr*extr(t) + 2*l4^2*extr(t)^2 - l8^4 + 2*l8^2*lr^2 + 4*l8^2*lr*extr(t) + 2*l8^2*extr(t)^2 - lr^4 - 4*lr^3*extr(t) - 6*lr^2*extr(t)^2 - 4*lr*extr(t)^3 - extr(t)^4)/(4*l8^2*(l3 + l4)^2)) + ((l4*cos(yaw(t))*sin(pitch(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (l4*cos(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)))*sqrt((- l3^4 - 4*l3^3*l4 - 6*l3^2*l4^2 + 2*l3^2*l8^2 + 2*l3^2*lr^2 + 4*l3^2*lr*extr(t) + 2*l3^2*extr(t)^2 - 4*l3*l4^3 + 4*l3*l4*l8^2 + 4*l3*l4*lr^2 + 8*l3*l4*lr*extr(t) + 4*l3*l4*extr(t)^2 - l4^4 + 2*l4^2*l8^2 + 2*l4^2*lr^2 + 4*l4^2*lr*extr(t) + 2*l4^2*extr(t)^2 - l8^4 + 2*l8^2*lr^2 + 4*l8^2*lr*extr(t) + 2*l8^2*extr(t)^2 - lr^4 - 4*lr^3*extr(t) - 6*lr^2*extr(t)^2 - 4*lr*extr(t)^3 - extr(t)^4)/(4*l8^2*(l3 + l4)^2)) - (((l3*cos(pitch(t))*cos(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (l3*sin(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)))*(l3^2 + 2*l3*l4 + l4^2 + l8^2 - lr^2 - 2*lr*extr(t) - extr(t)^2))/(2*l8*(l3 + l4)) - (((l4*cos(pitch(t))*cos(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (l4*sin(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)))*(l3^2 + 2*l3*l4 + l4^2 + l8^2 - lr^2 - 2*lr*extr(t) - extr(t)^2))/(2*l8*(l3 + l4)) - (Rr*sin(roll(t))*sin(yaw(t)))/(cos(yaw(t))^2 + sin(yaw(t))^2) - (l1*cos(pitch(t))*cos(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (l1*sin(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2));
Cry = y(t) + ((l3*cos(pitch(t))*cos(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (l3*sin(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)))*sqrt((- l3^4 - 4*l3^3*l4 - 6*l3^2*l4^2 + 2*l3^2*l8^2 + 2*l3^2*lr^2 + 4*l3^2*lr*extr(t) + 2*l3^2*extr(t)^2 - 4*l3*l4^3 + 4*l3*l4*l8^2 + 4*l3*l4*lr^2 + 8*l3*l4*lr*extr(t) + 4*l3*l4*extr(t)^2 - l4^4 + 2*l4^2*l8^2 + 2*l4^2*lr^2 + 4*l4^2*lr*extr(t) + 2*l4^2*extr(t)^2 - l8^4 + 2*l8^2*lr^2 + 4*l8^2*lr*extr(t) + 2*l8^2*extr(t)^2 - lr^4 - 4*lr^3*extr(t) - 6*lr^2*extr(t)^2 - 4*lr*extr(t)^3 - extr(t)^4)/(4*l8^2*(l3 + l4)^2)) + ((l4*cos(pitch(t))*cos(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (l4*sin(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)))*sqrt((- l3^4 - 4*l3^3*l4 - 6*l3^2*l4^2 + 2*l3^2*l8^2 + 2*l3^2*lr^2 + 4*l3^2*lr*extr(t) + 2*l3^2*extr(t)^2 - 4*l3*l4^3 + 4*l3*l4*l8^2 + 4*l3*l4*lr^2 + 8*l3*l4*lr*extr(t) + 4*l3*l4*extr(t)^2 - l4^4 + 2*l4^2*l8^2 + 2*l4^2*lr^2 + 4*l4^2*lr*extr(t) + 2*l4^2*extr(t)^2 - l8^4 + 2*l8^2*lr^2 + 4*l8^2*lr*extr(t) + 2*l8^2*extr(t)^2 - lr^4 - 4*lr^3*extr(t) - 6*lr^2*extr(t)^2 - 4*lr*extr(t)^3 - extr(t)^4)/(4*l8^2*(l3 + l4)^2)) - (Rr*cos(roll(t))^2)/(cos(roll(t))^2 + sin(roll(t))^2) + (((l3*cos(yaw(t))*sin(pitch(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (l3*cos(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)))*(l3^2 + 2*l3*l4 + l4^2 + l8^2 - lr^2 - 2*lr*extr(t) - extr(t)^2))/(2*l8*(l3 + l4)) + (((l4*cos(yaw(t))*sin(pitch(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (l4*cos(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)))*(l3^2 + 2*l3*l4 + l4^2 + l8^2 - lr^2 - 2*lr*extr(t) - extr(t)^2))/(2*l8*(l3 + l4)) + (Rr*cos(yaw(t))*sin(roll(t))^2)/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (l1*cos(yaw(t))*sin(pitch(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (l1*cos(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2));
Crz = z(t) + (Rr*cos(roll(t))*sin(roll(t)))/(cos(roll(t))^2 + sin(roll(t))^2) - (l1*cos(roll(t))*sin(yaw(t)))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (l3*cos(roll(t))*sin(yaw(t)))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (l4*cos(roll(t))*sin(yaw(t)))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (Rr*cos(roll(t))*cos(yaw(t))*sin(roll(t)))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2));
Cfx = x(t) - (cos(delta(t))*((l5*cos(pitch(t))*cos(roll(t))*sin(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(roll(t))*sin(pitch(t))*cos(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) - (cos(delta(t))*((cos(pitch(t))*cos(roll(t))*sin(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (cos(roll(t))*sin(pitch(t))*cos(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) - (sin(epsilon + pitch(t))*((sin(epsilon)*((sin(pitch(t))*((Rf*cos(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) - (Rf*cos(yaw(t))*sin(roll(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(pitch(t))^2 + sin(pitch(t))^2) + (Rf*cos(pitch(t))*sin(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(epsilon)^2 + sin(epsilon)^2) - (cos(epsilon)*((cos(pitch(t))*((Rf*cos(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) - (Rf*cos(yaw(t))*sin(roll(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(pitch(t))^2 + sin(pitch(t))^2) - (Rf*sin(pitch(t))*sin(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(epsilon)^2 + sin(epsilon)^2)))/(cos(epsilon + pitch(t))^2 + sin(epsilon + pitch(t))^2) - (cos(epsilon + pitch(t))*((cos(delta(t))*((sin(epsilon)*((cos(pitch(t))*((Rf*cos(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) - (Rf*cos(yaw(t))*sin(roll(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(pitch(t))^2 + sin(pitch(t))^2) - (Rf*sin(pitch(t))*sin(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(epsilon)^2 + sin(epsilon)^2) + (cos(epsilon)*((sin(pitch(t))*((Rf*cos(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) - (Rf*cos(yaw(t))*sin(roll(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(pitch(t))^2 + sin(pitch(t))^2) + (Rf*cos(pitch(t))*sin(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(epsilon)^2 + sin(epsilon)^2)))/(cos(delta(t))^2 + sin(delta(t))^2) + (sin(delta(t))*((Rf*sin(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) + (Rf*cos(roll(t))*cos(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2)))/(cos(epsilon + pitch(t))^2 + sin(epsilon + pitch(t))^2) - (sin(delta(t))*sin(roll(t))*(l6 + lf + extf(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) - (l5*sin(delta(t))*sin(roll(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l2*cos(pitch(t))*cos(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (l2*sin(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2));
Cfy = y(t) + (cos(epsilon + pitch(t))*((sin(epsilon)*((sin(pitch(t))*((Rf*cos(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) - (Rf*cos(yaw(t))*sin(roll(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(pitch(t))^2 + sin(pitch(t))^2) + (Rf*cos(pitch(t))*sin(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(epsilon)^2 + sin(epsilon)^2) - (cos(epsilon)*((cos(pitch(t))*((Rf*cos(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) - (Rf*cos(yaw(t))*sin(roll(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(pitch(t))^2 + sin(pitch(t))^2) - (Rf*sin(pitch(t))*sin(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(epsilon)^2 + sin(epsilon)^2)))/(cos(epsilon + pitch(t))^2 + sin(epsilon + pitch(t))^2) - (sin(epsilon + pitch(t))*((cos(delta(t))*((sin(epsilon)*((cos(pitch(t))*((Rf*cos(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) - (Rf*cos(yaw(t))*sin(roll(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(pitch(t))^2 + sin(pitch(t))^2) - (Rf*sin(pitch(t))*sin(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(epsilon)^2 + sin(epsilon)^2) + (cos(epsilon)*((sin(pitch(t))*((Rf*cos(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) - (Rf*cos(yaw(t))*sin(roll(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(pitch(t))^2 + sin(pitch(t))^2) + (Rf*cos(pitch(t))*sin(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(epsilon)^2 + sin(epsilon)^2)))/(cos(delta(t))^2 + sin(delta(t))^2) + (sin(delta(t))*((Rf*sin(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) + (Rf*cos(roll(t))*cos(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2)))/(cos(epsilon + pitch(t))^2 + sin(epsilon + pitch(t))^2) - (l2*cos(yaw(t))*sin(pitch(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) - (cos(pitch(t))*cos(roll(t))*cos(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l2*cos(pitch(t))*sin(roll(t))*sin(yaw(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2)) + (cos(roll(t))*sin(pitch(t))*sin(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) - (l5*cos(pitch(t))*cos(roll(t))*cos(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(roll(t))*sin(pitch(t))*sin(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2));
Cfz = z(t) + (cos(delta(t))*((Rf*sin(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) + (Rf*cos(roll(t))*cos(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) - (sin(delta(t))*((sin(epsilon)*((cos(pitch(t))*((Rf*cos(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) - (Rf*cos(yaw(t))*sin(roll(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(pitch(t))^2 + sin(pitch(t))^2) - (Rf*sin(pitch(t))*sin(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(epsilon)^2 + sin(epsilon)^2) + (cos(epsilon)*((sin(pitch(t))*((Rf*cos(roll(t))*sqrt(1 - (cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t)))^2))/(cos(roll(t))^2 + sin(roll(t))^2) - (Rf*cos(yaw(t))*sin(roll(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(pitch(t))^2 + sin(pitch(t))^2) + (Rf*cos(pitch(t))*sin(yaw(t))*(cos(delta(t))*sin(roll(t)) - sin(epsilon + pitch(t))*cos(roll(t))*sin(delta(t))))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2))))/(cos(epsilon)^2 + sin(epsilon)^2)))/(cos(delta(t))^2 + sin(delta(t))^2) - (sin(delta(t))*((l5*cos(pitch(t))*cos(roll(t))*sin(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(roll(t))*sin(pitch(t))*cos(epsilon))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) - (sin(delta(t))*((cos(pitch(t))*cos(roll(t))*sin(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (cos(roll(t))*sin(pitch(t))*cos(epsilon)*(l6 + lf + extf(t)))/((cos(pitch(t))^2 + sin(pitch(t))^2)*(cos(epsilon)^2 + sin(epsilon)^2)*(cos(roll(t))^2 + sin(roll(t))^2))))/(cos(delta(t))^2 + sin(delta(t))^2) + (cos(delta(t))*sin(roll(t))*(l6 + lf + extf(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l5*cos(delta(t))*sin(roll(t)))/((cos(delta(t))^2 + sin(delta(t))^2)*(cos(roll(t))^2 + sin(roll(t))^2)) + (l2*cos(roll(t))*sin(yaw(t)))/((cos(roll(t))^2 + sin(roll(t))^2)*(cos(yaw(t))^2 + sin(yaw(t))^2));

% center of mass
M = Ma+Mb+Mc+Md+Mg+Mf;
COMx = (Ma*Ax + Mb*Bx + Md*Dx + Mf*Fx + Mc*Cx + Mg*Gx)/M;
COMy = (Ma*Ay + Mb*By + Md*Dy + Mf*Fy + Mc*Cy + Mg*Gy)/M;
COMz = (Ma*Az + Mb*Bz + Md*Dz + Mf*Fz + Mc*Cz + Mg*Gz)/M;
COM = [COMx COMy COMz];

fw=0.02; %rolling resistance coeff
maxsideslipr=0.28;maxsideslipf = 0.28;
maxpneumatictrailf=0.03;maxpneumatictrailr=0.03;

% some more geometry parameters
b = abs(COMx - Crx);
p = abs(Crx - Cfx);
h = COMy;
dr = fw*Rr;
df = fw*Rf;
sr = tr*tan(roll);
sf = tf*tan(roll);
    
pneumatictrailr(t) = maxpneumatictrailr*(1 - (sideslipangr/maxsideslipr));
pneumatictrailf(t) = maxpneumatictrailf*(1 - (sideslipangf/maxsideslipf));

% FORCES & WORKFUNCTIONS %
% rear wheel
    % thrust acts in in +x dir in D frame
thrust = 400;
Wthrust = [thrust 0 0]*nm*md*[Crx;Cry;Crz];
    % normal reactional force +y dir in N frame acts at dr preceding Cr
Nr = M*g*(p-b)/p  + thrust*h/p;
WNr = [0 Nr 0]*[Crx;Cry;Crz];
    % longitudinal slip calc
Bk=1.2889;
Ck=0.8480;
Dk=1.2;
Ek=-0.6361;
kk = Dk*Ck*Bk;
longislip = thrust/(kk*Nr); % linear fit
% longislip = erfinv(thrust/(Dk*Nr)); % rather insulting way to find inverse to pacejka's magic formula
    %omega rear calc
omegar = sqrt(((diff(Dx,t))^2 + (diff(Dy,t))^2 + (diff(Dz,t))^2))*(1+longislip)/Rr;
    % Lateral force acts in +z dir in D frame
Ds=1.1;
Bl=7;
Cl=1.3;
El=-0.3;
kl = Ds*Cl*Bl;
Bp=0.714;
Cp=1.4;
Ep=-2;
kp = Ds*Cp*Bp;
latforcer = Ds*Nr*(sin(Cl*atan(Bl*sideslipangr - El*(Bl*sideslipangr-atan(Bl*sideslipangr)))) + sin(Cp*atan(Bp*roll - Ep*(Bp*roll-atan(Bp*roll))))); 
% latforcer = (kl*sideslipangr + kp*roll)*Nr;
spr = latforcer/ksr;
Wlatforcer = [0 0 latforcer]*nm*md*[Crx;Cry;Crz];
    % torque due to thrust/brake affects yaw
Tthrustr = (sr+spr)*thrust;
Mthrustr = Tthrustr*yaw;
    % torque due to pneumatic trail affects yaw
Ttrailr = pneumatictrailr*latforcer;
Mtrailr = Ttrailr*yaw;
    % torque due to rolling resistance affects swingarmangle
Trollresr = dr*Nr;
Mrollresr = Trollresr*swingarmang;
    % torque due to overturning affects roll
Toverturnr = -spr*Nr;
Moverturnr = Toverturnr*roll;


% front wheel
    %normal reactional force +y dir in N frame acts at df preceding Cf
Nf = M*g*b/p  - thrust*h/p;
WNf = [0 Nf 0]*[Cfx;Cfy;Cfz];
    % lateral force acts in +z dir in G frame
latforcef = Ds*Nr*(sin(Cl*atan(Bl*sideslipangf - El*(Bl*sideslipangf-atan(Bl*sideslipangf)))) + sin(Cp*atan(Bp*frontWcamber - Ep*(Bp*frontWcamber-atan(Bp*frontWcamber)))));
% latforcef = (kl*sideslipangr + kp*roll)*Nf;
spf = latforcef/ksf;
Wlatforcef = [0 0 latforcef]*nm*md*da*ae*ef*fg*[Cfx;Cfy;Cfz];
    % torque due to pneumatic trail affects steeringangle 
Ttrailf = pneumatictrailf*latforcef;
Mtrailf = Ttrailf*delta;
    % torque due to thrust/brake affects steeringangle
Tthrustf = (sr+spr)*thrust;
Mthrustf = Tthrustf*delta;
    % torque due to rolling resistance affects pitch
Trollresf = df*Nf;
Mrollresf = Trollresf*pitch;
    % torque due to overturning affects frontwheelcamber
Toverturnf = -spf*Nf;
Moverturnf = Toverturnf*frontWcamber;
omegaf = omegar;

Workfunc = Wthrust + WNr + Wlatforcer + Mthrustr + Mtrailr + Mrollresr + Moverturnr + WNf + Wlatforcef + Mthrustf + Mtrailf + Mrollresf + Moverturnf;

% dissipative functions
    % aerodynamic drag
Daero = 0.5*daero*(((diff(Ax,t))^2 + (diff(Ay,t))^2 + (diff(Az,t))^2));
    % rear susp
Dsuspr = 0.5*drear*diff(extr,t)^2;
    % front susp
Dsuspf = 0.5*dfront*diff(extf,t)^2;
Diss = Daero+Dsuspr+Dsuspf;

KE = 0.5*(Ma*((diff(Ax,t))^2 + (diff(Ay,t))^2 + (diff(Az,t))^2) + Mb*((diff(Bx,t))^2 + (diff(By,t))^2 + (diff(Bz,t))^2) + Md*((diff(Dx,t))^2 + (diff(Dy,t))^2 + (diff(Dz,t))^2) + Mf*((diff(Fx,t))^2 + (diff(Fy,t))^2 + (diff(Fz,t))^2) + Mc*((diff(Cx,t))^2 + (diff(Cy,t))^2 + (diff(Cz,t))^2) + Mg*((diff(Gx,t))^2 + (diff(Gy,t))^2 + (diff(Gz,t))^2) + Ir*omegar^2 + If*omegaf^2);
PE = g*(Ma*Ay + Mb*By + Md*Dy + Mf*Fy + Mc*Cy + Mg*Gy) + 0.5*(kr*extr^2 + kf*extf^2);
L = KE-PE;

for i=1:11
    Lag(i) = diff(diff(L,dq(i)),t) - diff(L,q(i)) + diff(Diss,dq(i)) == diff(Workfunc,q(i));
end

% [VF,Subs] = odeToVectorField(Lag(1));
% Subs
% odefcn = matlabFunction(VF, 'Vars',{t,Y});
% [X,Y] = ode45(odefcn, [0 5], [pi/2 0 pi/2 0]);





