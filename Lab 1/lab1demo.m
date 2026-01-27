%Numerical experiments with the harmonic oscillator
%Given d2x/dt2+omega^2 x=0, x(0)=1, dx/dt(0)=0
omega=2;
h=0.1;N=100;
t=[0:h:N*h];
x=cos(omega*t);
xp=-omega*sin(omega*t);
E=omega^2*(x.*x)+xp.*xp;
X=[x' xp'];
plot(t,x,t,xp)
xlabel('t'),ylabel('x,xp')
title('Harmonic oscillator, exact trajetories')
pause
plot(x,xp)
axis('equal')
xlabel('x'),ylabel('dx/dt')
title('Harmonic oscillator, exact phase portrait')
grid
pause
plot(t,E)
xlabel('t')
ylabel('E')
title('Harmonic oscillator, exact Energi vs time')
pause
% y=[1,0]';t=0;
% A=[0 1;-omega^2 0];
% y=[1,0]';t=0;
% Y=[y'];T=[t];E=[omega^2];
% for k=1:N
%    y=y+h*A*y;
%    Y=[Y;y'];
%    t=t+h;
%    T=[T;t];
%    E=[E;omega^2*y(1)*y(1)+y(2)*y(2)];
% end
% Yexeul=Y;
% plot(T,Y,'*')
% xlabel('t'),ylabel('y1,y2')
% title('Harmonic oscillator, Explicit Euler trajectories')
% pause
% plot(Y(:,1),Y(:,2),'*')
% axis('equal')
% xlabel('y1'),xlabel('y2')
% title('Harmonic oscillator, Explicit Euler phase portrait')
% pause
% plot(T,E,'*')
% xlabel('t'),ylabel('E')
% title('Harmonic oscillator, Explicit Euler Energi vs time')
% pause
% y=[1,0]';t=0;
% Y=[y'];T=[t];E=[omega^2];
% for k=1:N
%     y=(eye(2)-h*A)\y;
%     Y=[Y;y'];
%     t=t+h;
%     T=[T;t];
%     E=[E;omega^2*y(1)*y(1)+y(2)*y(2)];
% end
% Yimeul=Y;
% plot(T,Y,'o')
% xlabel('t'),ylabel('y1,y2')
% title('Harmonic Oscillator, Implicit Euler trajectories')
% pause
% plot(Y(:,1),Y(:,2),'o')
% axis('equal')
% xlabel('y1'),ylabel('y2')
% title('Harmonic Oscillator, Implicit Euler phase portrait')
% pause
% plot(T,E)
% xlabel('t'),ylabel('E')
% title('Harmonic Oscillator, Implicit Euler Energy vs time')
% pause
u1 = 1;  % Initial position
u2 = 0;  % Initial velocity
U1 = [u1];  % Will grow to store all positions
U2 = [u2];  % Will grow to store all velocities
t=0;
T=[t];
E=[omega^2]; %This is like this because init pos = 1 and vel = 0
for k=1:N
    u1 = U1(end) + h*U2(end);
    u2 = U2(end) - (h * (omega^2) * u1);
    U1=[U1;u1];
    U2=[U2;u2];
    t=t+h;
    T=[T;t];
    E=[E;(omega^2)*(u1^2) + (u2^2)];
end
plot(T,U1,T,U2,'o')
xlabel('t'),ylabel('u1,u2')
title('Harmonic Oscillator, Sympletic Euler trajectories')
pause
plot(U1, U2,'o')
axis('equal')
xlabel('y1'),ylabel('y2')
title('Harmonic Oscillator, Sympletic Euler phase portrait')
pause
plot(T,E)
xlabel('t'),ylabel('E')
title('Harmonic Oscillator, Sympletic Euler Energy vs time')
pause
%IMPLICIT MIDPOINT EULER
y=[1,0]';t=0;
A=[0 1;-omega^2 0];
Y=[y'];
T=[t];
E=[omega^2];
hh = h/2;
Anum = eye(2) + hh*A;
Aden = eye(2) - hh*A;
Amod = Aden\Anum;
for k=1:N
    y = Amod * Y(end,:)';
    Y=[Y;y'];
    t=t+h;
    T=[T;t];
    E=[E;omega^2*y(1)*y(1)+y(2)*y(2)];
end
Yimeul=Y;
plot(T,Y,'o')
xlabel('t'),ylabel('y1,y2')
title('Harmonic Oscillator, Implicit Midpoint trajectories')
pause
plot(Y(:,1),Y(:,2),'o')
axis('equal')
xlabel('y1'),ylabel('y2')
title('Harmonic Oscillator, Implicit Midpoint phase portrait')
pause
plot(T,E)
xlabel('t'),ylabel('E')
title('Harmonic Oscillator, Implicit Midpoint Energy vs time')
pause