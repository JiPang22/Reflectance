program aaa
implicit none
integer ii
complex i/(0.,1.)/,A,B,C,D,RR
real om,lam,N,n1,n2,aa,bb,degree,theta,beta,k1x,k2x,klam,R,pi
parameter(pi=acos(-1.))
parameter(N=10.)     !Number of Layer
parameter(n1=1.45, n2=2.)
parameter(aa=10.) !a = 10 nm
parameter(bb=10.) !b = 10 nm
open(1,file='aa')


!degree=0.
degree=65.

theta=degree*(2.*pi/360.)    !theta is radian
beta=sin(theta)

do ii = 4000,8000


lam=ii*1.e-3     !lam ~ 400 nm - 800 nm
om=2.*pi/lam     !omega_unit ~ c/lamda ~ 3 * 10**8 / 10**-9

k1x=sqrt((n1*(2.*pi/lam))**2-beta**2)
k2x=sqrt((n2*(2.*pi/lam))**2-beta**2)


A=exp(-i*k1x*aa)*(cos(k2x*bb)-(1./2.)*i*((k2x/k1x)+(k1x/k2x))*sin(k2x*bb))
B=exp(i*k1x*aa)*((-1./2.)*i*((k2x/k1x)-(k1x/k2x))*sin(k2x*bb))
C=exp(-i*k1x*aa)*((1./2.)*i*((k2x/k1x)-(k1x/k2x))*sin(k2x*bb))
D=exp(i*k1x*aa)*(cos(k2x*bb)+(1./2.)*i*((k2x/k1x)+(k1x/k2x))*sin(k2x*bb))

klam=acos((1./2.)*(A+B))

RR=(C*conjg(C))/(C*conjg(C)+(sin(klam)/sin(N*klam))**2)
R=sqrt(RR*conjg(RR))

write(1,*) om,R

enddo
end
