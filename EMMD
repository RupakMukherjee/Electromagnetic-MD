program electromagnetic_molecular_dynamics
implicit none
real pi,a,M,epsilon,wpd,L,Lx,Ly,Lz,Vxmax,Vymax,Vzmax,svx,svy,svz,nbar,r,Gamma,Temp,K,fx,fy,fz,KE,PE,TE,B,alpha
REAL t,tmax,tmic,tmag,dt,xdiff,ydiff,zdiff,scale,tau,sumvx,sumvy,sumvz,virial,f,Pressure
integer i,j,N,p
REAL, DIMENSION (1000):: x,y,z,vx,vy,vz,ax,ay,az,ux,uy,uz,avx,avy,avz,v0x,v0y,v0z,Bx,By,Bz
integer,parameter :: seed = 99999999
pi = 4.0*atan(1.0)
call srand(seed)

!====================== User Inputs =============================

!Coupling Parameter...
Gamma = 300.0
Temp = 1.0/Gamma

!Screening Parameter...
K = 1.0

!Mass of the Dust Particles...
M = 1.0

!Free Space Permittivity...
epsilon = 1.0 !8.8541878*10.0**(-12)

!Total Number of Particles...
N = 686

!Areal Density of Particles...
!nbar = float(N)/(2.0*Lx*2.0*Ly)
nbar = 3.0/(4.0*pi) !1.0/pi

!Interparticle Distance...
a = sqrt(1.0/(pi*nbar)) !1.0

!Dust Plasma Frequency...
wpd = sqrt(nbar/(M*epsilon))
!wpd = sqrt(nbar/(M*a*epsilon))   ! 2-D system -- 'a' is present in the denominator...

!Velocity of light
alpha = 1.0

!Normalized System Size...
L = (float(N)/(8.0*nbar))**(1.0/3.0)
Lx = L/a
Ly = L/a
Lz = L/a
!L = sqrt(float(N)/(4.0*nbar))
!Lx = L/a
!Ly = L/a

!Initial Normalized Maximum Velocity of the Particles
Vxmax = 1.0*sqrt(2.0)/a
Vymax = 1.0*sqrt(2.0)/a
Vzmax = 1.0*sqrt(2.0)/a

!Normalized Final Time and Time Steps...
tmax = 100.0!/sqrt(2.0)
dt = 0.01!/sqrt(2.0)
tmic = tmax/2.0
tmag = 3.0*tmax/4.0

svx = 0.0
svy = 0.0
svz = 0.0
virial = 0.0

!====================== Output Filenames ===========================

open(unit=1,file='Initial_Configuration.dat',status='unknown')
open(unit=2,file='Average_Velocity.dat',status='unknown')
open(unit=10,file='Energy.dat',status='unknown')

!====================== Definition of initial state =================

!Definition of the initial random positions and velocities in -Lx to Lx and -Ly to Ly rectangular box... 
DO i = 1, N, 1
x(i) = (rand())*2.0*Lx - Lx
y(i) = (rand())*2.0*Ly - Ly
z(i) = (rand())*2.0*Lz - Lz
vx(i) = (rand())*Vxmax - Vxmax/2.0
svx = svx + vx(i)                                               ! Center of Mass has Zero x Velocity...
vy(i) = (rand())*Vymax - Vymax/2.0
svy = svy + vy(i)                                               ! Center of Mass has Zero y Velocity...
vz(i) = (rand())*Vzmax - Vzmax/2.0
svz = svz + vz(i)                                               ! Center of Mass has Zero z Velocity...
ENDDO

!Definitions of corrected initial velocities...
do i = 1,N,1
vx(i) = vx(i) - svx/float(N)
vy(i) = vy(i) - svy/float(N)
vz(i) = vz(i) - svz/float(N)
enddo

!Calculating the initial accelerations of the particles...
do i = 1,N,1
do j = 1,N,1
ax(i) = 0.0
ay(i) = 0.0
az(i) = 0.0
avx(i) = 0.0
avy(i) = 0.0
avz(i) = 0.0
if (i .ne. j) then
xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0*Lx))*2.0*Lx        ! Minimum Image Convension Introduced...
ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0*Ly))*2.0*Ly        ! Minimum Image Convension Introduced...
zdiff = (z(i)-z(j)) - nint((z(i)-z(j))/(2.0*Lz))*2.0*Lz        ! Minimum Image Convension Introduced...
r = SQRT((xdiff)**2 + (ydiff)**2 + (zdiff)**2)
!r = sqrt((xdiff)**2 + (ydiff)**2)
f = (1.0+k*r)*exp(-k*r)/r**3
fx = (xdiff)*(1.0+k*r)*exp(-k*r)/r**3 !- alpha*(xdiff)*r
ax(i) = ax(i) + fx
fy = (ydiff)*(1.0+k*r)*exp(-k*r)/r**3 !- alpha*(ydiff)*r
ay(i) = ay(i) + fy
fz = (zdiff)*(1.0+k*r)*exp(-k*r)/r**3 !- alpha*(zdiff)*r
az(i) = az(i) + fz
virial = virial + r**2 * f
endif
enddo ! j
enddo ! i

do i = 1,N,1
write(1,*) 0,x(i),y(i),z(i),vx(i),vy(i),vz(i),ax(i),ay(i),az(i)
!write(1,*) 0,x(i),y(i),vx(i),vy(i),ax(i),ay(i)
enddo

!====== MD Time Evolution using velocity verlet method started =========

do t = 0.010,tmax,dt
sumvx = 0.0
sumvy = 0.0
sumvz = 0.0
virial = 0.0
Pressure = 0.0
PE = 0.0
KE = 0.0

!Calculating velocity after time dt/2.0 and position after time dt... 
if (t .le. tmag) then

do i = 1,N,1
ux(i) = vx(i) + dt*ax(i)/2.0
x(i) = x(i) + dt*ux(i)
x(i) = x(i) - (int(x(i)/Lx))*2.0*Lx                         ! Periodic Boundary Condition
uy(i) = vy(i) + dt*ay(i)/2.0
y(i) = y(i) + dt*uy(i)
y(i) = y(i) - (int(y(i)/Ly))*2.0*Ly                         ! Periodic Boundary Condition
uz(i) = vz(i) + dt*az(i)/2.0
z(i) = z(i) + dt*uz(i)
z(i) = z(i) - (int(z(i)/Lz))*2.0*Lz                         ! Periodic Boundary Condition
enddo ! i

!Calculating the acceleration after time dt...
do i = 1,N,1
ax(i) = 0.0
ay(i) = 0.0
az(i) = 0.0
do j = 1,N,1
if (i .ne. j) then                                          ! .ne. is used, PE should be halved...
xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0*Lx))*2.0*Lx     ! Minimum Image Convension Introduced...
ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0*Ly))*2.0*Ly     ! Minimum Image Convension Introduced...
zdiff = (z(i)-z(j)) - nint((z(i)-z(j))/(2.0*Lz))*2.0*Lz     ! Minimum Image Convension Introduced...
r = sqrt((xdiff)**2 + (ydiff)**2 + (zdiff)**2)
!r = sqrt((xdiff)**2 + (ydiff)**2)
f = (1.0+k*r)*EXP(-k*r)/r**3
fx = (xdiff)*(1.0+k*r)*exp(-k*r)/r**3 !- alpha*(xdiff)*r
ax(i) = ax(i) + fx
fy = (ydiff)*(1.0+k*r)*EXP(-k*r)/r**3 !- alpha*(ydiff)*r
ay(i) = ay(i) + fy
fz = (zdiff)*(1.0+k*r)*exp(-k*r)/r**3 !- alpha*(zdiff)*r
az(i) = az(i) + fz
virial = virial + r**2 * f
PE = PE + (exp(-k*r)/(2.0*r))                               ! Calculation of Potential Energy...
endif
enddo ! j
enddo ! i

!Calculating the velocity after time dt...
do i = 1,N,1
sumvx = sumvx + vx(i)                                       ! Check for average x-velocity...
sumvy = sumvy + vy(i)                                       ! Check for average y-velocity...
sumvz = sumvz + vz(i)                                       ! Check for average y-velocity...
vx(i) = ux(i) + ax(i)*dt/2.0
vy(i) = uy(i) + ay(i)*dt/2.0
vz(i) = uz(i) + az(i)*dt/2.0
KE = KE + (vx(i)**2 + vy(i)**2 + vz(i)**2)/2.0              ! Calculation of Kinetic Energy...
!KE = KE + (vx(i)**2 + vy(i)**2)/2.0 
enddo ! i

elseif (t .gt. tmag) then !=============== Magnetic Field Applied ==============================

do i = 1,N,1
ux(i) = vx(i) + dt*(ax(i) + avx(i))/2.0
x(i) = x(i) + dt*ux(i)
x(i) = x(i) - (int(x(i)/Lx))*2.0*Lx                         ! Periodic Boundary Condition
uy(i) = vy(i) + dt*(ay(i) + avy(i))/2.0
y(i) = y(i) + dt*uy(i)
y(i) = y(i) - (int(y(i)/Ly))*2.0*Ly                         ! Periodic Boundary Condition
uz(i) = vz(i) + dt*(az(i) + avz(i))/2.0
z(i) = z(i) + dt*uz(i)
z(i) = z(i) - (int(z(i)/Lz))*2.0*Lz                         ! Periodic Boundary Condition
enddo !i

!Calculating the acceleration after time dt...
do i = 1,N,1
ax(i) = 0.0
ay(i) = 0.0
az(i) = 0.0
Bx(i) = 0.0
By(i) = 0.0
Bz(i) = 0.0
avx(i) = 0.0
avy(i) = 0.0
avz(i) = 0.0
DO j = 1,N,1
if (i .ne. j) then                                          ! .ne. is used, PE should be halved...
xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0*Lx))*2.0*Lx     ! Minimum Image Convension Introduced...
ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0*Ly))*2.0*Ly     ! Minimum Image Convension Introduced...
zdiff = (z(i)-z(j)) - nint((z(i)-z(j))/(2.0*Lz))*2.0*Lz     ! Minimum Image Convension Introduced...
r = sqrt((xdiff)**2 + (ydiff)**2 + (zdiff)**2)
!r = sqrt((xdiff)**2 + (ydiff)**2)
f = (1.0+k*r)*exp(-k*r)/r**3
fx = (xdiff)*(1.0+k*r)*exp(-k*r)/r**3 !- alpha*(xdiff)*r
ax(i) = ax(i) + fx
Bx(i) = Bx(i) + (uy(j)*zdiff - uz(j)*ydiff)/r**3
fy = (ydiff)*(1.0+k*r)*exp(-k*r)/r**3 !- alpha*(ydiff)*r
ay(i) = ay(i) + fy
By(i) = By(i) + (uz(j)*xdiff - ux(j)*zdiff)/r**3
fz = (zdiff)*(1.0+k*r)*exp(-k*r)/r**3 !- alpha*(zdiff)*r
az(i) = az(i) + fz
Bz(i) = Bz(i) + (ux(j)*ydiff - uy(j)*xdiff)/r**3
virial = virial + r**2 * f
PE = PE + (exp(-k*r)/(2.0*r))                               ! Calculation of Potential Energy...
endif
ENDDO ! j
avx(i) = uy(i)*Bz(i) - uz(i)*By(i)
avy(i) = uz(i)*Bx(i) - ux(i)*Bz(i)
avz(i) = ux(i)*By(i) - uy(i)*Bx(i)
enddo ! i

!Calculating the velocity after time dt...
do i = 1,N,1
sumvx = sumvx + vx(i)                                       ! Check for average x-velocity...
sumvy = sumvy + vy(i)                                       ! Check for average y-velocity...
sumvz = sumvz + vz(i)                                       ! Check for average y-velocity...
v0x(i) = ux(i) + ax(i)*dt/2.0
v0y(i) = uy(i) + ay(i)*dt/2.0
v0z(i) = uz(i) + az(i)*dt/2.0
enddo ! i

do i = 1,N,1
vx(i) = v0x(i) + alpha*avx(i)*dt/2.0
vy(i) = v0y(i) + alpha*avy(i)*dt/2.0
vz(i) = v0z(i) + alpha*avz(i)*dt/2.0
enddo ! i

do i = 1,N,1
avx(i) = vy(i)*Bz(i) - vz(i)*By(i)
avy(i) = vz(i)*Bx(i) - vx(i)*Bz(i)
avz(i) = vx(i)*By(i) - vy(i)*Bx(i)
enddo ! i

DO i = 1,N,1
vx(i) = v0x(i) + alpha*avx(i)*dt/2.0
vy(i) = v0y(i) + alpha*avy(i)*dt/2.0
vz(i) = v0z(i) + alpha*avz(i)*dt/2.0
KE = KE + (vx(i)**2 + vy(i)**2 + vz(i)**2)/2.0              ! Calculation of Kinetic Energy...
!KE = KE + (vx(i)**2 + vy(i)**2)/2.0 
ENDDO ! i

endif

!Writing the instantaneous position, velocity, acceleration in files...
do i = 1,N,1
p = int(t/dt)
if (p  .ge. tmax/(2.0*dt) .and. mod(float(p),100.0) == 0.0) then
write(p,*) t,x(i),y(i),z(i),vx(i),vy(i),vz(i),ax(i),ay(i),az(i)
!write(p,*) t,x(i),y(i),vx(i),vy(i),ax(i),ay(i)
endif
enddo

!Total Energy and Pressure...
Pressure = (2.0*KE + virial/2.0)/(3.0*2.0*Lx*2.0*Ly*2.0*Lz) 
!Pressure = (2.0*KE + virial/2.0)/(3.0*2.0*Lx*2.0*Ly) 
TE = KE + PE
write(10,*) t,TE/float(N),KE/float(N),PE/float(N),Pressure


!====================  Thermostat =======================

tau = 10.0*dt
scale = sqrt(1.0 + (dt/tau)*((Temp/(2.0*KE/(3.0*float(N)) ))-1.0))
!scale = sqrt(1.0 + (dt/tau)*((Temp/(KE/(float(N)) ))-1.0))
if (t .le. tmic) then
do i = 1,N,1
vx(i) = scale*vx(i)
vy(i) = scale*vy(i)
vz(i) = scale*vz(i)
enddo
else
do i = 1,N,1
vx(i) = vx(i)
vy(i) = vy(i)
vz(i) = vz(i)
enddo
endif

enddo ! t

 close(1)
 close(10)

end
