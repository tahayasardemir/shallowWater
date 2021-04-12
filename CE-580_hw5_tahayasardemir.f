	  program WaterWaves
c..Taha Yaşar Demir /1881978
c..CE-580 - Homework #5
	  parameter(mx=20001)
	  common/para/ dx,dt,W,Pw(mx),rL,Cf,T,cycles,T_max,Hi,Vi,g,time
	  common/flow/ u(mx),h(mx),a(mx),b(mx),c(mx),Pwt(mx),pi,N
	  common/midstep/ Utmp(mx),Htmp(mx),Atmp(mx),Btmp(mx),Ctmp(mx)

	  open(11,file='cfl.dat')
	  open(12,file='h.dat')
	  open(13,file='surface.dat')

	  call init
	  time = 0.
	  do while(time.le.T_max)
	  	call output(1)
	  	time = time + 0.5*dt
	  	call mid_continuity
	  	call mid_momentum
	  	time = time + 0.5*dt
	  	call continuity
	  	call momentum
	  	print*, time
	  end do

	  call output(2)

	  close(11)
	  close(12)
	  close(13)

	  stop 
	  end
c-----------------------------------------------------------------------	  
	  subroutine init
	  parameter(mx=20001)
	  common/para/ dx,dt,W,Pw(mx),rL,Cf,T,cycles,T_max,Hi,Vi,g,time
	  common/flow/ u(mx),h(mx),a(mx),b(mx),c(mx),Pwt(mx),pi,N
	  common/midstep/ Utmp(mx),Htmp(mx),Atmp(mx),Btmp(mx),Ctmp(mx)

	  pi     = 22./7.
	  T      = 2. ! s
	  cycles = 10. 
	  T_max  = T*cycles ! s
	  dt     = 0.0008 ! s Change it according to grid number
	  W      = 1. ! m - width of channel
	  rL     = 20. ! m - lenght of channel
	  Hi     = 12. ! m - initial water heigth
	  Vi     = 0. ! m/s -initial velocity
	  N      = 2001 ! grid number / Try other Numbers
	  Cf     = 0.005
	  g      = 9.81 ! N.m/s^2 gravitational acceleration
	  dx     = rL/(N-1)
	  do i=1,N
	  	Pw(i)= 2*h(i) + W ! no need
	  	u(i) = Vi
	  	h(i) = Hi
	  	a(i) = u(i)*h(i)
	  	b(i) = h(i)*u(i)**2 + 0.5*g*h(i)**2 ! no need
	  	c(i) = (Cf*u(i)*abs(u(i))*Pw(i))/(2*W) ! no need
	  enddo

	  return 
	  end
c-----------------------------------------------------------------------
	  subroutine mid_continuity
	  parameter(mx=20001)
	  common/para/ dx,dt,W,Pw(mx),rL,Cf,T,cycles,T_max,Hi,Vi,g,time
	  common/flow/ u(mx),h(mx),a(mx),b(mx),c(mx),Pwt(mx),pi,N
	  common/midstep/ Utmp(mx),Htmp(mx),Atmp(mx),Btmp(mx),Ctmp(mx)

      Htmp(N) = Hi + 0.5*sin(2*pi*time/T) ! Resorvoir side is input wave function
 	  Htmp(1) = h(1) - (0.5*dt/dx)*(a(2)-a(1)) ! forward difference for first node
 	  do i=2,N-1
 	  	Htmp(i) = 0.5*(h(i+1)+h(i-1)) - (0.25*dt/dx)*(a(i+1)-a(i-1))
c 	  	print*, i,Htmp(i),h(i+1),h(i-1),a(i+1),a(i-1)
 	  enddo

	  return
	  end

c-----------------------------------------------------------------------
	  subroutine mid_momentum
      parameter(mx=20001)
	  common/para/ dx,dt,W,Pw(mx),rL,Cf,T,cycles,T_max,Hi,Vi,g,time
	  common/flow/ u(mx),h(mx),a(mx),b(mx),c(mx),Pwt(mx),pi,N
	  common/midstep/ Utmp(mx),Htmp(mx),Atmp(mx),Btmp(mx),Ctmp(mx)

	  Atmp(1) = 0. ! boundary condition
	  Utmp(1) = 0. ! boundary condition
	  Btmp(1) = 0.5*g*Htmp(1)**2
	  Pwt(1)  = 2*Htmp(1) + W
	  Ctmp(1) = 0.

	  Atmp(N) = a(N) - (0.5*dt/dx)*(b(N)-b(N-1)) - 0.5*dt*c(N)
	  Utmp(N) = Atmp(N) / Htmp(N)
	  Btmp(N) = Htmp(N)*Utmp(N)**2 + 0.5*g*Htmp(N)**2
	  Pwt(N)  = 2*Htmp(N) + W
	  Ctmp(N) = (Cf*Utmp(N)*abs(Utmp(N))*Pwt(N))/(2*W)
c	  print*, time,Ctmp(N),Utmp(N),Btmp(N),Htmp(N)
	  do i=2,N-1
	    Atmp(i) = 0.5*(a(i+1)+a(i-1))-(0.25*dt/dx)*(b(i+1)-b(i-1))
     +	         -0.5*dt*c(i)
	  	Utmp(i) = Atmp(i)/Htmp(i)
	  	Btmp(i) = Htmp(i)*Utmp(i)**2 + 0.5*g*Htmp(i)**2
	  	Pwt(i)  = 2*Htmp(i) + W
	  	Ctmp(i) = (Cf*Utmp(i)*abs(Utmp(i))*Pwt(i))/(2*W)
c	  	print*, time,i,Atmp(i),Utmp(i),Btmp(i),Ctmp(i)
	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine continuity
	  parameter(mx=20001)
	  common/para/ dx,dt,W,Pw(mx),rL,Cf,T,cycles,T_max,Hi,Vi,g,time
	  common/flow/ u(mx),h(mx),a(mx),b(mx),c(mx),Pwt(mx),pi,N
	  common/midstep/ Utmp(mx),Htmp(mx),Atmp(mx),Btmp(mx),Ctmp(mx)


	  h(1) = h(1) - (dt/dx)*(Atmp(2)-Atmp(1))
	  do i=2,N-1
	  	h(i) = h(i) - (0.5*dt/dx)*(Atmp(i+1)-Atmp(i-1))
	  enddo
	  h(N) = Hi + 0.5*sin(2*pi*time/T)



	  return
	  end
c-----------------------------------------------------------------------
	  subroutine momentum 
	  parameter(mx=20001)
	  common/para/ dx,dt,W,Pw(mx),rL,Cf,T,cycles,T_max,Hi,Vi,g,time
	  common/flow/ u(mx),h(mx),a(mx),b(mx),c(mx),Pwt(mx),pi,N
	  common/midstep/ Utmp(mx),Htmp(mx),Atmp(mx),Btmp(mx),Ctmp(mx)

	  a(1)  = 0.
	  u(1)  = 0.
	  b(1)  = 0.5*g*h(1)**2
	  c(1)  = 0.
	  Pw(1) = 2*h(1) + W 

	  a(N) = a(N) - (dt/dx)*(Btmp(N)-Btmp(N-1)) - dt*Ctmp(N)
	  u(N) = a(N)/h(N)
	  b(N) = h(N)*u(N)**2 + 0.5*g*h(N)**2
	  Pw(N)= 2*h(N) + W
	  c(N) = (Cf*u(N)*abs(u(N))*Pw(N-1))/(2*W)

c	  print*, time,a(N),b(N),u(N),h(N),Btmp(N),Btmp(N-1)

	  do i=2,N-1
	  	a(i) = a(i) - (0.5*dt/dx)*(Btmp(i+1)-Btmp(i-1)) - dt*Ctmp(i)
c	  	print*, Btmp(i+1),Btmp(i-1),Ctmp(i)
	  	u(i) = a(i) / h(i)
c	  	print*, i,time,u(i),a(i),h(i)
	  	b(i) = h(i)*u(i)**2 + 0.5*g*h(i)**2
	  	Pw(i)= 2*h(i) + W
	    c(i) = (Cf*u(i)*abs(u(i))*Pw(i))/(2*W)
	  enddo
	  
	  return
	  end
c-----------------------------------------------------------------------
	  subroutine output(m)
	  parameter(mx=20001)
	  common/para/ dx,dt,W,Pw(mx),rL,Cf,T,cycles,T_max,Hi,Vi,g,time
	  common/flow/ u(mx),h(mx),a(mx),b(mx),c(mx),Pwt(mx),pi,N
	  common/midstep/ Utmp(mx),Htmp(mx),Atmp(mx),Btmp(mx),Ctmp(mx)
	  real cfl(N),x

	  do i=1,N
	  	cfl(i) = dx/(abs(u(i))+sqrt(g*h(i)))
	  enddo

	  write(11,*) time,minval(cfl)
	  write(12,*) time,h(1),h(N)

	  x = 0.
	  if (m.eq.2) then
	  	do i=1,N
	  		write(13,*) x,h(i)
	  		x = x + dx 
	  	enddo
	  endif

	  return
	  end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------