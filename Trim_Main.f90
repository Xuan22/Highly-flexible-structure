!! 10/09/2017 Xuan
program Trim_Main
	use types
	use str_constants
	use array_inv
	use Coupling
	implicit none
	!-------------------------------------------------------------------------------------------
	!Trim variables
	!-------------------------------------------------------------------------------------------
	integer, parameter :: loops = 10
	integer, parameter :: model = 1
	real(dp) :: ang1
	real(dp) :: ang2
	real(dp) :: freq
	real(dp) :: dAng1
	real(dp) :: dAng2
	real(dp) :: dFreq
	
	real(dp) :: weit
	real(dp) :: distL
	real(dp) :: Xcg
	real(dp) :: Ycg
	
	real(dp) :: resF1
	real(dp) :: resF2
	real(dp) :: resF3
	real(dp) :: TdAng1
	real(dp) :: TdAng2
	real(dp) :: TdFreq

	real(dp), dimension(3,3) :: Jacob
	real(dp), dimension(3,3) :: invJacob
	real(dp) :: detJ

	integer :: it
	integer :: flag
	real(dp) :: omega
	real(dp) :: tStep
	real(dp), dimension(1, Na) :: totalF
	real(dp), dimension(loops, Na) :: lift
	real(dp), dimension(loops, 1) :: thrust
	real(dp), dimension(loops, 1) :: loopAng1
	real(dp), dimension(loops, 1) :: loopAng2
	real(dp), dimension(loops, 1) :: loopFreq

	!geometric parameters
	weit  = 0.62
	distL = 0.05
	Xcg   = 0.01
	Ycg   = 0.002
!!-------------------------------------------------------------------------------------------------
  if(model==1) then
	!set initial values
	ang1 = .2_dp
	ang2 = .2_dp
	freq = 15._dp
	thrust = 0._dp
	do it = 1, loops
	  print *, '--------------------------------', it,'----------------------------------'
		omega = 2.*pi*freq
		tStep = (1./freq)/20.
		call CouplingAna(it, loops, omega, tStep, totalF)
		
		loopAng1(it,1) = ang1
		loopAng2(it,1) = ang2
		loopFreq(it,1) = freq
		lift(it,:) = totalF(1,:)
		do tt = 1, Nts
			thrust(it,1) = thrust(it,1) + totalF(1,tt) *tStep *freq
		end do
		resF1 = 2.*thrust(it,1) -weit
		resF2 = thrust(it,1)*distL*(sin(ang1)+sin(ang2)) -weit*Xcg
		resF3 = thrust(it,1)*distL*(-cos(ang1)+cos(ang2)) -weit*Ycg
		
		! update angle and freq for next iteration
		if(it==1) then
			dAng1 = 0.01_dp
			dAng2 = 0.01_dp
			dFreq = 0.5_dp
		else
			TdAng1 = (thrust(it,1)-thrust(it-1,1)) /dAng1
			TdAng2 = (thrust(it,1)-thrust(it-1,1)) /dAng2
			TdFreq = (thrust(it,1)-thrust(it-1,1)) /dFreq
			Jacob(1,1) = 0._dp
			Jacob(2,1) = 2._dp*thrust(it,1)*distL*cos(ang1)
			Jacob(3,1) = 2._dp*thrust(it,1)*distL*sin(ang1)
			Jacob(1,2) = 0._dp
			Jacob(2,2) = 2._dp*thrust(it,1)*distL*cos(ang2)
			Jacob(3,2) =-2._dp*thrust(it,1)*distL*sin(ang2)
			Jacob(1,3) = 2._dp*TdFreq
			Jacob(2,3) = TdFreq*distL*(sin(ang1)+sin(ang2))
			Jacob(3,3) = TdFreq*distL*(-cos(ang1)+cos(ang2))
			call inv22(Jacob,idenMat3,invJacob)
!			detJ = -Jacob(1,2)*Jacob(2,1)
!			invJacob(1,1) = 2.*TdFreq*distL*sin(ang)/detJ
!			invJacob(2,1) =-2.*thrust(it,1)*distL*cos(ang)/detJ
!			invJacob(1,2) =-2.*TdFreq/detJ
!			invJacob(1,1) = 0._dp
			! 
			dAng1 = -invJacob(1,1)*resF1 -invJacob(1,2)*resF2 -invJacob(1,3)*resF3
			dAng2 = -invJacob(2,1)*resF1 -invJacob(2,2)*resF2 -invJacob(2,3)*resF3
			dFreq = -invJacob(3,1)*resF1 -invJacob(3,2)*resF2 -invJacob(3,3)*resF3
		endif
		print *, 'thrust:',thrust(it,1)
		print *, 'ang:',ang1,ang2,'freq:',freq
		print *, 'dAng:',dAng1,dAng2,'dFreq:',dFreq
		
		ang1 = ang1 +dAng1
		ang2 = ang2 +dAng2
		freq = freq +dFreq
	end do
  elseif(model==0) then
	! For Open Loop Analysis
	freq = 10._dp
	thrust = 0._dp
	do it = 1, loops
		print *, '--------------------------------', it,'----------------------------------'
		freq = freq + 0.5
		omega = 2.*pi*freq
		tStep = (1./freq)/20.
		call CouplingAna(it, loops, omega, tStep, totalF)
		
		loopFreq(it,1) = freq
		lift(it,:) = totalF(1,:)
		do tt = 1, Nts
			thrust(it,1) = thrust(it,1) + totalF(1,tt) *tStep*freq
		end do
		print *, 'thrust:',thrust(it,1)
	end do
  end if
	!------------------------------------------------------------
	open(unit=1000,file='finalLift.dat')
	do i = 1, Nts
		write(1000,*) i, totalF(1,i)
	end do
	close(unit=1000)
	!------------------------------------------------------------
	open(unit=2000,file='looplift.dat')
	do it = 1, loops
		do j = 1, Na
		write(2000,*) (it-1)*Na+j, lift(it,Na)
		end do
	end do
	close(unit=2000)
	!------------------------------------------------------------
	open(unit=3000,file='thrust.dat')
	do it = 1, loops
		write(3000,*) it, thrust(it,1)
	end do
	close(unit=3000)
	!------------------------------------------------------------
	open(unit=4000,file='angle1.dat')
	do it = 1, loops
		write(4000,*) it, loopAng1(it,1)
	end do
	close(unit=4000)
	!------------------------------------------------------------
	open(unit=4000,file='angle2.dat')
	do it = 1, loops
		write(4000,*) it, loopAng2(it,1)
	end do
	close(unit=4000)
	!------------------------------------------------------------
	open(unit=5000,file='frequency.dat')
	do it = 1, loops
		write(5000,*) it, loopFreq(it,1)
	end do
	close(unit=5000)
	!------------------------------------------------------------
end program Trim_Main



