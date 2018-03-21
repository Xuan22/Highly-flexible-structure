!! 03/5/2017 XUAN
module Aero_Dyn
	use types 
	use str_constants
	implicit none
	!! parameters for kinematics
	real(dp), dimension(1,Na) :: alpha10, beta10, alpha_dot, beta_dot 
	real(dp), dimension(nSecX,Na) :: AoA, AoA_dot
	real(dp), dimension(nSecX,Na) :: inAoA, updAoA, updAoA_dot
	real(dp), dimension(nSecX,Na) :: AoA_str
	real(dp), dimension(1,Nmax+1) :: strX
	real(dp), dimension(1,Mmax+1) :: strY
	
	!! parameters for Leading Edge Suction force
	real(dp), dimension(nSecX,Na) :: vPol
	real(dp), dimension(1,Na) :: L_LES
	real(dp), dimension(1,nSecX) :: AA0, AA1
	real(dp), dimension(1,nSecY) :: theta
	real(dp), dimension(1,nSecX) :: elemLocX
	real(dp), dimension(1,Nmax) :: eleStrX
	real(dp), dimension(1,Mmax) :: eleStrY
	
	real(dp), dimension(nSecY,nSecX,Na) :: gamma_theta
	real(dp), dimension(nSecY,nSecX,Na) :: pres_LES
	real(dp), dimension(nSecY,nSecX,Na) :: distrF_LES
	real(dp), dimension(nSecY,nSecX,Na) :: distrF_Wag
	
	!! parameters for Apparent mass
	real(dp), dimension(nSecX,Na) :: t_A0, t_A1, app_v
	real(dp), dimension(nSecX,Na) :: dt_A0V, dt_A1V
		
	real(dp), dimension(nSecY,nSecX,Na) :: gammapp
	real(dp), dimension(nSecY,nSecX,Na) :: pres_app
	real(dp), dimension(1,Na) :: L_app
	
	!! parameters for Wagner Function
	real(dp), dimension(1,Na)  :: L_Wag
	real(dp), dimension(1,Nmax)  :: V_circ
	real(dp), dimension(Nmax,Na) :: Xpart
	real(dp), dimension(Nmax,Na) :: Ypart
	real(dp), dimension(Nmax,Na) :: AoA_eqv
	real(dp), dimension(Nmax,Na) :: cl_Wag
	real(dp), dimension(Nmax,Na) :: Wag_Lift
	
	real(dp), dimension(nSecY,nSecX,Na) :: pres_Wag
	
	real(dp) :: deltaS

	real(dp), dimension(nSecY, nSecX, Nts) :: distrF02
	real(dp), dimension(nSecY, nSecX, Nts) :: distrF03
	real(dp), dimension(Mmax, Nmax, Nts) :: exterF0
!!======================================================================================================================	
contains
	subroutine AeroDyn(str_deflection, omega, tStep, totalF, exterF2, exterF3, L_LES, L_Wag)
		implicit none
		real(dp), dimension(Mmax+1,Nmax+1,Nts,3), intent(in) :: str_deflection
		real(dp), intent(in) :: omega
		real(dp), intent(in) :: tStep

		real(dp), dimension(Mmax,Nmax,Nts), intent(out) :: exterF2
		real(dp), dimension(Mmax,Nmax,Nts), intent(out) :: exterF3
		real(dp), dimension(1,Na), intent(out) :: totalF
		real(dp), dimension(1,Na), intent(out) :: L_LES
		real(dp), dimension(1,Na), intent(out) :: L_Wag
		
		real(dp), dimension(1,Nmax+1) :: aeNodeX
		real(dp), dimension(1,Mmax+1) :: aeNodeY
		
		real(dp), dimension(Mmax+1,Nmax+1,Nts) :: aeNodeZ
		real(dp), dimension(Mmax+1,Nmax+1,Na) :: zLoc
		
		integer :: i, j, k, tt
		
		!! set chord length
		do j = 1,Nmax
			chordLen(j) = Mv(j)*b
			bChord(j)   = 0.5*Mv(j)*b
		end do 
		!! set initial values
		
		!! control point location
		do j = 1, Nmax
			aeNodeX(1,j) = elemLenX *(j-1)
		end do
		do j = 1, Mmax
			aeNodeY(1,j) = elemLenY *(j-1)
		end do
		aeNodeZ = 0._dp
		zLoc    = 0._dp
!		print *, aeNodeX
		!!==============================================================================================================
		!! Kinematics
		do tt = 1, Na
			alpha10(1,tt)   = pi/2._dp
			alpha_dot(1,tt) = 0._dp
			beta10(1,tt)    = (5.5_dp/18._dp) *pi *sin(omega *(tt-1) *tStep -pi/2._dp) + (5.5_dp/18._dp) *pi
			beta_dot(1,tt)  = (5.5_dp/18._dp) *pi *omega *cos(omega *(tt-1) *tStep -pi/2._dp)
		end do
		
		do tt = 1, Na
			if (tt <= Na/2) then
				AoA(:,tt)     = alpha10(1,tt)
				AoA_dot(:,tt) = alpha_dot(1,tt)
			else
				AoA(:,tt)     = alpha10(1,tt)
				AoA_dot(:,tt) = -alpha_dot(1,tt)
			end if
		end do
!		print *, 'AoA', AoA
		!===============================================================================================================
		!! update new angles
		do j = 1, Nmax
			strX(1,j) = (j-1) *a
		end do
		do k = 1, Mmax
			strY(1,j) = (j-1) *b
		end do

!		do i = 1, Nts-1
!			call biInterp(strX, strY, str_deflection(:,:,i), aeNodeX, aeNodeY, aeNodeZ(:,:,i))
!		end do
		aeNodeZ = str_deflection(:,:,:,3)
!!!!!		
		do i = 1, Nts
			zLoc(:,:,i) = aeNodeZ(:,:,i)
		end do

!		print *, 'zLoc', zLoc
		do tt = 1, Na
			do j = 1, Nmax
				inAoA(j,tt)  = asin((-zLoc(1,j,tt) +zLoc(Mv(j),j,tt)) /chordLen(j))
				
				if (tt <= (Na+1)/2) then                         !! need to make sure the deflection and angle change in corresponding direction
				    updAoA(j,tt) = AoA(1,tt) -inAoA(j,tt)
				else
				    updAoA(j,tt) = AoA(1,tt) +inAoA(j,tt)
				end if
			end do
		end do

		updAoA_dot(:,1) = AoA_dot(:,1)
		do tt = 2, Na
			updAoA_dot(:,tt) = (updAoA(:,tt) -updAoA(:,tt-1)) /tStep
		end do
				
		AoA     = updAoA
		AoA_dot = updAoA_dot
		
!		do tt = 1, Nts-1
!			AoA_str(j, tt) = AoA(j, 2*tt-1)
!			AoA_str(j, tt) = 0.5_dp * (AoA(j,tt) + AoA(j,tt+1))
!		end do
!		AoA_str(j, Nts) = AoA(j, Nts)
!!!!!
!		print *, 'AoA'
!		print *, AoA(Ma-1,:)*180._dp/pi
!		print *, '-----------------------------------------------------------------------------------------------------'
!		print *, AoA(1,:)*180._dp/pi
		!!==============================================================================================================
		!! Leading Edge Suction Force
		do j = 1, nSecX
			elemLocX(1,j) = elemLenX/2._dp + elemLenX *(j-1)
		end do
		pres_LES = 0._dp
		
		do tt = 1, Na
			do j = 1, nSecX
				vPol(j,tt) = beta_dot(1,tt) * elemLocX(1,j)
				
				if (AoA(j,tt) <= pi/2) then
					AA0(1,j) = sin(AoA(j,tt)) + 0.5_dp * chordLen(j) * AoA_dot(j,tt) /vPol(j,tt)
					AA1(1,j) = 0.5_dp * chordLen(j) * AoA_dot(j,tt) * cos(AoA(j,tt)) /vPol(j,tt)
				else
					AA0(1,j) = sin(AoA(j,tt)) - 0.5_dp * chordLen(j) * AoA_dot(j,tt)  /vPol(j,tt)
					AA1(1,j) = -0.5_dp * chordLen(j) * AoA_dot(j,tt) * cos(AoA(j,tt)) /vPol(j,tt)
				end if

				do k = 1, Mv(j)
					theta(1,k) = 1._dp*(Mv(j)+1-k)/(Mv(j)+1) * pi
					gamma_theta(k,j,tt) = 2 * vPol(j,tt) * ( AA0(1,j) * (1+cos(theta(1,k))) /sin(theta(1,k)) + AA1(1,j) *sin(theta(1,k)) )
					
					pres_LES(k,j,tt) = rho_air * vPol(j,tt) * gamma_theta(k,j,tt)
				end do
			end do
		end do

		do tt = 1, Na
				L_LES(1,tt) = 0._dp
				do j = 1, nSecX
					do k = 1, Mv(j)
						distrF_LES(k,j,tt) = pres_LES(k,j,tt) * elemLenX * elemLenY * sin(AoA(j,tt)) * sin(AoA(j,tt))
						L_LES(1,tt)        = L_LES(1,tt) + distrF_LES(k,j,tt)
					end do
				end do
				if (tt > (Na-1)/2) then
					L_LES(1,tt) = L_LES(1,tt-(Na-1)/2)
				end if
		end do
		!!==============================================================================================================
		!! Apparent Mass
		L_app = 0._dp
		pres_app = 0._dp
!		do tt = 1, Na
!			do j = 1, nSecX
!				app_v(j,tt) = beta_dot(1,tt) * elemLocX(1,j)
! 
!				if (AoA(j,tt) <= pi/2) then
!					t_A0(j,tt) = sin(AoA(j,tt)) + 0.5_dp * chordLen(j) * AoA_dot(j,tt) * cos(AoA(j,tt)) / app_v(j,tt)
!					t_A1(j,tt) = 0.5_dp * chordLen(j) * AoA_dot(j,tt) * cos(AoA(j,tt)) / app_v(j,tt)
!				else
!					t_A0(j,tt) = sin(AoA(j,tt)) - 0.5_dp * chordLen(j) * AoA_dot(j,tt) * cos(AoA(j,tt)) / app_v(j,tt)
!					t_A1(j,tt) = -0.5_dp * chordLen(j) * AoA_dot(j,tt) * cos(AoA(j,tt)) / app_v(j,tt)
!				end if
!			end do
!		end do
		
!		do j = 1, nSecX
!			do tt = 1, Na-1
!				dt_A0V(j,tt) = (t_A0(j,tt+1)*app_v(j,tt+1) - t_A0(j,tt)*app_v(j,tt) ) /tStep
!				dt_A1V(j,tt) = (t_A1(j,tt+1)*app_v(j,tt+1) - t_A1(j,tt)*app_v(j,tt) ) /tStep
!			end do
!			dt_A0V(j,Na) = dt_A0V(j,Na-1)
!			dt_A1V(j,Na) = dt_A1V(j,Na-1)
!		end do
		
!		do tt = 1, Na
!			do j = 1, nSecX
!				L_app(1,tt) = L_app(1,tt) + rho_air * chordLen(j)**2 * pi * (1.5_dp * dt_A0V(j,tt) + 0.5_dp * dt_A1V(j,tt)) &
!																	 * cos(AoA(j,tt)) * elemLenX / 2._dp
!				
!				do k = 1, nSecY
!					theta(1,k)     = 1._dp*k/(1._dp*nSecY) * pi
!					gammapp(k,j,tt) = 2 * vPol(1,j) * ( t_A0(1,j) * (1+cos(theta(1,k))) / sin(theta(1,k)) + t_A1(1,j) &
!												   * sin(theta(1,k)) ) * chordLen(j)
!					pres_app(k,j,tt) = rho_air * beta_dot(1,tt) * elemLocX(1,j) * gammapp(k,j,tt)
!				end do
!			end do
!		end do
		
		!!==============================================================================================================
		!! Wagner Function 
		do j = 1, Nmax
			Xpart(j,1)        = 0._dp
			Ypart(j,1)        = 0._dp
			Xpart(j,(Na+1)/2) = 0._dp
			Ypart(j,(Na+1)/2) = 0._dp
		end do
		
		pres_Wag   = 0._dp
		Wag_Lift   = 0._dp
		L_Wag      = 0._dp
		
		do tt = 2, Na-1
			if (tt<=(Na-1)/2) then
				do j = 2, nSecX
					V_circ(1,j) = beta_dot(1,tt) *elemLocX(1,j)
					deltaS      = V_circ(1,j) *tStep /bChord(j)
					Xpart(j,tt)   = Xpart(j,tt-1) *exp(-0.041*deltaS) + 0.165 *(AoA(j,tt+1) - AoA(j,tt-1)) /2
					Ypart(j,tt)   = Ypart(j,tt-1) *exp(-0.32*deltaS)  + 0.335 *(AoA(j,tt+1) - AoA(j,tt-1)) /2
					AoA_eqv(j,tt) = AoA(j,tt) -Xpart(j,tt) -Ypart(j,tt)

					cl_Wag(j,tt)     = 2 *pi *sin(AoA_eqv(j,tt)) *cos(AoA_eqv(j,tt))
					pres_Wag(:,j,tt) = 0.5 *rho_air *V_circ(1,j)**2. *cl_Wag(j,tt)
					
					Wag_Lift(j,tt)   = 0.5 *rho_air *V_circ(1,j)**2. *cl_Wag(j,tt) *chordLen(j) *elemLenX
					L_Wag(1,tt)      = L_Wag(1,tt) + Wag_Lift(j,tt)
				end do
			else
				pres_Wag(:,j,tt) = pres_Wag(:,j, tt-(Na-1)/2)
				
				L_Wag(1,tt)      = L_Wag(1,tt-(Na-1)/2)
			end if
		end do
		
		!! Kusser Function 
!		do j = 1, Nmax
!			Xpart(j,1)        = 0._dp
!			Ypart(j,1)        = 0._dp
!			Xpart(j,(Na+1)/2) = 0._dp
!			Ypart(j,(Na+1)/2) = 0._dp
!		end do
		
!		pres_Wag   = 0._dp
!		Wag_Lift   = 0._dp
!		L_Wag      = 0._dp
		
!		do tt = 2, Na-1
!			if (tt<=(Na-1)/2) then
!				do j = 2, nSecX
!					V_circ(1,j)  = beta_dot(1,tt) * elemLocX(1,j)
!					deltaS       = V_circ(1,j) * tStep / bChord
!					Xpart(j,tt)   = Xpart(j,tt-1) * exp(-0.041*deltaS) + 0.165 * (AoA(j,tt+1) - AoA(j,tt-1)) / 2
!					Ypart(j,tt)   = Ypart(j,tt-1) * exp(-0.32*deltaS)  + 0.335 * (AoA(j,tt+1) - AoA(j,tt-1)) / 2
!					AoA_eqv(j,tt) = AoA(j,tt) - Xpart(j,tt) - Ypart(j,tt)

!					cl_Wag(j,tt)     = 2 * pi * sin(AoA_eqv(j,tt)) * cos(AoA_eqv(j,tt))
!					pres_Wag(:,j,tt) = (1./2.) * rho_air * V_circ(1,j)**2. * cl_Wag(j,tt)
					
!					Wag_Lift(j,tt)   = (1./2.) * rho_air * V_circ(1,j)**2. * cl_Wag(j,tt) * chordLen(j) * elemLenX
!					L_Wag(1,tt)      = L_Wag(1,tt) + Wag_Lift(j,tt)
!				end do
!			else
!				pres_Wag(:,j,tt) = pres_Wag(:,j, tt-(Na-1)/2)
			
!				L_Wag(1,tt)      = L_Wag(1,tt-(Na-1)/2)
!			end if
!		end do
	
	!!==================================================================================================================
		totalF = L_LES + L_Wag !+ L_app

		!! interpolation from aerodynamic mesh to structure mesh
		!!==============================================================================================================
		do tt = 1, Nts
			do j = 1, Nmax
				do k = 1, Mv(j)
					distrF02(k,j,tt) = (pres_LES(k,j,tt) +pres_Wag(k,j,tt)) * cos(AoA(j,tt))
					distrF03(k,j,tt) = (pres_LES(k,j,tt) +pres_Wag(k,j,tt)) * sin(AoA(j,tt))
				end do
			end do
		end do
!		do tt = 1, Nts-1
!			do j = 1, nSecX
!				do k = 1, nSecY
!				    exterF0(k,  j,  tt) = exterF0(k,  j,  tt) + distrF0(k,j,tt) / 4._dp
!				    exterF0(k+1,j,  tt) = exterF0(k+1,j,  tt) + distrF0(k,j,tt) / 4._dp
!				    exterF0(k,  j+1,tt) = exterF0(k,  j+1,tt) + distrF0(k,j,tt) / 4._dp
!				    exterF0(k+1,j+1,tt) = exterF0(k+1,j+1,tt) + distrF0(k,j,tt) / 4._dp
!				end do
!			end do
!		end do
		! eleStrX, eleStrY is the structure element center
		do j = 1, Nmax-1
			eleStrX(1,j) = 0.5_dp * a + (j-1) * a
		end do
		do k = 1, Mmax-1
			eleStrY(1,j) = 0.5_dp * b + (j-1) * b
		end do
!		do tt = 1, Nts-1
!			 call biInterp( aeNodeX, aeNodeY, exterF0(:,:,tt), eleStrX, eleStrY, exterF(:,:,tt))
!		end do
		exterF2 = distrF02
		exterF3 = distrF03

	end subroutine AeroDyn


end module Aero_Dyn






