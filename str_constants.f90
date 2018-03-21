! 10/07/2017 Xuan
module str_constants
	use types
	implicit none
	!! math constants
	real(dp),parameter :: pi = 3.1415926_dp
	
	integer,parameter :: gp = 5                                      ! Gauss points-- gp
	integer,parameter :: d  = 14                                     ! d = degree of freedom per node. Dont change d.
	
	real(dp),parameter :: ll = 0.12_dp
	real(dp),parameter :: bb = 0.05_dp                                ! Dimensions in meters

	real(dp),parameter :: alpha = 0.5_dp                             ! Constants for Newmark Beta method of nonlinear analysis
	real(dp),parameter :: beta  = 0.25_dp *(alpha+0.5_dp)**2	      !! what is beta???? //
	
	!! index
	integer,parameter :: Mmax = 8
	integer,parameter :: Nmax = 7
!	integer :: Nh(Mmax) = (/1,2,4,4,5,5/)
	integer :: Nh(Mmax) = (/1,2,4,5,6,6,7,7/)
!	integer :: Nh(Mmax) = (/1,3,4,6,6,7,7,8,8/)
!	integer,parameter :: eleN = 21
!	integer,parameter :: nodN = 33
	integer,parameter :: eleN = 38
	integer,parameter :: nodN = 54

	integer,dimension(Mmax) :: coeff
	integer :: ni
	integer :: Nx
	integer :: nNode
	integer :: ele
	integer :: indexNum

	!! simulation time
	real(dp),parameter :: t_init  = 0._dp
	real(dp),parameter :: t_final = 1.0_dp /10.0_dp

	! kinematics
	integer,parameter :: Nts = 20    				! number of time intervals
	real(dp),parameter :: dt = 0.1_dp      			! time step used in Newmark-beta method

	! element property
	real(dp),parameter :: a = ll /Nmax                          ! x direction element length
	real(dp),parameter :: b = bb /Mmax    				   	     ! y direction ''
!===============================================================================================
	! Aerodynamic parameters
	real(dp),parameter :: transFr_PI = 57.3_dp
	real(dp),parameter :: transTo_PI = 0.0175_dp
	
	integer,parameter :: Nstr = Nts
	integer,parameter :: Na   = Nts
	
	integer,parameter :: Ma = Nmax
	integer,parameter :: Ka = Mmax
	integer,parameter :: nSecX = Ma
	integer,parameter :: nSecY = Ka
	
	real(dp),parameter :: rho_air = 1.225_dp
	real(dp),parameter :: wingLen = ll
	
!	integer :: Mv(Nmax) = (/6,5,4,4,2/)
	integer :: Mv(Nmax) = (/8,7,6,6,5,4,2/)
!	integer :: Mv(Nmax) = (/9,8,8,7,6,6,4,2/)
	
	real(dp),parameter :: elemLenX = a
	real(dp),parameter :: elemLenY = b
	
	real(dp),dimension(Nmax) :: chordLen
	real(dp),dimension(Nmax) :: bChord
!===============================================================================================
	! Initial Curvature
	real(dp) :: k_50  = 0._dp
	real(dp) :: k_10  = 0._dp
	real(dp) :: k_610 = 0._dp
	real(dp) :: k_40  = 0._dp
	real(dp) :: k_620 = 0._dp
	real(dp) :: k_20  = 0._dp
	real(dp) :: k_60  = 0._dp
	!! Material properties
    real(dp),parameter :: Em1    = 5._dp *10**8
	real(dp),parameter :: nu1    = 0.3_dp
	real(dp),parameter :: thick1 = 0.0005
	real(dp),parameter :: rho1   = 1._dp *10**3
	
	real(dp),parameter :: Em2    = 5._dp *10**9
	real(dp),parameter :: nu2    = 0.3_dp
	real(dp),parameter :: thick2 = 0.0005
	real(dp),parameter :: rho2   = 1._dp *10**3
	
	real(dp) :: E
	real(dp) :: nu
	real(dp) :: rho
	real(dp) :: thick
    real(dp) :: dm
	
	!! GaussQuad
	real(dp) :: a0 = 0._dp
	
	!! identity matrix
	real(dp),dimension(3,3) :: idenMat3 = reshape((/1._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,1._dp/),(/3,3/))
	real(dp),dimension(3,3) :: T_0 = reshape((/1._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,1._dp/),(/3,3/))
	!! dimension of fixed boundary condition
	integer,parameter :: ADOF = 14
	integer,parameter :: BDOF = 3
end module str_constants




