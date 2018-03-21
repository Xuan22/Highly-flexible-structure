!-----------------------------------------------
! 03/08/2016 Xuan
! 10/07/2017 modified
!-----------------------------------------------
module global_var
	use types
	use str_constants
	implicit none
	integer :: i, j, ix, jy, k, iter, counter, tt, bi, JJ
	real(dp), dimension(gp) :: x, y, Wx, Wy
	
	real(dp) :: cc, cc1
	real(dp) :: r22
	real(dp) :: r33
	
	integer :: reminder
	
	real(dp) :: C1sx = 0._dp
	real(dp) :: C1sy = 0._dp
	real(dp) :: C2sx = 0._dp
	real(dp) :: C2sy = 0._dp
	real(dp) :: C3sx = 0._dp
	real(dp) :: C3sy = 0._dp
	real(dp) :: C4sx = 0._dp
	real(dp) :: C4sy = 0._dp
	
	! call shapeFunction2
	real(dp), dimension(5, d*4)  :: Ns                    !! Ns dimension, what is 5 for?? shapefunction, fixed
	real(dp), dimension(24, d*4) :: D_mat
	! call shapeFunctionM
	real(dp), dimension(5, d*4)  :: Nm                    !! dimensions
	real(dp), dimension(11, d*4) :: Dm_mat
	
	! call getTransformationMatrixMass
	real(dp), dimension(3,3) :: T1m
	real(dp) :: gammam
	real(dp) :: kg
	
	! call getMassM
	real(dp), dimension(18,18) :: MassM
	
	! update
	real(dp), dimension(56, 18) :: tmp1
	real(dp), dimension(18, 56) :: tmp2
	real(dp), dimension(56, 12) :: tmp3
	real(dp), dimension(12, 56) :: tmp4
	real(dp), dimension(56, 12) :: tmp5
	real(dp), dimension(12, 1)  :: tmp6
	real(dp), dimension(56, 1)  :: tmp7
	
	! call getTransformationMatrix
	real(dp), dimension(3,3) :: T1
	real(dp) :: gamma1
	
	! solve U_mat
	real(dp), dimension(1, d*4) :: temp_q2
	real(dp), dimension(24, 1)  :: U_mat
	
	real(dp) :: gam2d
	real(dp) :: gam1d
	real(dp) :: gam1dm = 0._dp
	real(dp) :: gam2dm = 0._dp
	
	! matInverse
	real(dp), dimension(3,3) :: Gam, Gamma0
	real(dp), dimension(3,3) :: Tr
	
	! hybrd nonlinear
	integer :: nf, info, lwa
	double precision :: xf(2), fvec(2), tol, wa(30)

	! getstrain
	real(dp) :: e1, e2
	! getChangeInCurvature
	real(dp), dimension(3,3) :: Tdx, Tdy

	! gatpsi
	real(dp), dimension(3,3) :: Kcap_1, Kcap_2
	real(dp), dimension(3,3) :: Kcap_10, Kcap_20
	real(dp) :: k_1, k_2, k_61, k_62, k_5, k_4, k_6
	real(dp), dimension(12,1) :: Psi
	! getphi
	real(dp), dimension(12,12) :: Phi

	! getZ, Zm
	real(dp), dimension(12,24) :: Z_mat = 0._dp
	real(dp), dimension(18,11) :: Zm    = 0._dp
	real(dp) :: C1, C2, C3, C4
	! getR
	real(dp), dimension(24,1) :: R1_cap
	
	! Condense global matrix and solve use MATINV
	integer, parameter :: NHBW= (Nmax+3)*d
	integer, parameter :: NBW = 2*NHBW
	integer, parameter :: NEQ = nodN*d-2*ADOF-(Nmax-1)*BDOF
	
	! Assembly
	real(dp), dimension(nodN*d, nodN*d) :: MT, MT_assem
	real(dp), dimension(nodN*d, nodN*d) :: KT, KT_assem
    real(dp), dimension(nodN*d, 1) :: RR_cap, RR_assem
    real(dp), dimension(nodN*d, 1) :: KQ_cap, KQ_assem
    
    !! Apply boundary conditions
	integer, dimension(ADOF) :: A1_mat
	integer, dimension(3)    :: A2_mat
	
	real(dp), dimension(NEQ, NEQ) :: newMT1, newKT1
	real(dp), dimension(NEQ, NEQ) :: KT2
	real(dp), dimension(NEQ, NBW) :: KTT							! Condensed matrix for banded global matrix
	
	real(dp), dimension(NEQ, 1) :: newRR1, newKQ1
	real(dp), dimension(NEQ, 1) :: q_dot1, q_ddot1
	real(dp), dimension(NEQ, 1) :: dQ_mat, dQ1, dQ_dot, dQ_ddot
	real(dp), dimension(NEQ, 1) :: RR2
	
	real(dp), dimension(nodN*d, 1) :: Qr_mat, Qr_dot, Qr_ddot
	
	real(dp), dimension(Mmax+1, Nmax+1, Nts+1) :: deform
end module global_var




