module ChangeInT
	use types
	use str_constants, only: a, b, gp, d, idenMat3
	use shapeFunction
	use initCurvature
	use Strain
	use global_var, only: ele, counter
	use modfcn
	use getTransformation
	use array_inv
	implicit none
contains
	subroutine getChangeInT(xi, yj, q2ele, T_00, dx, dy, T1dx)
		implicit none
		real(dp), intent(in) :: dx, dy
		real(dp), intent(in) :: xi, yj
		
		real(dp), dimension(1, d*4), intent(in) :: q2ele
		real(dp), dimension(3, 3), intent(in) :: T_00
		real(dp), dimension(3, 3), intent(out) :: T1dx

		real(dp), dimension(5, d*4)  :: Ns_vir
		real(dp), dimension(24, 1)  :: Uds
		real(dp), dimension(24, d*4) :: Dssd
	
		real(dp), dimension(3,3) :: Kcap_10_vir, Kcap_20_vir
		real(dp), dimension(3,3) :: T1sdx
		real(dp), dimension(3,3) :: Gamsdx, Gammas_capdx
		real(dp) :: e1sdx, e2sdx, kgsdx
		real(dp) :: gammasdx
		real(dp) :: gamma_1sdx, gamma_2sdx
	
		call shapeFunc2(xi+dx, yj+dy, Ns_vir, Dssd)       !! Dssd change from (:,:,ele) to (:,:)

		Uds = matmul(Dssd, transpose(q2ele))                  !! Dssd is 24 * 56 ??//

		call initCurv(T_00, Kcap_10_vir, Kcap_20_vir)

		call getStrain(Uds, e1sdx, e2sdx)                     !! in getstrain, e1, e2 is real number //

		call getTransform(Uds, e1sdx, e2sdx, T1sdx, gammasdx)

		nf    = 2
		xf(1) = 0._dp
		xf(2) = 0._dp
		tol   = 1.0e-10
		lwa   = 30		                
		call hybrd1(fcn, nf, xf, fvec, tol, info, wa, lwa)


		gamma_2sdx = xf(2)
		gamma_1sdx = gammasdx - gamma_2sdx
	
		Gamsdx = reshape((/cos(gamma_1sdx), sin(gamma_2sdx), 0._dp, &
	   			 		   sin(gamma_1sdx), cos(gamma_2sdx), 0._dp, &
						   0._dp,           0._dp,           1._dp/), (/3,3/) )                           !! Is this matrix form right ???? //
						   
	    call inv22(Gamsdx,idenMat3,Gammas_capdx)         	!	call inverse(Gamsdx, Gammas_capdx, 3)	  !! How to do this??//
	    
		T1dx = matmul(Gammas_capdx, T1sdx)
	
	end subroutine getChangeInT

end module ChangeInT


