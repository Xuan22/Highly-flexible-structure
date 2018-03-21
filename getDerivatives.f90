!! Xuan 02/29/2016
module Derivative
	use types
	use str_constants
	use shapeFunction
	use initCurvature
	use global_var, only : ele, counter
	use modfcn
	implicit none
contains
	subroutine getInpDeriv(xxi, yyj, dx, dy, T_01, q2ele, cc_1, cc_2, cc_3, cc_4)
		implicit none
		real(dp), intent(in) :: xxi, yyj
		real(dp), intent(in) :: dx, dy
	
		real(dp), dimension(3, 3), intent(in)   :: T_01
		real(dp), dimension(1, d*4), intent(in) :: q2ele
	
		real(dp), intent(out) :: cc_1, cc_2, cc_3, cc_4
	
		!! call shapefunction2
		real(dp), dimension(5,56)  :: Nsd
		real(dp), dimension(24,56) :: Dsd
		real(dp), dimension(24,1)  :: Ud
	
		!! call get strain
		real(dp) :: e1dx, e2dx

		!! call gettransformationmatrix
		real(dp), dimension(3,3) :: T1dx
		real(dp) :: gammadx
	
		real(dp) :: kgdx, gdx
		real(dp) :: gamma_1ddx, gamma_2ddx

		call shapeFunc2(xxi+dx, yyj+dy, Nsd, Dsd)                                                 !! Nsd could be Nsd(:,:)
	
		Ud = matmul(Dsd, transpose(q2ele))                                                                     !! here Ud is just one 24*1 array?? //

		call initCurv(T_01, Kcap_10, Kcap_20)

		call getStrain(Ud, e1dx, e2dx) 
		
		call getTransform(Ud, e1dx, e2dx, T1dx, gammadx) !! here Ud is just one 24*1 array?? eldx is 1 value?? //

		kgdx = (1+e2dx) / (1+e1dx)

		nf    = 2
		xf(1) = 0.01_dp
		xf(2) = 0.01_dp
		tol   = 1.0e-6
		lwa   = 20		                
						
		call hybrd1(fcn, nf, xf, fvec, tol, info, wa, lwa)            ! gdx  = fsolve(@(g) myfu(g,kgdx,gammadx),0)    !! this is the same as called from main program

		gamma_2ddx = xf(2)
		gamma_1ddx = gammadx - gamma_2ddx

		cc_1 = cos(gamma_1ddx) /(cos(gammadx) *(1.+e2dx))
		cc_2 = sin(gamma_2ddx) /(cos(gammadx) *(1.+e1dx))
		cc_3 = sin(gamma_1ddx) /(cos(gammadx) *(1.+e2dx))
		cc_4 = cos(gamma_2ddx) /(cos(gammadx) *(1.+e1dx))

	end subroutine getInpDeriv
	!! -----------------------------------------------------------------------------------------------------------------
	subroutine getDerivative(xxi, yyj, q2ele, T_0, C1sx, C1sy, C2sx, C2sy, C3sx, C3sy, C4sx, C4sy)
		implicit none
		real(dp), intent(in) :: xxi, yyj
		real(dp), dimension(3, 3), intent(in)   :: T_0
		real(dp), dimension(1, d*4), intent(in) :: q2ele               !! Is q2ele dimension right? row //
	
		real(dp), intent(out) :: C1sx, C1sy, C2sx, C2sy
		real(dp), intent(out) :: C3sx, C3sy, C4sx, C4sy

		real(dp) :: c1pdx, c2pdx, c3pdx, c4pdx                         !! it's real number right??
		real(dp) :: c1mdx, c2mdx, c3mdx, c4mdx
		real(dp) :: c1pdy, c2pdy, c3pdy, c4pdy
		real(dp) :: c1mdy, c2mdy, c3mdy, c4mdy

		real(dp) :: dx, dy

		dx = a /100.0_dp
		dy = 0.0_dp
		call getInpDeriv(xxi, yyj, dx, dy, T_0, q2ele, c1pdx, c2pdx, c3pdx, c4pdx)

		dx = - a /100.0_dp
		dy = 0.0_dp
		call getInpDeriv(xxi, yyj, dx, dy, T_0, q2ele, c1mdx, c2mdx, c3mdx, c4mdx)

		C1sx = (c1pdx-c1mdx) /(2.*dx)
		C2sx = (c2pdx-c2mdx) /(2.*dx)
		C3sx = (c3pdx-c3mdx) /(2.*dx)
		C4sx = (c4pdx-c4mdx) /(2.*dx)

		dx = 0.0_dp
		dy = b /100.0_dp
		call getInpDeriv(xxi, yyj, dx, dy, T_0, q2ele, c1pdy, c2pdy, c3pdy, c4pdy)

		dx = 0.0_dp
		dy = - b /100.0_dp
		call getInpDeriv(xxi, yyj, dx, dy, T_0, q2ele, c1mdy, c2mdy, c3mdy, c4mdy)

		C1sy = (c1pdy-c1mdy) / (2._dp*dy)
		C2sy = (c2pdy-c2mdy) / (2._dp*dy)
		C3sy = (c3pdy-c3mdy) / (2._dp*dy)
		C4sy = (c4pdy-c4mdy) / (2._dp*dy)

	end subroutine getDerivative
	
end module Derivative


