module ChangeInCurv
	use types
	use str_constants, only : a, b, gp, d
	use global_var, only : ele, counter
	use ChangeInT
	implicit none
contains
	subroutine getChangeInCurv(xi, yj, T_0, q2ele, Tdx, Tdy)
		implicit none
		real(dp), dimension(3,3), intent(in) :: T_0
		real(dp), intent(in) :: xi, yj
		real(dp), dimension(1, d*4), intent(in) :: q2ele                !! row 
	
		real(dp), dimension(3,3), intent(out) :: Tdx, Tdy				!! Is this dimension right??
		
		real(dp) :: dx, dy
		real(dp), dimension(3,3) :: T1pdx, T1mdx, T1pdy, T1mdy

		dx = a/100._dp
		dy = 0._dp
		call getChangeInT(xi, yj, q2ele, T_0, dx, dy, T1pdx)

		dx = - a/100._dp
		dy = 0._dp
		call getChangeInT(xi, yj, q2ele, T_0, dx, dy, T1mdx)

		dx = 0._dp
		dy = b/100._dp
		call getChangeInT(xi, yj, q2ele, T_0, dx, dy, T1pdy)

		dx = 0._dp
		dy = -b/100._dp
		call getChangeInT(xi, yj, q2ele, T_0, dx, dy, T1mdy)

		Tdx = (T1pdx-T1mdx) /(2.*a/100.)
		Tdy = (T1pdy-T1mdy) /(2.*b/100.)

	end subroutine getChangeInCurv

end module ChangeInCurv
