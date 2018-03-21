!! 03/06/2016 Xuan Yang
! module of quadrature points
module GaussQuad_mod
	use types
	use quad_int
	implicit none
	
contains
	subroutine GaussQuad(gp, a0, a, x1, Wx1)
		implicit none
	
		integer, intent(in) :: gp
		real(dp), intent(in) :: a0, a
	
		real(dp), dimension(gp), intent(out) :: x1, Wx1
	
		integer :: i
	
		call sub_quad_int(gp)
	
		if (a0==0) then
			do i = 1, gp
				x1(i)  = (quad_pt(i)+1._dp) * a / 2._dp
				Wx1(i) = weight(i) * a / 2._dp
			end do
		elseif (a0==-a) then
			do i = 1, gp
				x1(i)  = quad_pt(i) * a
				Wx1(i) = weight(i) * a
			end do
		else
			print *, "unexpected request for GaussQuad"
		end if
	
	end subroutine GaussQuad

end module GaussQuad_mod
