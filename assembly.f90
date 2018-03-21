!! Xuan 03/11/2016
module assembly_mod
	use types
	use str_constants
	implicit none
contains
	subroutine assembly(kk, k_e, index_ele, kk_upd)                !! input is kk, output is kk?? //
		implicit none
		real(dp), dimension(nodN*d, nodN*d), intent(in) :: kk
		real(dp), dimension(d*4, d*4), intent(in) :: k_e
		integer, dimension(1, d*4), intent(in) :: index_ele
		
		real(dp), dimension(nodN*d, nodN*d), intent(out) :: kk_upd
		integer :: si, sj, edof
		integer :: ii, jj
		
		edof = size(index_ele)
		
		do si = 1, edof                                            ! loop for i-th element
			ii = index_ele(1, si)                                  !! what is trying to do here?
			do sj = 1, edof
				
				jj = index_ele(1, sj)
				kk_upd(ii,jj) = kk(ii,jj) + k_e(si,sj)              !! Is this right??
			end do
		end do
	
	end subroutine assembly
!-----------------------------------------------------------------------------------------------------------------------
	subroutine assemblyr(ff, f_e, index_ele, ff_upd)               !! ff is input, and output // update ff value
		implicit none
		real(dp), dimension(nodN*d, 1), intent(in) :: ff
		real(dp), dimension(d*4, 1), intent(in) :: f_e
		integer, dimension(d*4, 1), intent(in) :: index_ele
		real(dp), dimension(nodN*d, 1), intent(out) :: ff_upd
	
		integer :: ai
		integer :: ii
		integer :: edof

		edof = size(index_ele)
		do ai = 1, edof
			ii = index_ele(ai,1)
			ff_upd(ii,1) = ff(ii,1)+f_e(ai,1)
		end do
	end subroutine assemblyr

end module assembly_mod




