! 03/14/2016
module addboundary
	use types
	use str_constants
	use global_var, only : A1_mat
	use ele_index
	implicit none
contains
	subroutine addboundfixpin(dQ_mat, A2_mat, Qr_mat)
		implicit none
		real(dp), dimension(NEQ, 1), intent(in) :: dQ_mat
		integer, dimension(BDOF), intent(in) :: A2_mat
		real(dp), dimension(nodN*d, 1), intent(out) :: Qr_mat
		
		real(dp), dimension(NEQ, 1) :: qr
		real(dp), dimension(nodN*d, 1) :: q1
		
		integer :: ii, ip, jp, ak, ite
		integer :: dofNum
		integer :: dofPrev
		
		qr = dQ_mat
		q1 = 0._dp

		call node_index(Mmax+1, coeff, Nx)
!		print*, coeff
		dofPrev = 0
		do ip = 1, Mmax
			dofPrev = dofPrev + (Nh(ip)+1) *d *coeff(ip)
		end do
		do ii = 1, dofPrev
			q1(ii, 1) = qr(ii,1)
		end do
		do ii = dofPrev+1, dofPrev+2*d
			q1(ii, 1) = 0._dp
		end do
		dofPrev = dofPrev +2*d
!		print*, dofPrev
		do i = 1, (Nmax-1)
			do ii = dofPrev+1, dofPrev+d
				if (ii==dofPrev+A2_mat(1)) then
					q1(ii,1) = 0._dp
				elseif (ii>dofPrev+A2_mat(1).and.ii<dofPrev+A2_mat(2)) then
					q1(ii,1) = qr(ii-2*ADOF-(i-1)*BDOF-1, 1)
				elseif (ii==dofPrev+A2_mat(2)) then
					q1(ii,1) = 0._dp
				elseif (ii>dofPrev+A2_mat(2).and.ii<dofPrev+A2_mat(3)) then
					q1(ii,1) = qr(ii-2*ADOF-(i-1)*BDOF-2, 1)
				elseif (ii==dofPrev+A2_mat(3)) then
					q1(ii,1) = 0._dp
				else
					q1(ii, 1) = qr(ii-2*ADOF-(i-1)*BDOF-3, 1)
				endif
			end do
			dofPrev = dofPrev + d
		end do
		Qr_mat = q1
	end subroutine addboundfixpin
end module addboundary

