!! Xuan 02/27/2016
module applyboundary
	use types
	use str_constants
	use ele_index
	implicit none
contains
!---------------------------------------------------------------------------------------
	subroutine applyboundfixpin(kk, A1_mat, A2_mat, kkc)
		implicit none
		real(dp), dimension(nodN*d, nodN*d), intent(in) :: kk
		integer, dimension(ADOF), intent(in) :: A1_mat
		integer, dimension(BDOF), intent(in) :: A2_mat
		real(dp), dimension(NEQ, NEQ), intent(out) :: kkc
		
		real(dp), dimension(nodN*d, nodN*d) :: kk_copy
		real(dp), dimension(NEQ, nodN*d) :: kkb
		integer, dimension(2*ADOF+(Nmax-1)*BDOF) :: collect
		integer :: ii, ai, aj, ak, ite
		integer :: dofPrev
		integer :: dofNum
		
		kk_copy = kk
		!-------------------------------------------------------------------------------
		! find the leading edge start point
		call node_index(Mmax+1, coeff, Nx)
		dofPrev = 0
		do aj = 1, Mmax
			dofPrev = dofPrev + (Nh(aj)+1) *d *coeff(aj)
		end do
!		print *, dofPrev
		! boundary control degree in leading edge first element
		ii = 0
		do ite = dofPrev+1, dofPrev+d
		do ak = 1, size(A1_mat)
			dofNum = dofPrev +A1_mat(ak)
			if (ite == dofNum) then
				ii = ii + 1
				collect(ii) = ite
			end if
		end do
		end do

		do ite = dofPrev+d+1, dofPrev+2*d
		do ak = 1, size(A1_mat)
			dofNum = dofPrev +d +A1_mat(ak)
			if (ite == dofNum) then
				ii = ii + 1
				collect(ii) = ite
			end if
		end do
		end do
		dofPrev = dofPrev + 2*d
		! BC degree for pin condition
		do i = 1, Nmax-1
			do ite = dofPrev+1, dofPrev+d
			do ak = 1, BDOF
				dofNum = dofPrev +A2_mat(ak)
				if (ite == dofNum) then
					ii = ii +1
					collect(ii) = ite
				end if
			end do
			end do
			dofPrev = dofPrev +d
		end do
!		print *, collect
		!--------------------------------------------------------------------------------
		! copy rows
		ii = 0
		do ite = 1, nodN*d
			if (ite == collect(ii+1)) then
				ii = ii + 1
				if (ii > 2*ADOF+(Nmax-1)*BDOF) then
					print *, 'error in applyboundary'
				end if
			else
				kkb(ite-ii,:) = kk_copy(ite,:)
			end if
		end do
		!---------------------------------------------------------------------------------------------------------------
		! copy collums
		ii = 0
		do ite = 1, nodN*d
			if (ite == collect(ii+1)) then
				ii = ii + 1
				if (ii > 2*ADOF+(Nmax-1)*BDOF) then
					print *, 'error in applyboundary'
				end if
			else
				kkc(:,ite-ii) = kkb(:,ite)
			end if
		end do

	end subroutine applyboundfixpin
!-----------------------------------------------------------------------------------------------------------------------
	subroutine applyboundfixpinr1(rr, A1_mat, A2_mat, rrb)
		implicit none
		integer, dimension(ADOF), intent(in) :: A1_mat
		real(dp), dimension(nodN*d,1), intent(in) :: rr
		integer, dimension(BDOF), intent(in) :: A2_mat
		real(dp), dimension(nodN*d-2*ADOF-(Nmax-1)*BDOF,1), intent(out) :: rrb
		
		real(dp), dimension(nodN*d,1) :: rr_copy
		
		integer, dimension(2*ADOF+(Nmax-1)*BDOF) :: collect
		integer :: ii, ai, aj, ak, ite
		integer :: dofPrev
		integer :: dofNum
		
		rr_copy = rr
		!---------------------------------------------------------------------------------------------------------------
		call node_index(Mmax+1, coeff, Nx)
		dofPrev = 0
		do aj = 1, Mmax
			dofPrev = dofPrev + (Nh(aj)+1) *d *coeff(aj)
		end do
!		print*, dofPrev
		! boundary control degree in leading edge first element
		ii = 0
		do ite = dofPrev+1, dofPrev+d
			do ak = 1, size(A1_mat)
				dofNum = dofPrev +A1_mat(ak)
				if (ite == dofNum) then
					ii = ii + 1
					collect(ii) = ite
				end if
			end do
		end do
		
		do ite = dofPrev+d+1, dofPrev+2*d
			do ak = 1, size(A1_mat)
				dofNum = dofPrev +d +A1_mat(ak)
				if (ite == dofNum) then
					ii = ii + 1
					collect(ii) = ite
				end if
			end do
		end do
		dofPrev = dofPrev + 2*d
		! BC degree for pin condition
		do i = 1, Nmax-1
			do ite = dofPrev+1, dofPrev+d
			do ak = 1, BDOF
				dofNum = dofPrev +A2_mat(ak)
				if (ite == dofNum) then
					ii = ii + 1
					collect(ii) = ite
				end if
			end do
			end do
			dofPrev = dofPrev + d
		end do
!		print *, collect
		!---------------------------------------------------------------------------------------------------------------
		ii = 0
		do ite = 1, nodN*d
			if (ite == collect(ii+1)) then
				ii = ii + 1
				if (ii > 2*ADOF+(Nmax-1)*BDOF) then
					print *, 'error in applyboundary'
				end if
			else
				rrb(ite-ii,1) = rr_copy(ite,1)
			end if
		end do
!		print *, ii
		
	end subroutine applyboundfixpinr1
	
end module applyboundary



