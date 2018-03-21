!! 03/04/2016 Xuan
module ele_index
	use types
	use str_constants
	use global_var
	implicit none
	integer, dimension(eleN, 1)  :: nel
	integer, dimension(eleN, d*4) :: index1
	integer :: num1
contains
	subroutine sub_ele_index()
		implicit none
		do j = 1,Mmax
			call element_index(j, coeff, Nx)
			do i = 1,Nx
				ele = 0
				do ni = 1,Mmax-1
					ele = ele + Nh(ni) *coeff(ni)
				end do
				ele = ele + i
				nel(ele,1) = ele
				
				do iter = 1,2
			    	do k = 1,14
			    		num1 = 0
			    		do ni = 1,j
			    	    	num1 = num1 +(Nh(ni)+1) *d *coeff(ni)
			    	    end do
			    	    index1(ele, k+d*(iter-1)) = num1 +(i-1+iter-1)*d +k
			   		end do
				end do
				do iter = 3,4
			    	do k = 1,14
			    		num1 = 0
			    		do ni = 1,j
			    	    	num1 = num1 + (Nh(ni)+1) *d *coeff(ni)
			    	    end do
			    	    index1(ele, k+d*(iter-1)) = num1 +(Nx+1)*d +(i-1+iter-3)*d +k
			    	end do
				end do
			end do !end i loop
		end do     !end j loop
	end subroutine sub_ele_index
!!--------------------------------------------------------------------------------------------
	subroutine element_index(jj, coeff, Nx)
		implicit none
		integer, intent(in) :: jj
		integer, dimension(Mmax), intent(out) :: coeff
		integer, intent(out) :: Nx
		integer :: ic, ik
		do ic = 1,Mmax
			coeff(ic) = 0
		end do
		if(jj==1) then
			Nx = Nh(1)
			do ik = 1, Mmax
				coeff(ik) = 0
			end do
		else
			do ik = 2,jj
				coeff(ik-1) = 1
			end do
			Nx = Nh(jj)
		endif
	end subroutine element_index
!!--------------------------------------------------------------------------------------------
	subroutine node_index(jj, nCoeff, Nx)
		implicit none
		integer, intent(in) :: jj
		integer, dimension(Mmax), intent(out) :: nCoeff
		integer, intent(out) :: Nx
		integer :: ic, ik
		do ic = 1,Mmax
			nCoeff(ic) = 0
		end do
		if(jj==1) then
			Nx = Nh(1)
		else if(jj==Mmax+1) then
			do ik = 2,jj
				nCoeff(ik-1) = 1
			end do
			Nx = Nh(Mmax)
		else
			do ik = 2,jj
				nCoeff(ik-1) = 1
			end do
			Nx = Nh(jj)
		endif
	end subroutine node_index
!---------------------------------------------------------------------------------------------
end module ele_index


