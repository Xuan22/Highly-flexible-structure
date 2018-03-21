module modfcn
	use types
	use global_var
	use strain
	use getTransformation
	implicit none
contains
	subroutine fcn(nf,x,fvec,iflag)
		implicit none
		    integer :: nf
		    integer :: iflag
		    double precision :: x(nf), fvec(nf)
		    
		    real(dp) :: e1_vir, e2_vir
		
			call getStrain(U_mat, e1_vir, e2_vir)
			
			call getTransform(U_mat, e1_vir, e2_vir, T1, gamma1)

			fvec(1) = x(1) + x(2) - gamma1
			fvec(2) = (1._dp+e1_vir) * sin(x(1)) - (1._dp+e2_vir) * sin(x(2))

	end subroutine fcn

end module modfcn
