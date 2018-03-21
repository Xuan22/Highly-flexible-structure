!! Xuan 02/29/2016
module massM_mod
	use types
	use GaussQuad_mod
	implicit none

contains

	subroutine getMreal(xx, thick, rho, M_mat)
		implicit none

		real(dp), intent(in) :: xx, thick, rho
		real(dp), dimension(18,18), intent(out) :: M_mat
	
		real(dp), dimension(6,6) :: I_mat
		real(dp) :: g14, g25
	
		real(dp) :: I_0,  I_1,  I_2,  I_5,  I_6,  I_7, I_8
		real(dp) :: I_51, I_61, I_71, I_81
		real(dp) :: I_55, I_56, I_57, I_58
		real(dp) :: I_66, I_67, I_68, I_77, I_78, I_88
	
		g14 = 0._dp
		g25 = 0._dp
	
		I_0 = rho
		I_1 = rho * xx
		I_2 = rho * xx**2
		I_5 = rho * g14
		I_6 = rho * (xx - 4._dp*xx**3/(3._dp*thick**2))				! deleted unnecessary ()
		I_7 = rho * (xx - 4._dp*xx**3/(3._dp*thick**2))
		I_8 = rho * g25

		I_51 = rho * g14 * xx
		I_61 = rho * xx * (xx - 4._dp*xx**3/(3._dp*thick**2))
		I_71 = rho * xx * (xx - 4._dp*xx**3/(3._dp*thick**2))
		I_81 = rho * g25 * xx

		I_55 = rho * g14 * g14 * xx**0                                              !! x^0 ??? //
		I_56 = rho * g14 * (xx - 4._dp*xx**3/(3*thick**2))
		I_57 = rho * g14 * (xx - 4._dp*xx**3/(3*thick**2))
		I_58 = rho * g14 *g25

		I_66 = rho * (xx - 4._dp*xx**3/(3._dp*thick**2)) * (xx - 4._dp*xx**3/(3._dp*thick**2))      !! why use .* ?
		I_67 = rho * (xx - 4._dp*xx**3/(3._dp*thick**2)) * (xx - 4._dp*xx**3/(3._dp*thick**2))
		I_68 = rho * g25 * (xx - 4._dp*xx**3/(3._dp*thick**2))
		I_77 = rho * (xx - 4._dp*xx**3/(3._dp*thick**2)) * (xx - 4._dp*xx**3/(3._dp*thick**2))
		I_78 = rho * g25 * (xx - 4._dp*xx**3/(3._dp*thick**2))
		I_88 = rho * g25 * g25 * xx**0


		I_mat = reshape((/I_0, I_1,  I_5,  I_6,  I_7,  I_8, &
			   			  I_1, I_2,  I_51, I_61, I_71, I_81, &
			 		  	  I_5, I_51, I_55, I_56, I_57, I_58, &
			   			  I_6, I_61, I_56, I_66, I_67, I_68, &
			   			  I_7, I_71, I_57, I_67, I_77, I_78, &
			   			  I_8, I_81, I_58, I_68, I_78, I_88/), (/6, 6/))
	
   M_mat = reshape( (/I_0,  I_1,  I_5,  I_6,  I_7,  I_8, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
		  	   I_1,  I_2,  I_51, I_61, I_71, I_81, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
	 		   I_5,  I_51, I_55, I_56, I_57, I_58, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
			   I_6,  I_61, I_56, I_66, I_67, I_68, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
		  	   I_7,  I_71, I_57, I_67, I_77, I_78, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
		 	   I_8,  I_81, I_58, I_68, I_78, I_88, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
	 		   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_0,  I_1,  I_5,  I_6,  I_7,  I_8,  0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
	 		   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_1,  I_2,  I_51, I_61, I_71, I_81, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
	  		   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_5,  I_51, I_55, I_56, I_57, I_58, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
	 		   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_6,  I_61, I_56, I_66, I_67, I_68, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
			   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_7,  I_71, I_57, I_67, I_77, I_78, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&    
			   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_8,  I_81, I_58, I_68, I_78, I_88, 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
			   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_0,  I_1,  I_5,  I_6,  I_7,  I_8,  &
			   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_1,  I_2,  I_51, I_61, I_71, I_81, &
			   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_5,  I_51, I_55, I_56, I_57, I_58, &
	 		   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_6,  I_61, I_56, I_66, I_67, I_68, &
	 		   0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_7,  I_71, I_57, I_67, I_77, I_78, &
	0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,I_8,  I_81, I_58, I_68, I_78, I_88 /), (/18, 18/) )

	end subroutine getMreal
!-----------------------------------------------------------------------------------------------------------------------
	subroutine getMassM(thick,rho,MassM)
		implicit none

		real(dp), intent(in) :: thick, rho
	
		real(dp), dimension(18,18), intent(out) :: MassM

		integer, parameter :: NN = 5
		integer :: ii
	
		real(dp), dimension(18,18) :: temp
		real(dp), dimension(18,18) :: M_mat
		real(dp), dimension(NN) :: si, Wsi

		temp = 0._dp

		call GaussQuad(NN, -thick/2._dp, thick/2._dp, si, Wsi)

		do ii = 1, size(si)
		
			call getMreal(si(ii), thick, rho, M_mat)
		
			temp = temp + M_mat * Wsi(ii)
		end do

		MassM = temp

	end subroutine getMassM
!-----------------------------------------------------------------------------------------------------------------------	
	
end module massM_mod


