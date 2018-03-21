! Xuan 02/28/2016
module Strain
	use types
	use str_constants
	implicit none
contains
	subroutine getStrain(U_mat, e1, e2) 
		implicit none
		real(dp), dimension(24,1), intent(in)  :: U_mat
		real(dp), intent(out) :: e1, e2

	e1 = sqrt((1._dp +U_mat(2,1) -U_mat(7,1)*k_50 +U_mat(13,1)*k_10)**2 + (U_mat(8,1)+U_mat(1,1)*k_50+U_mat(13,1)*k_610)**2+&
			(U_mat(14,1) -U_mat(1,1)*k_10 -U_mat(7,1)*k_610)**2) - 1._dp
				  
	e2 = sqrt((U_mat(3,1) -U_mat(7,1)*k_40 +U_mat(13,1)*k_620)**2 + (1._dp+U_mat(9,1)+U_mat(1,1)*k_40+U_mat(13,1)*k_20)**2 +&
		    (U_mat(15,1) -U_mat(1,1)*k_620 -U_mat(7,1)*k_20)**2) - 1._dp

	end subroutine getstrain

end module Strain
