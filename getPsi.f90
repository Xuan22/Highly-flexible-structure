!! Xuan 02/29/2016
module Psi_mod
	use types
	use str_constants, only: k_10, k_20, k_60
	implicit none
contains
	subroutine getPsi(e1, e2, gamma_1d, gamma_2d, U_mat, Kcap_1, Kcap_2, psi)
		implicit none
		real(dp), intent(in) :: e1, e2, gamma_1d, gamma_2d
		real(dp), dimension(3,3), intent(in) :: Kcap_1, Kcap_2
		real(dp), dimension(24,1), intent(in) :: U_mat
		real(dp), dimension(12,1), intent(out) :: psi
		
		real(dp), dimension(12) :: psi1
		real(dp) :: k_1, k_2, k_61, k_62, k_6
		real(dp) :: term1, term2, term3, term4, term5, term6
		
		k_1  =-Kcap_1(1,3)
		k_61 =-Kcap_1(2,3)
		k_62 =-Kcap_2(1,3)
		k_2  =-Kcap_2(2,3)
		k_6  = k_61 + k_62
		
		term1 = (1._dp+e1) *cos(gamma_1d) - 1._dp
		term2 = (1._dp+e2) *cos(gamma_2d) - 1._dp
		term3 = (1._dp+e1) *sin(gamma_1d) + (1._dp+e2) *sin(gamma_2d)
		term4 = k_1 -k_10
		term5 = k_2 -k_20
		term6 = k_6 -k_60
		
		psi1 = (/term1, term2, term3, term4, term5, term6, U_mat(20,1), U_mat(21,1), U_mat(23,1), U_mat(24,1), U_mat(19,1), U_mat(22,1)/)
	 	psi = reshape(psi1, (/12,1/))
	end subroutine getPsi
end module Psi_mod
