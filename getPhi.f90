!! Xuan 02/29/2016
module Phi_mod
	use types
	use GaussQuad_mod
	implicit none
contains	
	subroutine getQ(E, nu, Q1_mat)
		implicit none
		real(dp), intent(in) :: E, nu
		real(dp), dimension(5, 5), intent(out) :: Q1_mat

		real(dp) :: Gg
		real(dp) :: c1, c2
		
		Q1_mat = 0._dp

		Gg = E/(2._dp*(1._dp+nu))

		Q1_mat(1,1) = E /(1._dp-nu*nu)
		Q1_mat(2,2) = Q1_mat(1,1)
		Q1_mat(1,2) = nu *Q1_mat(1,1)
		Q1_mat(2,1) = Q1_mat(1,2)

		c1 = 5.0_dp /6.0_dp
		c2 = 5.0_dp /6.0_dp                              !! c1, c2 is the same?? //

		Q1_mat(4,4) = Gg *c1
		Q1_mat(5,5) = Gg *c2
		Q1_mat(3,3) = Gg
	end subroutine getQ
!-----------------------------------------------------------------------------------------------------------------------
	subroutine getS(s, Kcap_1, Kcap_2, E, nu, matS)           !! need to change one of the s to another name
		implicit none
		real(dp), intent(in) :: E, nu
		real(dp), intent(in) :: s
		real(dp), dimension(3,3), intent(in) :: Kcap_1, Kcap_2
		real(dp), dimension(12,12), intent(out) :: matS
			
		real(dp) :: k_1, k_2, k_4, k_5, k_61, k_62, k_6
		real(dp), dimension(60) :: tempS
		real(dp), dimension(5,12) :: S1
		real(dp), dimension(5,5) :: Q1_mat
		real(dp), dimension(12,5) :: MatTemp

		k_5  = Kcap_1(1,2)
		k_1  =-Kcap_1(1,3)
		k_61 =-Kcap_1(2,3)
		k_4  = Kcap_2(1,2)
		k_62 =-Kcap_2(1,3)
		k_2  =-Kcap_2(2,3)
		k_6  = k_61 + k_62

		tempS = (/1._dp,  0._dp,  0._dp,  0._dp,       0._dp,   0._dp, 1._dp, 0._dp, 0._dp,   0._dp,   &
				  0._dp,  0._dp,  1._dp,  0._dp,       0._dp,   s,     0._dp, 0._dp, 0._dp,   0._dp,   &
				  0._dp,  s,      0._dp,  0._dp,       0._dp,   0._dp, 0._dp, s,     0._dp,   0._dp,   &
				  0._dp,  0._dp,  s,      0._dp,       0._dp,   0._dp, s,     0._dp, 0._dp,   0._dp,   &
				  s,      0._dp,  0._dp,  0._dp,       0._dp,   0._dp, 0._dp, s,     0._dp,   0._dp,   &
				  -k_5*s, 0._dp,  -k_4*s, 1._dp-k_2*s, -k_61*s, 0._dp, k_4*s, k_5*s, -k_62*s, 1._dp-k_1*s /)
		S1 = reshape(tempS,(/5, 12/))
		call getQ(E, nu, Q1_mat)
	
		MatTemp = matmul(transpose(S1), Q1_mat)
		matS = matmul(MatTemp, S1)
	end subroutine getS
!----------------------------------------------------------------------------------------------------------------------
	subroutine getPhi(E, nu, thick, Kcap_1, Kcap_2, phi)
		implicit none
		real(dp), intent(in) :: E, nu, thick
		real(dp), dimension(3,3), intent(in) :: Kcap_1, Kcap_2
		real(dp), dimension(12, 12), intent(out) :: phi

		integer, parameter :: NN = 5
		real(dp), dimension(12, 12) :: matS
		real(dp), dimension(NN) :: si, wsi                              !! the output of GaussQuad
		integer :: ii
		call GaussQuad(NN, -thick/2, thick/2, si, wsi)
		
		phi = 0._dp
		do ii = 1, NN
			call getS(si(ii), Kcap_1, Kcap_2, E, nu, matS)
			phi = phi + matS *wsi(ii)                             !! S_mat is 12*12, Wsi(sii) is one value ?? //
		end do
	end subroutine getPhi
	
end module Phi_mod
