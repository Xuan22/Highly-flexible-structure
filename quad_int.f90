! module of quadrature points
!
module quad_int
	use types
    implicit none

    integer, parameter  :: nqp  = 11    						 ! no of quadrature point
    real(dp), dimension(nqp) :: quad_pt, weight
    
contains

    subroutine sub_quad_int (nqp1)
        implicit none
        integer, intent(in) :: nqp1
        if (nqp1 == 5) then
    
            quad_pt(1)  = -0.9061798459386640_dp
            quad_pt(2)  = -0.5384693101056831_dp
            quad_pt(3)  =  0.0000000000000000_dp
            quad_pt(4)  =  0.5384693101056831_dp
            quad_pt(5)  =  0.9061798459386640_dp
    
            weight(1)   = 0.2369268850561891_dp     ! original one
            weight(2)   = 0.4786286704993665_dp
            weight(3)   = 0.5688888888888889_dp
            weight(4)   = 0.4786286704993665_dp
            weight(5)   = 0.2369268850561891_dp
        
        else if (nqp1 == 11) then
    
            quad_pt(1)  = -0.9782286581460570_dp
            quad_pt(2)  = -0.8870625997680953_dp
            quad_pt(3)  = -0.7301520055740494_dp
            quad_pt(4)  = -0.5190961292068118_dp
            quad_pt(5)  = -0.2695431559523450_dp
            quad_pt(6)  =  0.0000000000000000_dp
            quad_pt(7)  =  0.2695431559523450_dp
            quad_pt(8)  =  0.5190961292068118_dp
            quad_pt(9)  =  0.7301520055740494_dp
            quad_pt(10) =  0.8870625997680953_dp
            quad_pt(11) =  0.9782286581460570_dp
    
            weight(1)   = 0.0556685671161737_dp
            weight(2)   = 0.1255803694649046_dp
            weight(3)   = 0.1862902109277343_dp
            weight(4)   = 0.2331937645919905_dp
            weight(5)   = 0.2628045445102467_dp
            weight(6)   = 0.2729250867779006_dp
            weight(7)   = 0.2628045445102467_dp
            weight(8)   = 0.2331937645919905_dp
            weight(9)   = 0.1862902109277343_dp
            weight(10)  = 0.1255803694649046_dp
            weight(11)  = 0.0556685671161737_dp 
            
        end if

    end subroutine sub_quad_int
    
    
end module quad_int


