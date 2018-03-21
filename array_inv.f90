module array_inv
    use types, only :dp
    implicit none
    
contains
    
    subroutine inv21 (p,q,r)
        implicit none
        real(dp), intent (in)  :: p(:,:), q(:)
        real(dp), intent (out) :: r( size(q) )
        
        real(dp) :: a( size(q), size(q) ), c( size(q) )
        
        integer :: lda, n,info ,ipvt( size(q) )
        
        a = p
        c = q
        
        lda = size(a, dim=1)
        n = size(q)
        
        call dgefa(a,lda,n,ipvt,info)
        call dgesl(a,lda,n,ipvt,c,0)
        
        r = c
        
    end subroutine inv21
    
    subroutine inv22 (p,q,r)
        implicit none
        real(dp), intent (in)  :: p(:,:), q(:,:)
        real(dp), intent (out) :: r( size(q,dim=1), size(q,dim=2) )
        
        real(dp), dimension ( size(q,dim=1) ) :: q1, r1
        
        integer :: lda, n,info ,ipvt( size(q,dim=1) ), jj
        
        lda = size(p, dim=1)
        n   = size(q, dim=1)
        
        do  jj = 1, size(q,dim=2)
            q1(:) = q(:,jj)
            
            call inv21(p, q1, r1)
            
            r(:,jj) = r1(:)
        end do
        
    end subroutine inv22 
       
end module array_inv
