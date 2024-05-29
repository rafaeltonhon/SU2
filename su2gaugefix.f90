module su2gaugefix
    use su2lattice
    use omp_lib
    implicit none
    contains
    ! subroutine that makes a random gauge transform
    subroutine random_gauge(u,ugauge)
        real(kind=r8), dimension(nr,nr,nr,nr,4,4) :: u, ugauge
        real(kind=r8), dimension(nr,nr,nr,nr,4) :: omega
        real(kind=r8) :: r(4)
        integer(kind=i4) :: e1, e2, e3, e4, mi, n(4)
        ! create the random matrix
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            call random_number(r)
            omega(e1,e2,e3,e4,:)=r/dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3)+r(4)*r(4))
        enddo
        enddo
        enddo
        enddo

        ! do the gauge transformation
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            do mi=1,4
                n=(/e1,e2,e3,e4/)
                n(mi)=n(mi)+1
                if(n(mi).gt.nr) n(mi)=1
                ugauge(e1,e2,e3,e4,mi,:)=linkmult(omega(e1,e2,e3,e4,:),&
                linkmult(u(e1,e2,e3,e4,mi,:),dlink(omega(n(1),n(2),n(3),n(4),:))))
            enddo
        enddo
        enddo
        enddo
        enddo
    endsubroutine random_gauge

    ! subroutine that fixes the direct maximal 
    subroutine maximal_center_gauge(u,ug,tol)
        real(kind=r8), dimension(nr,nr,nr,nr,4,4) :: u, ug
        real(kind=r8), dimension(nr,nr,nr,nr,4) :: g
        real(kind=r8) :: d(4,8), bigd(4,4), x(4,1), norm, tracenew, traceold, ngs, omega, tol
        integer(kind=i4) :: e1, e2, e3, e4, mi, n(4)
        integer(kind=i4) :: i, j, l, igs, cont
        ! initialize the gauge fixed configuration
        ug=u
        traceold=0.0_r8
        g=1.0_r8
        omega=1.70_r8
        cont=0
        do igs=1,100000
            ! visit each lattice site
            !$OMP PARALLEL DO PRIVATE(e1,e2,e3,e4,mi,l,d,i,j,x,norm,n) SHARED(ug,g)
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                ! we must pick the values of d(i,l)
                do l=1,4
                    n=(/e1,e2,e3,e4/)
                    d(1,l)=-ug(n(1),n(2),n(3),n(4),l,1)
                    d(2,l)=-ug(n(1),n(2),n(3),n(4),l,2)
                    d(3,l)=-ug(n(1),n(2),n(3),n(4),l,3)
                    d(4,l)=ug(n(1),n(2),n(3),n(4),l,4)
                    n(l)=n(l)-1
                    if(n(l).lt.1) n(l)=nr
                    d(1,l+4)=ug(n(1),n(2),n(3),n(4),l,1)
                    d(2,l+4)=ug(n(1),n(2),n(3),n(4),l,2)
                    d(3,l+4)=ug(n(1),n(2),n(3),n(4),l,3)
                    d(4,l+4)=ug(n(1),n(2),n(3),n(4),l,4)
                enddo

                ! now we can find the matrix bigd, tha must be diagonalized
                do i=1,4
                do j=1,4
                    bigd(i,j)=0.0_r8
                    do l=1,8
                    bigd(i,j)=bigd(i,j)+d(i,l)*d(j,l)
                    enddo
                enddo
                enddo

                ! now we must find the eigenvector associated with the bigest eigenvalue of this matrix
                ! we will use the power method
                x=1.0_r8
                do i=1,50
                    x=matmul(bigd,x)
                    ! we must find the bigger component of x
                    norm=x(1,1)
                    do j=2,4
                        if(x(j,1).gt.norm) norm=x(j,1)
                    enddo
                    x=x/norm
                enddo

                ! we found the gauge transformation
                g(e1,e2,e3,e4,1)=x(1,1)
                g(e1,e2,e3,e4,2)=x(2,1)
                g(e1,e2,e3,e4,3)=x(3,1)
                g(e1,e2,e3,e4,4)=x(4,1)

                ! overrelaxation step
                g(e1,e2,e3,e4,:)=g(e1,e2,e3,e4,:)-ident()
                g(e1,e2,e3,e4,:)=g(e1,e2,e3,e4,:)*omega
                g(e1,e2,e3,e4,:)=g(e1,e2,e3,e4,:)+ident()
                norm=dsqrt(detlink(g(e1,e2,e3,e4,:)))
                g(e1,e2,e3,e4,:)=g(e1,e2,e3,e4,:)/norm
                !g(e1,e2,e3,e4,1)=g(e1,e2,e3,e4,1)/norm
                !g(e1,e2,e3,e4,2)=g(e1,e2,e3,e4,2)/norm
                !g(e1,e2,e3,e4,3)=g(e1,e2,e3,e4,3)/norm
                !g(e1,e2,e3,e4,4)=g(e1,e2,e3,e4,4)/norm
                
                do mi=1,4
                n=(/e1,e2,e3,e4/)
                n(mi)=n(mi)-1
                if(n(mi).lt.1) n(mi)=nr
                ug(e1,e2,e3,e4,mi,:)=linkmult(g(e1,e2,e3,e4,:),ug(e1,e2,e3,e4,mi,:))
                ug(n(1),n(2),n(3),n(4),mi,:)=linkmult(ug(n(1),n(2),n(3),n(4),mi,:),dlink(g(e1,e2,e3,e4,:)))
                enddo
            enddo
            enddo
            enddo
            enddo
            !$OMP END PARALLEL DO
            
            ! we gauge transform the configuration
            tracenew=0.0_r8
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                do mi=1,4
                    ! compute the trace of the transformed configuration
                    tracenew=tracenew+ug(e1,e2,e3,e4,mi,4)**2.0_r8/4.0_r8
                enddo
            enddo
            enddo
            enddo
            enddo
            if(abs(tracenew-traceold).lt.tol) cont=cont+1
            if(abs(tracenew-traceold).gt.tol) cont=0
            if(cont.eq.50)then
                print*,'Maximal Center Gauge fixed, tolerance=',abs(tracenew-traceold)
                exit
            endif
            traceold=tracenew
        enddo
        if(cont.ne.50) print*,'MCG NOT fixed, tolerance=',abs(tracenew-traceold)
    endsubroutine maximal_center_gauge

    ! subroutine that fixes the  landau gauge
    subroutine landau_gauge(u,ug,tol)
        real(kind=r8) :: u(nr,nr,nr,nr,4,4), ug(nr,nr,nr,nr,4,4)
        real(kind=r8) :: g(nr,nr,nr,nr,4), a(4), lambda, nabla, sum, tol
        integer(kind=i4) :: e1, e2, e3, e4, mi, ac, int, cont, n(4), nm(4)
        ug=u
        cont=0
        do int=1,100000
            !$OMP PARALLEL DO PRIVATE(e1,e2,e3,e4,mi,a,n,nm,lambda) SHARED(ug,g)
            ! sweep the lattice
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                ! find the marix g on this site
                a=0.0_r8
                do mi=1,4
                    n=(/e1,e2,e3,e4/)
                    nm=n
                    nm(mi)=nm(mi)-1
                    if(nm(mi).lt.1) nm(mi)=nr

                    a(1)=a(1)+ug(nm(1),nm(2),nm(3),nm(4),mi,1)-ug(n(1),n(2),n(3),n(4),mi,1)
                    a(2)=a(2)+ug(nm(1),nm(2),nm(3),nm(4),mi,2)-ug(n(1),n(2),n(3),n(4),mi,2)
                    a(3)=a(3)+ug(nm(1),nm(2),nm(3),nm(4),mi,3)-ug(n(1),n(2),n(3),n(4),mi,3)
                    a(4)=a(4)+ug(nm(1),nm(2),nm(3),nm(4),mi,4)+ug(n(1),n(2),n(3),n(4),mi,4)
                enddo
                lambda=a(1)**2.0_r8+a(2)**2.0_r8+a(3)**2.0_r8+a(4)**2.0_r8
                lambda=-0.5_r8*dsqrt(lambda)
                g(e1,e2,e3,e4,:)=0.5_r8*a/lambda
                g(e1,e2,e3,e4,:)=g(e1,e2,e3,e4,:)/dsqrt(detlink(g(e1,e2,e3,e4,:)))

                ! gauge transform
                do mi=1,4
                    n=(/e1,e2,e3,e4/)
                    nm=n
                    nm(mi)=nm(mi)-1
                    if(nm(mi).lt.1) nm(mi)=nr
                    ug(n(1),n(2),n(3),n(4),mi,:)=linkmult(g(e1,e2,e3,e4,:),ug(n(1),n(2),n(3),n(4),mi,:))
                    ug(nm(1),nm(2),nm(3),nm(4),mi,:)=linkmult(ug(nm(1),nm(2),nm(3),nm(4),mi,:),dlink(g(e1,e2,e3,e4,:)))
                enddo
            enddo
            enddo
            enddo
            enddo
            !$OMP END PARALLEL DO
            ! compute the convergence
            sum=0.0_r8
            !$OMP PARALLEL DO PRIVATE(e1,e2,e3,e4,mi,ac,n,nm,nabla) SHARED(sum)
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                ! compute the divergence of each color
                do ac=1,3
                    nabla=0.0_r8
                    do mi=1,4
                        n=(/e1,e2,e3,e4/)
                        nm=n
                        nm(mi)=nm(mi)-1
                        if(nm(mi).lt.1) nm(mi)=nr

                        nabla=nabla+ug(n(1),n(2),n(3),n(4),mi,ac)-ug(nm(1),nm(2),nm(3),nm(4),mi,ac)
                    enddo
                    sum=sum+nabla**2.0_r8
                enddo
            enddo
            enddo
            enddo
            enddo
            !$OMP END PARALLEL DO
            sum=sum/(1.0_r8*nr**4.0_r8)
            !print*,sum
            if(sum.lt.tol) cont=cont+1
            if(sum.gt.tol) cont=0
            if(cont.eq.50)then
                print*,'Landau Gauge Convergence, tolerance=', sum
                exit
            endif
        enddo
        if(cont.ne.50) print*,'Landau Gauge DO NOT CONVERGED, tol= ',sum
    endsubroutine landau_gauge

    ! subroutine tht fixes the maximal abelian gauge
    subroutine maximal_abelian_gauge(u,ug,tol)
        real(kind=r8) :: u(nr,nr,nr,nr,4,4), ug(nr,nr,nr,nr,4,4), g(nr,nr,nr,nr,4)
        real(kind=r8) :: mag, magold, tol, a, b, c, aux(4)
        integer(kind=i4) :: e1, e2, e3, e4, mi, int, cont, n(4)
        ug=u
        magold=0.0_r8
        g=1.0_r8
        do int=1,10000
            ! apply the gauge trasformation
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                ! compute the coefficientes
                a=0.0_r8
                b=0.0_r8
                c=0.0_r8
                do mi=1,4
                    n=(/e1,e2,e3,e4/)
                    aux=ug(n(1),n(2),n(3),n(4),mi,:)
                    a=a+aux(1)*aux(4)+aux(2)*aux(3)
                    b=b+aux(2)*aux(4)-aux(1)*aux(3)
                    c=c-aux(3)*aux(3)-aux(4)*aux(4)+0.50_r8

                    n(mi)=n(mi)-1
                    if(n(mi).lt.1) n(mi)=nr
                    aux=dlink(ug(n(1),n(2),n(3),n(4),mi,:))
                    a=a+aux(1)*aux(4)+aux(2)*aux(3)
                    b=b+aux(2)*aux(4)-aux(1)*aux(3)
                    c=c-aux(3)*aux(3)-aux(4)*aux(4)+0.50_r8
                enddo
                a=2.0_r8*a
                b=2.0_r8*b
                c=2.0_r8*c

                ! the eigenvector of the biggest eigenvalue
                g(e1,e2,e3,e4,1)=-a*c+a*dsqrt(a*a+b*b+c*c)
                g(e1,e2,e3,e4,2)=-b*c+b*dsqrt(a*a+b*b+c*c)
                g(e1,e2,e3,e4,3)=0.0_r8
                g(e1,e2,e3,e4,4)=a*a+b*b
                g(e1,e2,e3,e4,:)=g(e1,e2,e3,e4,:)/dsqrt(detlink(g(e1,e2,e3,e4,:)))

                ! gauge transform the configuration
                do mi=1,4
                    n=(/e1,e2,e3,e4/)
                    ug(n(1),n(2),n(3),n(4),mi,:)=linkmult(g(e1,e2,e3,e4,:),ug(n(1),n(2),n(3),n(4),mi,:))
                    n(mi)=n(mi)-1
                    if(n(mi).lt.1) n(mi)=nr
                    ug(n(1),n(2),n(3),n(4),mi,:)=linkmult(ug(n(1),n(2),n(3),n(4),mi,:),dlink(g(e1,e2,e3,e4,:)))
                enddo
            enddo
            enddo
            enddo
            enddo

            ! verify the convergence
            mag=0.0_r8
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                do mi=1,4
                    mag=mag+ug(e1,e2,e3,e4,mi,1)**2.0_r8+ug(e1,e2,e3,e4,mi,2)**2.0_r8
                enddo
            enddo
            enddo
            enddo
            enddo
            mag=mag/(8.0_r8*nr**4.0_r8) ! if the out diagonal mean vanishs, the field is abelian!?
            if(abs(mag).lt.tol) cont=cont+1
            if(abs(mag).gt.tol) cont=0
            if(cont.eq.50) then
                print*,'Maximal Abelian Gauge convergence, tolerance= ',mag
                exit
            endif
        enddo
        if(cont.ne.50) print*,'Maximal Abelian Gauge DO NOT CONVERGED, tolerance= ',mag
    endsubroutine maximal_abelian_gauge
endmodule su2gaugefix