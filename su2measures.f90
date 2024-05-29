module su2measures
    use su2lattice
    use omp_lib
    implicit none
    contains
    function wilsonplanar(u,n,mi,ni,a,b)
        real(kind=r8) :: u(nr,nr,nr,nr,4,4), w(4), sum, wilsonplanar
        integer(kind=i4) :: mi, ni, n(4), a, b, i
        sum=0.0_r8
        w=ident()
        do i=1,a
            w=linkmult(w,u(n(1),n(2),n(3),n(4),mi,:))
            !if((a.eq.nr).and.(b.eq.nr)) print*,u(n(1),n(2),n(3),n(4),mi,:)
            !if((a.eq.nr).and.(b.eq.nr)) print*,n
            !if((a.eq.nr).and.(b.eq.nr)) print*,w
            n(mi)=n(mi)+1
            if(n(mi).gt.nr) n(mi)=1
        enddo
                
        do i=1,b
            w=linkmult(w,u(n(1),n(2),n(3),n(4),ni,:))
            !if((a.eq.nr).and.(b.eq.nr)) print*,u(n(1),n(2),n(3),n(4),ni,:)
            !if((a.eq.nr).and.(b.eq.nr)) print*,n
            !if((a.eq.nr).and.(b.eq.nr)) print*,w
            n(ni)=n(ni)+1
            if(n(ni).gt.nr) n(ni)=1
        enddo
        do i=1,a
            n(mi)=n(mi)-1
            if(n(mi).lt.1) n(mi)=nr
            w=linkmult(w,dlink(u(n(1),n(2),n(3),n(4),mi,:)))
            !if((a.eq.nr).and.(b.eq.nr)) print*,u(n(1),n(2),n(3),n(4),mi,:)
            !if((a.eq.nr).and.(b.eq.nr)) print*,n
            !if((a.eq.nr).and.(b.eq.nr)) print*,w
        enddo
        do i=1,b
            n(ni)=n(ni)-1
            if(n(ni).lt.1) n(ni)=nr
            w=linkmult(w,dlink(u(n(1),n(2),n(3),n(4),ni,:)))
            !if((a.eq.nr).and.(b.eq.nr)) print*,u(n(1),n(2),n(3),n(4),ni,:)
            !if((a.eq.nr).and.(b.eq.nr)) print*,n
            !if((a.eq.nr).and.(b.eq.nr)) print*,w
        enddo
        wilsonplanar=w(4)
    endfunction wilsonplanar

    ! function that makes the measure of the wilson loop in the lattice
    subroutine measurewilson(u,a,b,w)
        real(kind=r8), intent(in) :: u(nr,nr,nr,nt,4,4)
        real(kind=r8), intent(out) :: w(a,b)
        real(kind=r8) :: sum
        integer(kind=i4) :: a,b, i, j, n(4), e1, e2, e3, e4, mu, nu
        w=0.0_r8
        ! we compute the loops w(a,b) and save the result
        !$OMP PARALLEL DO PRIVATE(e1,e2,e3,e4,mu,nu,i,j) SHARED(w)
        do i=1,a
        do j=1,b
            !sum=0.0_r8
            ! compute the mean wilson loop
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                n=(/e1,e2,e3,e4/)
                do mu=1,3
                do nu=mu+1,4
                    w(i,j)=w(i,j)+wilsonplanar(u,n,mu,nu,i,j)
                    ! the lattice is simmetric, so we
                    ! must consider the case a->b and b->a
                    ! now we interchange a and b
                    if(i.ne.j)then 
                        w(i,j)=w(i,j)+wilsonplanar(u,n,mu,nu,j,i)                    
                    endif
                    ! the loop for this site is calculated!
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
            
        enddo
        enddo
        !$OMP END PARALLEL DO
        do i=1,a
        do j=1,b
            if(i.eq.j)then
                w(i,j)=w(i,j)/(6.0_r8*nt*nr**3)    ! quadratic loops
            else
                w(i,j)=w(i,j)/(12_r8*nt*nr**3) ! rectangular loops
            endif
        enddo
        enddo
    endsubroutine measurewilson

    !+=============================================================================!+
    ! subrotuine that computes the gluon propagator
    ! compute the propagator
    subroutine gluonpropagator(ugauge,flag)
        real(kind=r8) :: ugauge(nr,nr,nr,nr,4,4)
        real(kind=r8) :: gluonp(0:nr/2), p(0:nr/2), lambda, fouriercos, fouriersin
        integer(kind=i4) :: x1, x2, x3, x4, mi, ac, k, flag
        gluonp=0.0_r8
        lambda=nr**4.0_r8
        ! for zero momentum
        do ac=1,3
        do mi=1,4
            fouriercos=0.0_r8
            do x1=1,nr
            do x2=1,nr
            do x3=1,nr
            do x4=1,nr
                fouriercos=fouriercos+ugauge(x1,x2,x3,x4,mi,ac)
            enddo
            enddo
            enddo
            enddo
            gluonp(0)=gluonp(0)+fouriercos**2.0_r8
        enddo
        enddo
        gluonp(0)=gluonp(0)/(12.0_r8*lambda)

        ! for non zero momentum
        do k=1,nr/2
            do ac=1,3
            do mi=1,4
            fouriercos=0.0_r8
            fouriersin=0.0_r8
                do x1=1,nr
                do x2=1,nr
                do x3=1,nr
                do x4=1,nr
                    !print*,x1,x2,x3,x4,mi,ac
                    fouriercos=fouriercos+dcos(k*x4*2.0_r8*pi/nr)*ugauge(x1,x2,x3,x4,mi,ac)
                    fouriersin=fouriersin+dsin(k*x4*2.0_r8*pi/nr)*ugauge(x1,x2,x3,x4,mi,ac)
                enddo
                enddo
                enddo
                enddo
                gluonp(k)=gluonp(k)+fouriercos**2.0_r8+fouriersin**2.0_r8
            enddo
            enddo
        enddo

        do k=0,nr/2
            p(k)=2.0_r8*dsin(pi*k/nr)
        enddo
        gluonp(1:nr/2)=gluonp(1:nr/2)/(9.0_r8*lambda)
        write(flag+1,*) gluonp(:)
        write(10,*) p(:)
        do k=0,nr/2
            write(flag+2,*) p(k),gluonp(k)
        enddo
    endsubroutine gluonpropagator

    ! gluon propagator
    subroutine gluon_formfactor(ug,flag)
        real(kind=r8) :: ug(nr,nr,nr,nr,4,4)
        real(kind=r8) :: f(0:nr/2), p(0:nr/2), lambda, deltatug, fouriercos, fouriersin
        integer(kind=i4) :: x1, x2, x3, x4, mi, ac, k, flag, n(4)
        f=0.0_r8
        lambda=nr**3.0_r8
        ! for zero momentum
        do ac=1,3
        do mi=1,4
            fouriercos=0.0_r8
            do x1=1,nr
            do x2=1,nr
            do x3=1,nr
            do x4=1,nr
                n=(/x1,x2,x3,nr/)
                n(4)=n(4)+1
                if(n(4).gt.nr) n(4)=1
                deltatug=ug(n(1),n(2),n(3),n(4),mi,ac)-ug(x1,x2,x3,x4,mi,ac)
                fouriercos=fouriercos+deltatug
            enddo
            enddo
            enddo
            enddo
            f(0)=f(0)+fouriercos**2.0_r8
        enddo
        enddo
        f(0)=f(0)/(9.0_r8*lambda)

        ! for non zero momentum
        do k=1,nr/2
            do ac=1,3
            do mi=1,4
                fouriercos=0.0_r8
                fouriersin=0.0_r8
                do x1=1,nr
                do x2=1,nr
                do x3=1,nr
                do x4=1,nr
                    n=(/x1,x2,x3,x4/)
                    n(4)=n(4)+1
                    if(n(4).gt.nr) n(4)=1
                    deltatug=ug(n(1),n(2),n(3),n(4),mi,ac)-ug(x1,x2,x3,x4,mi,ac)
                    fouriercos=fouriercos+dcos(k*x4*2.0_r8*pi/nr)*deltatug
                    fouriersin=fouriersin+dsin(k*x4*2.0_r8*pi/nr)*deltatug
                enddo
                enddo
                enddo
                enddo
                f(k)=f(k)+fouriercos**2.0_r8+fouriersin**2.0_r8
            enddo
            enddo
        enddo

        do k=0,nr/2
            p(k)=2.0_r8*dsin(pi*k/nr)
        enddo
        f(1:nr/2)=f(1:nr/2)/(6.0_r8*lambda)
        write(flag+1,*) f(:)
        write(10,*) p(:)
        do k=0,nr/2
            write(flag+2,*) p(k),f(k)
        enddo
    endsubroutine gluon_formfactor
endmodule su2measures