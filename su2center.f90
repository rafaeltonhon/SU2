module su2center
    use su2lattice
    use su2measures
    use omp_lib
    implicit none
    !real(kind=r8), dimension(:,:), allocatable :: we, wo
    contains
    ! subroutine that does the center projection
    subroutine centerprojection(u,z)
        real(kind=r8) :: u(nr,nr,nr,nt,4,4), z(nr,nr,nr,nt,4,4)
        integer(kind=i4) :: e1,e2,e3,e4,mi!, z(nr,nr,nr,nt,4)
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            do mi=1,4
                if(u(e1,e2,e3,e4,mi,4).gt.0.0_r8)then
                    z(e1,e2,e3,e4,mi,:)=(/0.0_r8,0.0_r8,0.0_r8,1.0_r8/)
                else
                    z(e1,e2,e3,e4,mi,:)=(/0.0_r8,0.0_r8,0.0_r8,-1.0_r8/)
                endif
            enddo
        enddo
        enddo
        enddo
        enddo
    endsubroutine centerprojection

    ! subroutine that remove the vortices
    subroutine centerremotion(u,ur)
        real(kind=r8) :: u(nr,nr,nr,nt,4,4),ur(nr,nr,nr,nt,4,4)
        integer(kind=i4) :: e1,e2,e3,e4,mi
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            do mi=1,4
                if(u(e1,e2,e3,e4,mi,4).gt.0.0_r8) ur(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)
                if(u(e1,e2,e3,e4,mi,4).lt.0.0_r8) ur(e1,e2,e3,e4,mi,:)=-u(e1,e2,e3,e4,mi,:)
            enddo
        enddo
        enddo
        enddo
        enddo
    endsubroutine centerremotion

    ! function that says if we have a vortex in a given plaquette
    function vortexinplaq(z,n,mi,ni)
        real(kind=r8) ::  vortexinplaq,z(nr,nr,nr,nt,4,4)
        integer(kind=i4) :: mi, ni
        integer(kind=i4), dimension(4) :: n,npz,nmz,nzp,nzm,npp,npm,nmp

        call neighborhood(n,npz,nmz,nzp,nzm,npp,npm,nmp,mi,ni)

        vortexinplaq=z(n(1),n(2),n(3),n(4),mi,4)*z(npz(1),npz(2),npz(3),npz(4),ni,4)*&
        &z(nzp(1),nzp(2),nzp(3),nzp(4),mi,4)*z(n(1),n(2),n(3),n(4),ni,4)
    endfunction vortexinplaq

    ! function that returns the number of vortices in a flat surface in the lattice
    function numberofvortices(z,n,mi,ni,lmi,lni)
        real(kind=r8) :: z(nr,nr,nr,nt,4,4), nv,  numberofvortices
        integer(kind=i4) ::  mi, ni, lmi, lni, ini, imi
        integer(kind=i4) :: n(4), m(4)
        nv=0
        m=n
        ! we need to check plaquete per plaquette of the lattice
        do imi=1,lmi
        ! move in the ni direction
        do ini=1,lni
            if(vortexinplaq(z,m,mi,ni).lt.0.0_r8) nv=nv+1.0_r8
            m(ni)=m(ni)+1
            if(m(ni).gt.nr) m(ni)=1
        enddo
        ! move in the mi direction for all the ni's
        m(mi)=m(mi)+1
        if(m(mi).gt.nr) m(mi)=1
        m(ni)=n(ni)
        enddo

        numberofvortices=nv
    endfunction numberofvortices
    
    subroutine pvortexdensity(z,a,b)
        real(kind=r8) :: z(nr,nr,nr,nt,4,4), probeven, probodd, nv
        integer(kind=i4) :: i, j, a, b, lamba, area
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        lamba=nr*nr*nr*nt
        do i=1,a
        do j=1,b
            probeven=0.0_r8
            probodd=0.0_r8
            area=i*j
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                do mi=1,3
                do ni=mi+1,4
                    ! count the number of vortices
                    nv=numberofvortices(z,(/e1,e2,e3,e4/),mi,ni,i,j)
                    ! nv is odd
                    if(mod(nv,2.0_r8).gt.0.0_r8)then
                        probodd=probodd+1.0_r8
                    else
                        probeven=probeven+1.0_r8
                    endif
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo

            ! we have the probability for this area, get out with it
            probeven=probeven/(6.0_r8*lamba)
            probodd=probodd/(6.0_r8*lamba)
            write(400+area,*) area, probeven, probodd, probeven+probodd
        enddo
        enddo
        
    endsubroutine pvortexdensity

    subroutine vortexsperatedwilson(u,z,a,b,we,wo)
        real(kind=r8), intent(in) :: u(nr,nr,nr,nt,4,4),z(nr,nr,nr,nt,4,4)
        real(kind=r8), intent(out) :: we(a,b), wo(a,b)
        real(kind=r8) :: sumo, sume
        integer(kind=i4) :: a,b, i, j, n(4), ne, no
        integer(kind=i4) :: e1, e2, e3, e4, mu, nu
        ! we compute the loops w(a,b) and save the result
        do i=1,a
        do j=1,b
            sumo=0.0_r8
            sume=0.0_r8
            ne=0
            no=0
            ! compute the mean wilson loop
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                n=(/e1,e2,e3,e4/)
                do mu=1,3
                do nu=mu+1,4
                    if(mod(numberofvortices(z,n,mu,nu,i,j),2.0_r8).gt.0.5_r8)then
                        sumo=sumo+wilsonplanar(u,n,mu,nu,i,j)
                        no=no+1
                    else
                        sume=sume+wilsonplanar(u,n,mu,nu,i,j)
                        ne=ne+1
                    endif
                    ! the lattice is simmetric, so we
                    ! must consider the case a->b and b->a
                    ! now we interchange a and b
                    if(i.ne.j)then 
                        if(mod(numberofvortices(z,n,mu,nu,j,i),2.0_r8).gt.0.5_r8)then
                            sumo=sumo+wilsonplanar(u,n,mu,nu,j,i)
                            no=no+1
                        else
                            sume=sume+wilsonplanar(u,n,mu,nu,j,i)
                            ne=ne+1
                        endif
                    endif
                    ! the loop for this site is calculated!
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
            we(i,j)=sume/ne
            wo(i,j)=sumo/no
        enddo
        enddo
    endsubroutine vortexsperatedwilson

    subroutine vortexcomponets(wf,we,wo,a,b)
        real(kind=r8), dimension(a,b) :: wf, wo, we
        integer(kind=i4) :: a, b, i, j, area
        do i=1,a
        do j=1,b
            area=i*j
            write(500+area,*) area, wf(i,j), we(i,j), wo(i,j)
        enddo
        enddo
    endsubroutine vortexcomponets

    subroutine vortexlimitedwilson(u,z,a,b)
        real(kind=r8), dimension(nr,nr,nr,nt,4,4) :: u, z
        real(kind=r8) :: w0, w1, w2
        integer(kind=i4) :: ia, ib, a, b, n0, n1, n2, nv
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        do ia=1,a
        do ib=1,b
            w0=0.0_r8
            w1=0.0_r8
            w2=0.0_r8
            n0=0
            n1=0
            n2=0
            ! we sweep the lattice
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                do mi=1,3
                do ni=mi+1,4
                    ! we select the data
                    nv=int(numberofvortices(z,(/e1,e2,e3,e4/),mi,ni,ia,ib))
                    if(nv.eq.0)then
                        w0=w0+wilsonplanar(u,(/e1,e2,e3,e4/),mi,ni,ia,ib)
                        n0=n0+1
                    elseif(nv.eq.1)then
                        w1=w1+wilsonplanar(u,(/e1,e2,e3,e4/),mi,ni,ia,ib)
                        n1=n1+1
                    elseif(nv.eq.2)then
                        w2=w2+wilsonplanar(u,(/e1,e2,e3,e4/),mi,ni,ia,ib)
                        n2=n2+1
                    endif
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
            write(600+ia*ib,*) ia*ib, w0, w1, w2
        enddo
        enddo
    endsubroutine vortexlimitedwilson

    ! subroutine that measures all kind of wilson loops
    subroutine wilsonvortex(u,ur,z,w,wr,wv,we,wo,a,b)
        real(kind=r8), dimension(nr,nr,nr,nt,4,4), intent(in) :: u, ur, z
        real(kind=r8), dimension(a,b), intent(out) :: w, wr, wv, we, wo
        real(kind=r8) :: sum, sumr, sumv, sume, sumo, w0, w1, w2
        real(kind=r8) :: probeven, probodd, nv
        integer(kind=i4) :: a, b, i, j, ne, no, n(4)
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        ! we compute the loops w(a,b) and save the result
        w=0.0_r8
        wr=0.0_r8
        wv=0.0_r8
        we=0.0_r8
        wo=0.0_r8
        !$OMP PARALLEL DO SHARED(w,wr,wv) PRIVATE(e1,e2,e3,e4,mi,ni,i,j,n)
        do i=1,a
        do j=1,b
            ! sweep the lattice
            ! compute the mean all the kinds of wilson loops
            do mi=1,3
            do ni=mi+1,4
                do e1=1,nr
                do e2=1,nr
                do e3=1,nr
                do e4=1,nt
                n=(/e1,e2,e3,e4/)
                
                    ! =============================================================== !
                    ! compute the normal wilson loops
                    ! full  wilson
                    w(i,j)=w(i,j)+wilsonplanar(u,n,mi,ni,i,j)
                    ! vortex projected wilson
                    wv(i,j)=wv(i,j)+wilsonplanar(z,n,mi,ni,i,j)
                    ! vortex removed wilson
                    wr(i,j)=wr(i,j)+wilsonplanar(ur,n,mi,ni,i,j)

                    ! the lattice is simmetric, so we
                    ! must consider the case a->b and b->a
                    ! now we interchange a and b
                    if(i.ne.j)then 
                        ! full  wilson
                        w(i,j)=w(i,j)+wilsonplanar(u,n,mi,ni,i,j)
                        ! vortex projected wilson
                        wv(i,j)=wv(i,j)+wilsonplanar(z,n,mi,ni,i,j)
                        ! vortex removed wilson
                        wr(i,j)=wr(i,j)+wilsonplanar(ur,n,mi,ni,i,j)
                    endif
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
                wv(i,j)=wv(i,j)/(6.0_r8*nt*nr**3)    ! quadratic loops
                wr(i,j)=wr(i,j)/(6.0_r8*nt*nr**3)    ! quadratic loops
            else
                w(i,j)=w(i,j)/(12.0_r8*nt*nr**3)    ! quadratic loops
                wv(i,j)=wv(i,j)/(12.0_r8*nt*nr**3)    ! quadratic loops
                wr(i,j)=wr(i,j)/(12.0_r8*nt*nr**3)    ! quadratic loops
            endif
        enddo
        enddo

        ! get the data out
        do i=1,a
            write(100+i,*) i, w(i,:)
            write(200+i,*) i, wv(i,:)
            write(300+i,*) i, wr(i,:)
            !do j=1,b
            !    write(5000+i*j,*) i*j, w(i,j), we(i,j), wo(i,j)
            !enddo
        enddo
    endsubroutine wilsonvortex

    ! center vortex vizualization
    subroutine vortex_plot_ss(z,i)
        real(kind=r8), dimension(nr,nr,nr,nr,4,4) :: z
        real(kind=r8) :: sum
        integer(kind=i4) :: i, e1, e2, e3, e4, mi, ni, n(3), nn(3)
        do e4=1,nr
            ! run this spacial section and find the center vortices
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
                do mi=1,2
                do ni=mi+1,3
                    ! we have a vortex
                    if(vortexinplaq(z,(/e1,e2,e3,e4/),mi,ni).lt.0.0_r8)then
                        n=(/e1,e2,e3/)
                        n(mi)=n(mi)+0.5_r8
                        n(ni)=n(ni)+0.5_r8
                        nn=n
                        if((mi==1).and.(ni==2)) n(3)=n(3)-0.5_r8
                        if((mi==1).and.(ni==3)) n(2)=n(2)-0.5_r8
                        if((mi==2).and.(ni==3)) n(1)=n(1)-0.5_r8
                        if((mi==2).and.(ni==1)) n(3)=n(3)-0.5_r8
                        if((mi==3).and.(ni==2)) n(1)=n(1)-0.5_r8
                        if((mi==3).and.(ni==1)) n(2)=n(2)-0.5_r8
                        write(10000,*) n(:),nn(:)-n(:)
                    endif
                enddo
                enddo
            enddo
            enddo
            enddo
            write(10000,*) ' '
            write(10000,*) ' '
        enddo
        write(10000,*) ' '
        write(10000,*) ' '
    endsubroutine vortex_plot_ss
endmodule su2center