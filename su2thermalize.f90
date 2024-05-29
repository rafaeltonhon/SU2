module su2thermalize
    use su2lattice
    use omp_lib
    implicit none
    contains
    function hbmatrix(a,beta)
        real(kind=r8) :: hbmatrix(4), x(4), r(3)
        real(kind=r8) :: lambda2, rr, a, beta
        integer(kind=i4) :: i, j

        do j=1,200
            do i=1,3
                call random_number(r(i))
                r(i)=1.0_r8-r(i)
            enddo
            lambda2=-0.5_r8*(dlog(r(1))+(dcos(2.0_r8*pi*r(2)))**2.0_r8*dlog(r(3)))/(a*beta)
            call random_number(rr)
            rr=1.0_r8-rr
            if(rr*rr.lt.(1.0_r8-lambda2))exit
        enddo
        x(4)=1.0_r8-2.0_r8*lambda2

        do j=1,2000
            do i=1,3
                call random_number(r(i))
                r(i)=1.0_r8-2.0_r8*r(i)
            enddo
            if((r(1)*r(1)+r(2)*r(2)+r(3)*r(3)).lt.1.0_r8)exit
        enddo

        do i=1,3
            x(i)=dsqrt((1.0_r8-x(4)*x(4)))*r(i)/dsqrt((r(1)*r(1)+r(2)*r(2)+r(3)*r(3)))
        enddo
        hbmatrix=x
    endfunction hbmatrix

    ! function that makes the heatbath steap
    subroutine hbstep(u,beta)
        real(kind=r8), intent(inout) :: u(nr,nr,nr,nt,4,4)
        real(kind=r8) ::  sums(4), a, beta
        integer(kind=i4) :: e1, e2, e3, e4, mi, i
        !$OMP PARALLEL DO PRIVATE(e1,e2,e3,e4,mi,a,sums) SHARED(u)
        do mi=1,4
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                sums=sumstamples(u,e1,e2,e3,e4,mi)
                a=dsqrt(detlink(sums))
                if(a.gt.5e-3)then   ! det sum_stamples > 0
                    u(e1,e2,e3,e4,mi,:)=linkmult(hbmatrix(a,beta),dlink(sums))/a
                else
                    do i=1,4
                        call random_number(u(e1,e2,e3,e4,mi,i))
                    enddo
                endif
                u(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)/dsqrt(detlink(u(e1,e2,e3,e4,mi,:)))
            enddo
            enddo
            enddo
            enddo
        enddo
        !$END OMP PARALLEL DO
    endsubroutine hbstep

    subroutine init(u,beta,type)
        real(kind=r8), intent(inout) :: u(nr,nr,nr,nt,4,4)
        real(kind=r8) :: beta
        integer(kind=i4) :: e1, e2, e3, e4, mi, type
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nt
            do mi=1,4
                if(type==1)then         ! cold init
                    u(e1,e2,e3,e4,mi,1:3)=0.0_r8
                    u(e1,e2,e3,e4,mi,4)=1.0_r8
                elseif(type==2)then         ! hot init
                    u(e1,e2,e3,e4,mi,:)=hbmatrix(1.0_r8,beta)
                !elseif(type==3)then    ! continue from a previous configuration
                endif
            enddo
        enddo
        enddo
        enddo
        enddo
    endsubroutine init

    ! subroutine that makes one step of smearing
    subroutine smearing(u,usmeared)
        real(kind=r8), dimension(nr,nr,nr,nr,4,4) :: u, usmeared
        real(kind=r8) :: sums(4), temp(4)
        integer(kind=i4) :: x, y, z, t, mi, ni
        integer(kind=i4), dimension(4) :: nzz, npz, nmz, nzp, nzm, npp, npm, nmp
        usmeared=0.0_r8
        ! run over the lattice
        do x=1,nr
        do y=1,nr
        do z=1,nr
        do t=1,nr
            do mi=1,4
                ! compute the sum of the stamples
                nzz=(/x,y,z,t/)
                sums=0.0_r8
                do ni=1,4
                    if (ni.ne.mi) then
                        call neighborhood(nzz,npz,nmz,nzp,nzm,npp,npm,nmp,mi,ni)
                        temp=ident()
                        temp=linkmult(temp,u(nzz(1),nzz(2),nzz(3),nzz(4),ni,:))
                        temp=linkmult(temp,u(nzp(1),nzp(2),nzp(3),nzp(4),mi,:))
                        temp=linkmult(temp,dlink(u(npz(1),npz(2),npz(3),npz(4),ni,:)))
                        
                        sums=sums+temp
                        temp=ident()
                        temp=linkmult(temp,dlink(u(nzm(1),nzm(2),nzm(3),nzm(4),ni,:)))
                        temp=linkmult(temp,u(nzm(1),nzm(2),nzm(3),nzm(4),mi,:))
                        temp=linkmult(temp,u(npm(1),npm(2),npm(3),npm(4),ni,:))
                        sums=sums+temp
                    endif
                enddo
                !usmeared(x,y,z,t,mi,:)=ulink(x,y,z,t,mi,:)+sums
                usmeared(x,y,z,t,mi,:)=usmeared(x,y,z,t,mi,:)/dsqrt(detlink(usmeared(x,y,z,t,mi,:)))
            enddo
        enddo
        enddo
        enddo
        enddo
    endsubroutine smearing
endmodule su2thermalize