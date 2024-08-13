module su2lattice
    use omp_lib
    implicit none
    ! define the precision
    integer, parameter :: r8=selected_real_kind(8,15)
    integer, parameter :: i4=selected_int_kind(8)

    ! define the parameters
    real(kind=r8), parameter :: pi=acos(-1.0_r8)
    real(kind=r8), dimension(4), parameter :: sigma3=(/0.0_r8,0.0_r8,1.0_r8,0.0_r8/)
    integer(kind=i4), parameter :: nr=20 ! number of sites in space
    integer(kind=i4), parameter :: nt=20 !  number of sites in time
    integer(kind=i4), parameter :: a=10 ! maximum size of the wilson loops we want
    integer(kind=i4), parameter :: b=10 ! maximum size of the wilson loops we want
    real(kind=r8), parameter :: beta=2.30_r8              ! beta
    integer(kind=i4), parameter :: nmc=1000  ! number of monte carlo sweeps number 
    integer(kind=i4), parameter :: nterm=1000        ! number of thermalization sweeps
    integer(kind=i4), parameter :: iprint=2            ! if 1 we print the configurations on the file suN-lattice
    integer(kind=i4), parameter :: ncorr=200             ! number of discarted configurations
    real(kind=r8), parameter :: tol=2e-5               ! gauge fix tolerance
    integer(kind=i4), parameter :: initflag=2          ! init flag 1-> cold init, 2-> hot init
    character(len=1), parameter :: function='b'          ! function a ->  termilize the lattice and find the correlation time
                                ! b -> termilize the lattice and generate uncorrelated configurations
                                ! c -> read the configurations and make measures

    ! define the matrices
    !real(kind=r8), allocatable, dimension(:,:,:,:,:,:) :: u, ug, r, z
    !real(kind=r8), allocatable, dimension(:,:) :: w, wr, wv, we, wo
    real(kind=r8), dimension(nr,nr,nr,nr,4,4) :: u, ug, r, z
    real(kind=r8), dimension(a,b) :: w, wr, wv, we, wo

    ! defining our variables
    !real(kind=r8) :: beta, nr, nt, a, b, initflag, tol
    integer(kind=i4) :: rate, itimes1, itimes2
    integer(kind=i4) :: i, j
    ! nmc, nterm, ncorr, iprint
    !character(len=1) :: function
    !integer(kind=i4) :: a, b
    contains
    function ident()
        real(kind=r8) :: ident(4)
        ident=(/0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8/)
    endfunction ident

    function crossprod(a,b)
        real(kind=r8), dimension(4) :: a, b
        real(kind=r8), dimension(3) ::  c, crossprod
        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)
        crossprod=c
    endfunction crossprod

    function dot3vec(a,b)
        real(kind=r8), dimension(4) :: a, b
        real(kind=r8):: dot3vec
        dot3vec=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
    endfunction dot3vec

    function linkmult(a,b)
        real(kind=r8), intent(in), dimension(4) :: a, b
        real(kind=r8), dimension(4) :: linkmult
        real(kind=r8) :: ab(3), c(4)
        ab=crossprod(a,b)
        c(4)=a(4)*b(4)-dot3vec(a,b)

        c(1)=b(4)*a(1)+a(4)*b(1)-ab(1)
        c(2)=b(4)*a(2)+a(4)*b(2)-ab(2)
        c(3)=b(4)*a(3)+a(4)*b(3)-ab(3)
        linkmult=c
    endfunction linkmult

    function trlink(u)
        real(kind=r8) :: trlink, u(4)
        trlink=2.0_r8*u(4)
    endfunction trlink

    ! function that computes the determinant of su(2) matrices
    function detlink(u)
        real(kind=r8), intent(in) :: u(4)
        real(kind=r8) :: detlink
        detlink=u(1)*u(1)+u(2)*u(2)+u(3)*u(3)+u(4)*u(4)
    endfunction detlink

    ! function that returns the dagger of the link
    function dlink(u)
        real(kind=r8), dimension(4) :: dlink
        real(kind=r8), intent(in), dimension(4) :: u!, v
        !v(1)=-u(1)
        !v(2)=-u(2)
        !v(3)=-u(3)
        !v(4)=u(4)
        dlink=(/-u(1),-u(2),-u(3),u(4)/)
    endfunction dlink

    ! subroutine that returns the neighbors of a given link
    ! ONLY FOR SQUARE LATTICES
    subroutine neighborhood(nzz,npz,nmz,nzp,nzm,npp,npm,nmp,mi,ni)
        integer(kind=i4), dimension(4) :: nzz,npz,nmz,nzp,nzm,npp,npm,nmp
        integer(kind=i4) :: mi, ni
        ! start the neighboors vectors
        npz=nzz
        nmz=nzz
        nzp=nzz
        nzm=nzz
        npp=nzz
        npm=nzz
        nmp=nzz
        ! in mi direction
        ! plus mi
        npz(mi)=npz(mi)+1
        npp(mi)=npp(mi)+1
        npm(mi)=npm(mi)+1
        if(npz(mi).gt.nr)then
            npz(mi)=1
            npp(mi)=1
            npm(mi)=1
        endif
        ! minus mi
        nmz(mi)=nmz(mi)-1
        nmp(mi)=nmp(mi)-1
        if(nmz(mi).lt.1)then
            nmz(mi)=nr
            nmp(mi)=nr
        endif

        ! in ni direction
        ! plus ni
        nzp(ni)=nzp(ni)+1
        npp(ni)=npp(ni)+1
        nmp(ni)=nmp(ni)+1
        if(nzp(ni).gt.nr)then
            nzp(ni)=1
            npp(ni)=1
            nmp(ni)=1
        endif

        ! minus ni
        nzm(ni)=nzm(ni)-1
        npm(ni)=npm(ni)-1
        if(nzm(ni).lt.1)then
            nzm(ni)=nr
            npm(ni)=nr
        endif
    endsubroutine neighborhood

    ! function that computes the plaquette
    function plaquette(u)
        real(kind=r8) :: u(nr,nr,nr,nt,4,4), sum, w(4), plaquette
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        integer(kind=i4), dimension(4) :: nzz,npz,nmz,nzp,nzm,npp,npm,nmp
        sum=0.0_r8
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nt
            nzz=(/e1,e2,e3,e4/)
            do mi=1,3
            do ni=mi+1,4
                call neighborhood(nzz,npz,nmz,nzp,nzm,npp,npm,nmp,mi,ni)
                w=ident()
                !w=linkmult(w,dlink(u(nzz(1),nzz(2),nzz(3),nzz(4),ni,:)))
                !w=linkmult(w,dlink(u(nzp(1),nzp(2),nzp(3),nzp(4),mi,:)))
                !w=linkmult(w,u(npz(1),npz(2),npz(3),npz(4),ni,:))
                !w=linkmult(w,u(nzz(1),nzz(2),nzz(3),nzz(4),mi,:))
                w=linkmult(w,u(nzz(1),nzz(2),nzz(3),nzz(4),mi,:))
                w=linkmult(w,u(npz(1),npz(2),npz(3),npz(4),ni,:))
                w=linkmult(w,dlink(u(nzp(1),nzp(2),nzp(3),nzp(4),mi,:)))
                w=linkmult(w,dlink(u(nzz(1),nzz(2),nzz(3),nzz(4),ni,:)))
                sum=sum+w(4)
            enddo
            enddo
        enddo
        enddo
        enddo
        enddo
        plaquette=sum/(6.0_r8*nr*nr*nr*nt)
    endfunction plaquette

    ! funtion that compute the sum of the stamples
    function sumstamples(u,e1,e2,e3,e4,mi)
        real(kind=r8) :: u(nr,nr,nr,nt,4,4)
        real(kind=r8), dimension(4) :: sumstamples, temp, sum
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        integer(kind=i4), dimension(4) :: nzz, npz, nmz, nzp, nzm, npp, npm, nmp

        nzz=(/e1,e2,e3,e4/)
        sum=0.0_r8
        do ni=1,4
            if (ni.ne.mi) then
                call neighborhood(nzz,npz,nmz,nzp,nzm,npp,npm,nmp,mi,ni)
                temp=ident()
                temp=linkmult(temp,u(npz(1),npz(2),npz(3),npz(4),ni,:))
                temp=linkmult(temp,dlink(u(nzp(1),nzp(2),nzp(3),nzp(4),mi,:)))
                temp=linkmult(temp,dlink(u(nzz(1),nzz(2),nzz(3),nzz(4),ni,:)))
                
                sum=sum+temp
                temp=ident()
                temp=linkmult(temp,dlink(u(npm(1),npm(2),npm(3),npm(4),ni,:)))
                temp=linkmult(temp,dlink(u(nzm(1),nzm(2),nzm(3),nzm(4),mi,:)))
                temp=linkmult(temp,u(nzm(1),nzm(2),nzm(3),nzm(4),ni,:))
                sum=sum+temp
            endif
        enddo
        sumstamples=sum
    endfunction sumstamples

    subroutine save_lattice(u,nr,iconf)
        real(kind=r8) :: u(nr,nr,nr,nr,4,4)
        integer(kind=i4) :: nr, iconf, e1, e2, e3, e4, mi, seed(1)
        character(len=60) :: conf, str, b

        ! create the file of the form
        ! conf_nre4_sweep=i.dat
        write(b,'(f4.2)') beta
        conf='configurations/beta='
        conf=trim(conf)//trim(b)
        conf=trim(conf)//'_lattice='
        if(nr.lt.10) write(str,'(i1)') nr
        if(nr.ge.10) write(str,'(i2)') nr

        conf=trim(conf)//trim(str)
        conf=trim(conf)//'e4_sweep='

        if(iconf.lt.10) write(str,'(i1)') iconf
        if((iconf.ge.10).and.(iconf.lt.100)) write(str,'(i2)') iconf
        if((iconf.ge.100).and.(iconf.lt.1000)) write(str,'(i3)') iconf
        if((iconf.ge.1000).and.(iconf.lt.10000)) write(str,'(i4)') iconf

        conf=trim(conf)//trim(str)
        conf=trim(conf)//'.dat'
        print*,conf
        open(unit=10,file=conf)

        ! write the configuration
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            do mi=1,4
                write(10,*) e1, e2, e3, e4, mi, u(e1,e2,e3,e4,mi,:)
            enddo
        enddo
        enddo
        enddo
        enddo
        close(10)
    endsubroutine save_lattice

    subroutine read_lattice(u,nr,iconf)
        real(kind=r8), intent(inout) :: u(nr,nr,nr,nr,4,4)
        real(kind=r8) :: uaux(4)
        integer(kind=i4) :: nr, iconf, e1, e2, e3, e4, mi, i, seed(1)
        character(len=60) :: conf, str, b

        ! create the file of the form
        ! conf_nre4_sweep=i.dat
        write(b,'(f4.2)') beta
        conf='configurations/beta='
        conf=trim(conf)//trim(b)
        conf=trim(conf)//'_lattice='
        if(nr.lt.10) write(str,'(i1)') nr
        if(nr.ge.10) write(str,'(i2)') nr

        conf=trim(conf)//trim(str)
        conf=trim(conf)//'e4_sweep='

        if(iconf.lt.10) write(str,'(i1)') iconf
        if((iconf.ge.10).and.(iconf.lt.100)) write(str,'(i2)') iconf
        if((iconf.ge.100).and.(iconf.lt.1000)) write(str,'(i3)') iconf
        if((iconf.ge.1000).and.(iconf.lt.10000)) write(str,'(i4)') iconf

        conf=trim(conf)//trim(str)
        conf=trim(conf)//'.dat'
        print*,conf
        open(unit=10,file=conf)

        ! write the configuration
        do i=1,nr**3*nt*4
            read(10,*) e1,e2,e3,e4,mi,uaux(:)
            u(e1,e2,e3,e4,mi,:)=uaux
        enddo
        close(10)
    endsubroutine read_lattice
endmodule su2lattice