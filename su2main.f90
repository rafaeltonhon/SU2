program su2main
    use su2lattice
    use su2thermalize
    use su2measures
    use su2center
    use su2gaugefix
    use omp_lib
    implicit none

    call system_clock(count_rate=rate)
    call system_clock(itimes1)
!    open(unit=1,file='su2-in.dat')
!    read(1,*) nr, nt            ! lattice size space time
!    read(1,*) beta              ! beta
!    read(1,*) nmc, nterm        ! number of monte carlo sweeps number of thermalization sweeps
!    read(1,*) a, b              ! maximum size of the wilson loops we want
!!    read(1,*) iprint            ! if 1 we print the configurations on the file suN-lattice
!    read(1,*) ncorr             ! number of discarted configurations
!    read(1,*) tol               ! gauge fix tolerance
!    read(1,*) initflag          ! init flag 1-> cold init, 2-> hot init
!    read(1,*) function          ! function a ->  termilize the lattice and find the correlation time
                                ! b -> termilize the lattice and generate uncorrelated configurations
                                ! c -> read the configurations and make measures
    close(1)

    !allocate(u(nr,nr,nr,nt,4,4),r(nr,nr,nr,nt,4,4))
    !allocate(ug(nr,nr,nr,nt,4,4),z(nr,nr,nr,nt,4,4))
    !allocate(w(a,b),wv(a,b),wr(a,b),we(a,b),wo(a,b))
    call init(u,beta,initflag)

    select case(function)
    case('a') ! termilize the lattice and find the correlation time
        ! termilize the lattice
        do i=1,nterm
            call hbstep(u,beta)
        enddo
        do i=1,nmc
            call hbstep(u,beta)
            call measurewilson(u,1,1,w)
            write(101,*) i, w(1,1)
        enddo
    case('b') ! termilize the lattice and generate uncorrelated configurations
        do i=1,nterm
            call hbstep(u,beta)
        enddo
        do i=1,nmc
            ! save the actual configuration
            print*,'Configuration',i
            if(iprint.eq.1)call save_lattice(u,nr,i)
            
            ! compute the center vortex observables
            call maximal_center_gauge(u,ug,tol)
            call centerprojection(ug,z)
            call centerremotion(ug,r)
            call wilsonvortex(ug,r,z,w,wr,wv,a,b)
            call vortexsperatedwilson(u,z,a,b,we,wo)
            call vortexcomponets(w,we,wo,a,b)
            call pvortexdensity(z,a,b)
           !call vortexlimitedwilson(u,z,a,b)
           !call vortex_plot_ss(z,i)
            
            ! compute the propagators in the landau gauge
            call landau_gauge(u,ug,tol/1e11) ! full ensemble
            call gluonpropagator(ug,1000) ! complete
            call gluon_formfactor(ug,1100)

            !call centerlandau_gauge(z,z,tol/1e11)
            !call centergluonpropagator(z,2000)
            !call centergluon_formfactor(ug,2100)

            call landau_gauge(r,r,tol/1e11) ! removed ensenmble
            call gluonpropagator(r,3000)
            call gluon_formfactor(r,3100)
        
          	!call measurewilson(u,a,b,w)
	        !do j=1,a
	            !write(100+j,*)	i, w(j,:)
            !enddo
            ! get the next uncorrelated configurations
            do j=1,ncorr ! number of correlated configurations
                ! discart ncorr configurations
                call hbstep(u,beta)
            enddo
        enddo
    case('c')   ! read the configurations and make measures
        do i=1,nmc-1
            ! save the actual configuration
            call read_lattice(u,nr,i)
            call maximal_center_gauge(u,ug,tol)
            !call maximal_oldabelian_gauge(ulink,ugauge,tol)
            !call centerprojection(ugauge,z)
            !call centerremotion(ugauge,ur)
            !call wilsonvortex(ugauge,ur,z,w,wr,wv,we,wo,a,b)
        enddo
    case('d') 
        do i=1,nterm
            call hbstep(u,beta)
        enddo
            do j=1,100 ! number of correlated configurations
                call cooling(u,15)
                call hbstep(u,beta)
            enddo
            do i=1,nmc
                ! save the actual configuration
                print*,'Configuration',i
                if(iprint.eq.1)call save_lattice(u,nr,i)
                
                call cooling(u,10)
                call measurewilson_asym(u,a,b,w)
                do j=1,a
                    write(100+j,*) w(j,:)
                enddo
                ! get the next uncorrelated configurations
                do j=1,ncorr ! number of correlated configurations
                    ! discart ncorr configurations
                    call hbstep(u,beta)
                enddo
            enddo
    case default
        stop 'no valid case was selected'
    endselect

    !deallocate(u,ug,z,r,wr,w,wv)
    call system_clock(itimes2)
    print*,"cpu time [seconds]: ", (itimes2-itimes1)/float(rate)!(end-start)/60
    write(9000,*) "cpu time [seconds]: ", (itimes2-itimes1)/float(rate)
    write(9000,*) "cpu time [minutes]: ", (itimes2-itimes1)/float(60*rate)
    write(9000,*) "cpu time [hours]: ", (itimes2-itimes1)/float(3600*rate)
end program su2main
