program su2statistics
    ! define the precision
    integer, parameter :: r8=selected_real_kind(8,15)
    integer, parameter :: i4=selected_int_kind(8)

    ! 
    contains
    ! function that coputes the mean of a set of values a
    function meana(a,np)
        real(kind=r8) :: a(np), mia, meana
        integer(kind=i4) :: np, i
        mia=0.0_r8
        do i=1,np
            mia=mia+a(i)
        enddo
        meana=mia/np
    endfunction meana

    ! function that computes the standart deviation of an uncorrelated set of values 
    function sigmaa(a,mia,np)
        real(kind=r8) :: a(np), mia, siga, sigmaa
        integer(kind=i4) :: np, i
        sigmaa=0.0_r8
        do i=1,np
            siga=siga+(mia-a(i))**2.0_r8
        enddo
        sigmaa=dsqrt(siga/(np*(np-1)))
    endfunction sigmaa

    ! routine that applies the bining method to a set of values a
    subroutine bining(nf,bmax,filename)
        ! defining ou arrays
        real(kind=r8), dimension(nf) :: a, sigmab, ab
        real(kind=r8) :: abar, aux
        integer(kind=i4) :: nf, i, ib, jb, nb, b, bmax
        character(len=20) :: filename

        ! read the input values
        open(unit=2,file=filename)
        abar=0.0_r8
        do i=1,nf
            read(2,*) aux, a(i)
            abar=abar+a(i)
        enddo
        abar=abar/nf

        close(2)
        open(unit=3,file='bining.dat')
        ! make the block method
        do b=1,bmax
            ! we are making nb bins, then each bin has b points
            ab=0.0_r8
            nb=nf/b
            ! construct the nb blocks
            do ib=0,nb-1
                do jb=1,b
                    ab(ib+1)=ab(ib+1)+a(jb+ib*b)
                enddo
            enddo
            ab=ab/b

            ! compute the standar deviation of this set of size
            sigmab(nb)=0.0_r8
            do ib=1,nb
            
                sigmab(nb)=sigmab(nb)+(ab(ib)-abar)**2.0_r8
            enddo
            sigmab(nb)=sigmab(nb)/(nb-1)
            write(3,*) b,sigmab(nb)
        enddo
        close(3)
    endsubroutine bining

    ! subroutine that computes the autocorrelation (n=nfile-ttrans)
    subroutine autocorrelation(ttrans,tmax,n,filename)
        real(kind=r8), dimension(n) :: x, chi
        integer(kind=i4), dimension(n) :: y
        real(kind=r8) :: chi0, xbar, mu, sum, mu1, mu2, tau, a, b
        integer(kind=i4) :: nfile, ttrans, i, j, n, t, tmax
        character(len=20) :: filename
        open(unit=2,file=filename)

        do i=1,nfile
            read(2,*) a, b
            if(i.gt.ttrans)then
                y(j)=a
                x(j)=b
                mu=mu+x(j)
                j=j+1
            endif
        enddo
        close(2)
        mu=mu/(j-1)

        ! now we calculate the autocorrelation function
        chi=0.0_r8
        do t=1,n!tmax
        !print*,t
            chi0=0.0_r8
            mu1=0.0_r8
            mu2=0.0_r8
            do i=1,n-t!tmax-t
                chi0=chi0+x(i)*x(i+t)
                mu1=mu1+x(i)
                mu2=mu2+x(i+t)
            enddo
            chi(t)=chi0/(n-t)-mu1*mu2/(n-t)**2.0_r8
        enddo
        chi=chi/chi(1)

        ! get out with the data
        open(unit=3,file='chi.dat')
        tau=0.50_r8
        do t=1,tmax
            write(3,*) t-1, chi(t), log(chi(t))
            tau=tau+chi(t)
        enddo
        open(unit=4,file='tauint.dat')
        write(4,*) tau  
        close(3)
    endsubroutine autocorrelation
endprogram su2statistics