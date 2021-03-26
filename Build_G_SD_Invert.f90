      subroutine Build_G_SD_Invert(G_C,Sigma_l,Sigma_r,energy,size_l,size_c,size_lc,size_lcr,norb,inputfile,numocc,numvirt,iter,B1data,B2data,mo_ener,mo_coeff,mo_coeff2,doubles,currentflag,energy_values,ener_val,G_S,corr_ener,numatomic,B0_Coeff,use_b0,gamess,maple,numfcore,numfvirt,b0_type)
      use InterfaceMod2
      use TypeMod
      use FunctionMod
      implicit none
      character(len=40) :: inputfile,Bfile,Bfile2,readtemp,cstring,b0_type
      integer :: size_l,size_lc,size_c,size_lcr,i,j,ioerror,norb,iter,a,b,counter
      integer :: numocc,numvirt,k,p,q,z,y,x,r,s,t,energy_values,ener_val,numatomic
      integer :: aa,bb,pp,qq,ii,jj,tt,rr,ss,xx,yy,zz
      integer :: tempx, tempy,numvirt_act,numocc_act,numfvirt,numfcore
      real(8) :: temp,corr_ener,hf_energy,energy,tempval,dft_energy,use_dft_mo
      complex(8), allocatable, dimension(:,:) :: G_C,G_temp,Sigma_L,Sigma_R
      type(B2) :: B2data
      type(B1) :: B1data
      real(8) :: cisd_b0
      real(8), allocatable, dimension(:,:) :: mo_coeff,mo_coeff2,temp_gf,temp_gf2
      real(8), allocatable, dimension(:) :: mo_ener,B0_coeff
      logical :: singles,doubles,currentflag,use_b0,maple,gamess
      type(energr) :: G_S
      integer(8) :: rt1,rt2,rt3,rt4,rt5,io 
      numocc_act=numocc-numfcore
      numvirt_act=numvirt-numfvirt
      if(doubles.eqv..true.) then
        singles = .false.
       else if(doubles.eqv..false.) then
        singles = .true.
      end if

      10 format(A)
      if(iter.eq.1) then
       if(doubles.eqv..true.) then
       allocate(B1data%a%o(numocc+1:numocc+numvirt_act,numfcore+1:numocc))
       allocate(B1data%b%o(numocc+1:numocc+numvirt_act,numfcore+1:numocc))
       allocate(B2data%aa%m(numocc+1:numocc+numvirt_act,numfcore+1:numocc))
       allocate(B2data%ab%m(numocc+1:numocc+numvirt_act,numfcore+1:numocc))
       B1data%a%o=0
       B1data%b%o=0
       end if
       allocate(mo_ener(1:norb))
       allocate(mo_coeff(1:numatomic,1:norb))
       allocate(mo_coeff2(1:norb,1:numatomic))
       mo_coeff = 0
       mo_ener = 0
       mo_coeff2 = 0
       if(doubles.eqv..true.) then
       counter=1
       do j=numfcore+1,numocc
        do a = numocc+1,numocc+numvirt_act
          allocate(B2data%aa%m(a,j)%n(numocc+1:numocc+numvirt_act,numfcore+1:numocc))
          allocate(B2data%ab%m(a,j)%n(numocc+1:numocc+numvirt_act,numfcore+1:numocc))
            B2data%aa%m(a,j)%n=0
            B2data%ab%m(a,j)%n=0
         end do
        end do
       if(gamess.eqv..true.) then
         Bfile = trim(inputfile) // "T2"
         open(unit=4,file=Bfile,action='READ',iostat = ioerror)
         read(4,*) readtemp
         write(*,*) "GAMESS version does not work with frozen orbitals."
         write(*,*) "Do not use any frozen or core orbitals. ECP only."
         write(*,*) 'Total Number of T2(ab) and T2(aa) Elements',readtemp
        do i=1,numocc_act
         do a=numocc_act+1,numocc_act+numvirt_act
          do j=1,numocc_act
            do b=numocc_act+1,numvirt_act+numocc_act
              read(4,*) B2data%aa%m(a,i)%n(b,j)
              counter=counter+1
             end do
            end do
           end do
         end do
        do i=1,numocc_act
         do a=numocc_act+1,numocc_act+numvirt_act
          do j=1,numocc_act
            do b=numocc_act+1,numvirt_act+numocc_act
              read(4,*) B2data%ab%m(a,i)%n(b,j)
              counter=counter+1
             end do
            end do
           end do
         end do

      ! do i=1,numocc
      !   do a=numocc+1,numocc+numvirt
      !    do j=1,numocc
      !     do b=numocc+1,numvirt+numocc
             !write(*,*) i,a,j,b
       !      read(4,*) readtemp!B2data%bb%m(a,i)%n(b,j)
       !     end do
       !    end do
       !   end do
       ! end do


       do i=1,numocc_act
        do a=numocc_act+1,numvirt_act+numocc_act
          read(4,*) B1data%a%o(a,i)
             counter=counter+1
         end do
       end do
       write(*,*) 'Number of T2(ab) and T1(a) Values in T2 file',counter-1

        elseif(maple.eqv..true.) then
         Bfile = trim(inputfile) // ".T2"
         open(unit=4,file=Bfile,action='READ',iostat = ioerror)
         read(4,*) cstring,readtemp,corr_ener,use_dft_mo,hf_energy,dft_energy
         read(4,*) use_dft_mo,hf_energy,dft_energy

          write(*,*) 'Max Number of T2(ab) and T2(aa) Elements',readtemp
          if(use_dft_mo.eq.1) then
            write(*,*) "Using DFT MO Energies"
            corr_ener=hf_energy+corr_ener-dft_energy
           end if
          write(*,*) 'Difference in Energy',corr_ener

          counter=0

         do
           counter=counter+1
           read(4,*,IOSTAT=io) rt5,rt1,rt2,rt3,rt4,tempval
           if(io.eq.0..and.rt5.ne.-1) then
            if(rt5.eq.2.and.rt1.ne.0) then
              B2data%aa%m(rt4+numfcore,rt1+numfcore)%n(rt3+numfcore,rt2+numfcore)=tempval
            elseif(rt5.eq.1.and.rt1.ne.0) then
              B2data%ab%m(rt4+numfcore,rt1+numfcore)%n(rt3+numfcore,rt2+numfcore)=tempval
            elseif(rt5.eq.0.and.rt1.ne.0) then
              !add T2(bb) here once open shell is complete
            elseif(rt5.eq.1.and.rt1.eq.0) then
              B1data%a%o(rt4+numfcore,rt3+numfcore)=tempval
            elseif(rt5.eq.0.and.rt1.eq.0) then
              B1data%b%o(rt4+numfcore,rt3+numfcore)=tempval
            end if
           else
            write(*,*) "Done reading T2 file. ",counter," terms read."
            exit
          end if
         end do



       end if

       close(4)

       end if
       if(doubles.eqv..true.) then
        if(use_b0.eqv..false.) then
          write(*,*) "Using CEPA conventation for B0"
         elseif(use_b0.eqv..true.) then
         if(trim(b0_type).eq."cisd") then
          write(*,*) "Using CISD conventation for B0"
          Call Build_B0_CISD(cisd_b0,norb,numfcore,numvirt_act,numocc,numvirt,B1data,B2data)
         end if
        end if
       end if

       if(gamess.eqv..true.) then
       write(*,*) 'Reading MO info'
       Bfile2 = trim(inputfile) // ".mo_dat"
       open(unit=5,file=Bfile2,action='READ',iostat = ioerror)
        do i=1,numatomic
         do j=1,norb
           read(5,*) mo_coeff(i,j)
           mo_coeff2(j,i) = mo_coeff(i,j)
          end do
         end do

        do i=1,norb
          read(5,*) mo_ener(i)
         end do      
        read(5,*) corr_ener
        write(*,*) 'Correlation Energy',corr_ener

        elseif(maple.eqv..true.) then
         write(*,*) "Reading MO data from .scf_dat file"
         Bfile2 = trim(inputfile) // ".scf_dat"
         open(unit=5,file=Bfile2,action='READ',iostat = ioerror)
         do while(trim(readtemp)/="Molecular Orbital Coefficients")
            read(5,'(A)') readtemp
           end do
         do i=1,numatomic
          do j=1,norb
            read(5,*) tempx,tempy,tempval
            mo_coeff(tempx,tempy)=tempval
            mo_coeff2(j,i) = mo_coeff(i,j)
           end do
          end do
        read(5,'(A)')
        do i=1,norb
          read(5,*) tempx,mo_ener(i)
         end do

       end if

      allocate(G_S%en(1:energy_values))
      do i = 1,energy_values
      allocate(G_S%en(i)%gf(1:numatomic,1:numatomic))
       G_S%en(i)%gf = 0.0
      end do
      end if

      if(currentflag.eqv..false.) then
        allocate(temp_gf(1:norb,1:norb))
       if(singles.eqv..true.) then
        corr_ener = 0
       end if
       temp_gf = 0.0
      if(doubles.eqv..true.) then
       do i=numfcore+1,numocc
          if(use_b0.eqv..true..and.trim(b0_type).eq."cisd") then
           temp_gf(i,i) = cisd_b0*(energy-(corr_ener+mo_ener(i))*27.211396132)**(-1)
            else if(use_b0.eqv..true..and.trim(b0_type).eq."rdm") then
            temp_gf(i,i) = b0_coeff(i)**(2)*((energy-(corr_ener+mo_ener(i))*27.211396132)**(-1))
            else if(use_b0.eqv..false.) then
           temp_gf(i,i) =(energy-(corr_ener+mo_ener(i))*27.211396132)**(-1)
          end if
        end do 
       do a=numocc+1,numvirt_act+numocc
          if(use_b0.eqv..true..and.trim(b0_type).eq."cisd") then
           temp_gf(a,a) = cisd_b0*(energy+(corr_ener-mo_ener(a))*27.211396132)**(-1)
           else if(use_b0.eqv..true..and.trim(b0_type).eq."rdm") then
            temp_gf(a,a) = b0_coeff(a)**(2)*((energy+(corr_ener-mo_ener(a))*27.211396132)**(-1))
           else if(use_b0.eqv..false.) then
           temp_gf(a,a) = (energy+(corr_ener-mo_ener(a))*27.211396132)**(-1)
          end if
       end do 

       if(numfcore.gt.0) then
        do i=1,numfcore
           temp_gf(i,i) = (energy-(corr_ener+mo_ener(i))*27.211396132)**(-1)
         end do
       end if

       if(numfvirt.gt.0) then
       do a=numocc+1+numvirt_act,numocc+numvirt
          temp_gf(a,a) = (energy+(corr_ener-mo_ener(a))*27.211396132)**(-1)
        end do
       end if

       else

        do i=1,numocc
           temp_gf(i,i) = (energy-(mo_ener(i))*27.211396132)**(-1)
         end do

        do a=numocc+1,numvirt+numocc
          temp_gf(a,a) = (energy-(mo_ener(a))*27.211396132)**(-1)
        end do

      end if

      if(singles.eqv..true.) then
        goto 100
       else if(doubles.eqv..true.) then

       do a=numocc+1,numocc+numvirt_act
        do b=numocc+1,numocc+numvirt_act
         do q=numfcore+1,numocc
 
           temp_gf(a,b) = temp_gf(a,b) + B1data%a%o(a,q)*B1data%a%o(b,q)*(energy-(corr_ener+mo_ener(q))*27.211396132)**(-1)

          end do
         end do
        end do


        do a=numfcore+1,numocc
        do b=numfcore+1,numocc
         do p=numocc+1,numocc+numvirt_act

           temp_gf(a,b) = temp_gf(a,b) + B1data%a%o(p,a)*B1data%a%o(p,b)*(energy+(corr_ener-mo_ener(p))*27.211396132)**(-1)

          end do
         end do
        end do

       do a=numfcore+1,numocc
        do y=numfcore+1,numocc
         do z=numocc+1,numocc+numvirt_act

           temp_gf(a,a) = temp_gf(a,a) + B1data%a%o(z,y)*B1data%a%o(z,y)*(energy+(-corr_ener-mo_ener(a)-mo_ener(y)+mo_ener(z))*27.211396132)**(-1)

          end do
         end do
        end do

       do a=numocc+1,numocc+numvirt_act
        do b=numocc+1,numocc+numvirt_act
         do z=numocc+1,numocc+numvirt_act
          do x=numfcore+1,numocc
           do y=numfcore+1,numocc

           temp_gf(a,b) = temp_gf(a,b) + (B2data%ab%m(a,x)%n(z,y)*B2data%ab%m(b,x)%n(z,y))*((energy+(-corr_ener-mo_ener(x)-mo_ener(y)+mo_ener(z))*27.211396132)**(-1))
           if(x.gt.y.and.a.gt.z.and.b.gt.z) then
             temp_gf(a,b) = temp_gf(a,b) + (B2data%aa%m(a,x)%n(z,y)*B2data%aa%m(b,x)%n(z,y))*(energy+(-corr_ener-mo_ener(x)-mo_ener(y)+mo_ener(z))*27.211396132)**(-1)
            end if
            
          end do
         end do
        end do
       end do
      end do

      do a=numocc+1,numocc+numvirt_act
       do t=numfcore+1,numocc
         do r=numocc+1,numocc+numvirt_act

           temp_gf(a,a) = temp_gf(a,a) + B1data%a%o(r,t)*B1data%a%o(r,t)*(energy+(corr_ener-mo_ener(r)-mo_ener(a)+mo_ener(t))*27.211396132)**(-1)

          end do
         end do
        end do

       do a=numfcore+1,numocc
        do b=numfcore+1,numocc
         do t=numfcore+1,numocc
          do r=numocc+1,numocc+numvirt_act
           do s=numocc+1,numocc+numvirt_act

           if(r.gt.s.and.a.gt.t.and.b.gt.t) then
            temp_gf(a,b) = temp_gf(a,b) + (B2data%aa%m(r,t)%n(s,b)*B2data%aa%m(r,t)%n(s,a))*(energy+(corr_ener-mo_ener(r)-mo_ener(s)+mo_ener(t))*27.211396132)**(-1)
            end if

          end do
         end do
        end do
       end do
      end do

        end if
100     allocate(temp_gf2(1:norb,1:numatomic))

        temp_gf2 = 0
        call matmul_dgemm2(temp_gf,mo_coeff2,temp_gf2)
        call matmul_dgemm2(mo_coeff,temp_gf2,G_S%en(ener_val)%gf)
        G_S%en(ener_val)%gf = inv_real(G_S%en(ener_val)%gf)
        deallocate(temp_gf)
        deallocate(temp_gf2)

       end if

       G_C = CMPLX(0,0)
       G_C=G_S%en(ener_val)%gf(size_l+1:size_lc,size_l+1:size_lc)
       G_C = G_C - Sigma_l - Sigma_r
       G_C = inv(G_C)

      end subroutine
