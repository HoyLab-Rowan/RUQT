      subroutine Build_B0_CISD(cisd_b0,norb,numfcore,numfvirt,numocc,numvirt,B1data,B2data)
      !use InterfaceMod
      use TypeMod
      use FunctionMod
      implicit none
      integer :: size_l,size_lc,size_c,size_lcr,i,j,ioerror,norb,iter,a,b,numfvirt
      integer :: numocc,numvirt,p,k,q,z,y,x,r,s,t,energy_values,numatomic,numfcore
      type(B2) :: B2data
      type(B1) :: B1data
      real(8) :: cisd_b0
      real(8) :: sum_B0,sum_B0_2,sum_B0_1

      write(*,*) 'Building B0'

      sum_B0_1=0

         do b=numocc+1,numocc+numfvirt
           do k=numfcore+1,numocc
 
           sum_B0_1 = sum_B0_1 - B1data%a%o(b,k)*B1data%a%o(b,k)

          end do
         end do


        sum_B0_2=0
        do a=numocc+1,numocc+numfvirt
         do k=numfcore+1,numocc
           do b=a,numfvirt+numocc
            do j=k,numocc
             sum_B0_2 = sum_B0_2 - B2data%ab%m(a,k)%n(b,j)*B2data%ab%m(a,k)%n(b,j)
              !write(*,*) a,k,b,j 
             if(k.ne.j.and.a.ne.b) then
              sum_B0_2 = sum_B0_2 - B2data%aa%m(a,k)%n(b,j)*B2data%aa%m(a,k)%n(b,j)
             end if

             end do
            end do
           end do
          end do

          cisd_b0 =1+sum_B0_1+sum_B0_2

          write(*,*) "B0 Final",cisd_b0
       end subroutine
