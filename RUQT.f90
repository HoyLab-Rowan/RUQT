      program RUQT
      Use InterfaceMod
      Use FunctionMod
      !Use InterfaceMod2
      Use TypeMod
      implicit none

      real(8),allocatable,dimension(:,:) :: H_one,H_Two,OneInts,H_Two_le,H_Two_re,H_Two_cen,H_Two_le_trans,H_Two_re_trans
      real(8), allocatable, dimension(:) :: TwoIntsCompact,transm,transm_curr,current,energy_list,mo_ener,B0_coeff
      real(8), allocatable, dimension(:,:) :: Smat_le,Smat_re,Smat_cen,Smat,coupling_mat,mo_coeff,mo_coeff2
      complex(8), allocatable, dimension(:,:) :: gfc_r,gfc_a,current_temp,Sigma_l,Sigma_r,Gamma_L,Gamma_R
      complex(8),allocatable,dimension(:) :: voltage
      real(8) :: coupling_r,coupling_l,energy,energy_start,energy_end,delta_en,volt_start,volt_end,delta_volt,corr_ener
      integer :: size_l,size_r,size_c,size_lc,size_lcr,numatomic,num_threads,numfcore,numfvirt
      character(len=100) :: inputfile,outfile,option
      logical :: libint,gamess,rdm_flag,invert,doubles,currentflag,cisd_flag,qchem,hf_flag
      logical :: dft_flag,pyscf,maple,use_b0,molcas,write_ruqt_data
      integer :: i,j,k,counter,counter2,current_values,norb,numact,energy_val,volt_val,ioerror,numocc,numvirt
      character(len=40) :: ElectrodeType,CalcType,functional,inputcode,b0_type
      real(8) :: KT,current_con,Fermi_enl,Fermi_enR,localden_fermi,localden_fermi_l,localden_fermi_r,temp
      real(8) :: fermi_l,fermi_r
      complex(8), allocatable, dimension(:,:) :: test
      type(B1) :: B1data,l1data
      type(B2) :: B2data,l2data
      type(energr) :: G_S
      real(8) :: time_start,time_end
      integer :: state_num

      call cpu_time(time_start)
      invert=.true.
      currentflag=.false.     
      write(*,*) "Just getting started"
      Call Get_Command_Argument(1,inputfile)

     Call ReadInput(inputfile,norb,numfcore,numfvirt,numocc,numvirt,size_l,size_r,size_c,energy_start,energy_end,delta_en,volt_start,volt_end,delta_volt,inputcode,KT,ElectrodeType,Fermi_enl,Fermi_enR,CalcType,localden_fermi_l,localden_fermi_r,doubles,numatomic,functional,num_threads,use_b0,b0_type,write_ruqt_data,state_num)
     write(*,*) "Input File Read"
     write(*,*) "Using the following parameters for the transport calculation"
     write(*,*) "Number of OpenMP Threads:",num_threads
     write(*,*) "Number of  Molecular and Atomic Orbitals:",norb,numatomic
     write(*,*) "Number of Active Orbitals:",norb-numfcore-numfvirt
     write(*,*) "Number of Occupied Orbitals:",numocc
     write(*,*) "Number of Virtual Orbitals:",numvirt
     write(*,*) "Orbitals in left electrode:",size_l
     write(*,*) "Orbitals in device region:",size_c
     write(*,*) "Orbitals in right electrode:",size_r
     write(*,*) "Transmission Energy Window(eV) and dE:",energy_start,energy_end,delta_en
     write(*,*) "Voltage Window and dV:",volt_start,volt_end,delta_volt
     write(*,*) "Fermi Density:",localden_fermi_l,localden_fermi_r
     write(*,*) "KT:",KT
     write(*,*) "Transport Calculated for State ",state_num
      size_lc = size_l + size_c
      size_lcr = size_l + size_c + size_r
      libint=.false.
      Call Flag_set(inputcode,functional,cisd_flag,rdm_flag,hf_flag,dft_flag,qchem,gamess,pyscf,maple,molcas)

      if(qchem.eqv..true.) then
        write(*,*) 'Using Qchem data for this run'
        Call Get_HF_Qchem(inputfile,norb,H_Two,Smat)
       elseif(molcas.eqv..true.) then
        Call Get_HF_Molcas2("MolEl.dat",norb,H_Two,Smat,state_num)
       elseif(libint.eqv..true.) then
        Call Get_HF_libint(inputfile,norb,numact,H_one,Smat,mo_coeff,OneInts,TwoIntsCompact)
       elseif(gamess.eqv..true.) then
        Call Get_HF_GAMESS(inputfile,numatomic,H_Two,Smat,norb)
        write(*,*) 'Using GAMESS data for this run'
       elseif(pyscf.eqv..true.) then
        write(*,*) 'Using PySCF data for this run'
        Call Get_HF_PySCF(inputfile,numatomic,H_Two,Smat,norb)
        !Call Get_HF_Psi4()
       elseif(maple.eqv..true.) then
        write(*,*) "Using Maple+QuantumChemistry data for this run"
        Call Get_HF_PySCF(inputfile,numatomic,H_Two,Smat,norb)
       end if

      write(*,*) 'This run using:'
      if(rdm_flag.eqv..true.) then
        write(*,*) 'The Lehmann representation of a'

         if(doubles.eqv..true.) then
            write(*,*) 'p2-RDM Greens function'
           elseif(doubles.eqv..false.) then
            write(*,*) 'HF Greens function'
          end if

        elseif(qchem.eqv..true.) then
           write(*,*) 'QCHEM HF/DFT Greens function'

        elseif(molcas.eqv..true.) then
           write(*,*) 'Using Molcas Fock Matrix: FOCK_AO'
      
        elseif(pyscf.eqv..true.) then
           write(*,*) 'PYSCF HF/DFT Greens Function'

       end if

      !CALL OMP_SET_NUM_THREADS(num_threads)
       
!Here we want to parition the Smatrix and H matrix into electrodes and
!the device
      if(trim(ElectrodeType).eq."Metal_WBL") then
       write(*,*) 'Starting Metal WBL calculation'
       allocate(Sigma_l(1:size_c,1:size_c))
       allocate(Sigma_r(1:size_c,1:size_c))
  
       Call PartitionHS_MetalWBL(Smat,H_Two,size_l,size_r,size_c,size_lc,size_lcr,Smat_le,Smat_re,Smat_cen,H_Two_le,H_Two_re,H_Two_cen,H_Two_le_trans,H_Two_re_trans,write_ruqt_data,inputfile)

       write(*,*) 'Getting Electrodes'
       Call Electrodes_MetalWBL(Sigma_l,Sigma_r,Smat_re,Smat_le,H_Two_le,H_Two_re,localden_fermi_l,localden_fermi_r,size_c,size_lc,size_lcr,size_l,size_r,H_Two_le_trans,H_Two_re_trans,write_ruqt_data,inputfile)

       allocate(Gamma_L(1:size_c,1:size_c))
       allocate(Gamma_R(1:size_c,1:size_c))
       Gamma_L=-(DIMAG(Sigma_L)-DIMAG(adjoint(Sigma_L,size_c)))!,1E-12,8)
       Gamma_R=-(DIMAG(Sigma_R)-DIMAG(adjoint(Sigma_R,size_c)))!,1E-12,8)



      if(trim(CalcType).eq."current".or.trim(CalcType).eq."Current") then
       write(*,*) "Starting Current Calculation"
       volt_val = abs(int((volt_end - volt_start)/delta_volt)+1)
       energy_val = abs(int((energy_start-energy_end)/delta_en)+1)
       current_values = volt_val 
       allocate(gfc_r(1:size_c,1:size_c))
       allocate(gfc_a(1:size_c,1:size_c))
       allocate(current(1:current_values))
       allocate(voltage(1:current_values))
       allocate(transm(1:energy_val))
       allocate(transm_curr(1:energy_val))
       allocate(current_temp(1:size_c,1:size_c))
       allocate(energy_list(1:energy_val))
       transm=0
       transm_curr=0
       current_temp=0

     
       energy = energy_start
       current_con = 2*1.6021766E-19*(4.135667E-15)**(-1)
       counter = 1
      do k=1,energy_val
          transm(k) = 0
          current_temp=0
          energy = energy_start + (k-1)*delta_en
          energy_list(k) = energy
         if((hf_flag.eqv..true.).or.(dft_flag.eqv..true.)) then
          gfc_r = 0
          gfc_a = 0
          gfc_r = energy*Smat_cen-H_Two_cen
          gfc_r = gfc_r - Sigma_l - Sigma_r
          gfc_r = inv(gfc_r)
          gfc_a = adjoint(gfc_r,size_c)
         else if((cisd_flag.eqv..true.).or.(rdm_flag.eqv..true.)) then
            gfc_r = 0
            gfc_a = 0
            Call Build_G_SD_Invert(gfc_r,Sigma_l,Sigma_r,energy,size_l,size_c,size_lc,size_lcr,norb,inputfile,numocc,numvirt,counter,B1data,B2data,mo_ener,mo_coeff,mo_coeff2,doubles,currentflag,energy_val,k,G_S,corr_ener,numatomic,B0_coeff,use_b0,gamess,maple,numfcore,numfvirt,b0_type)
            gfc_a = adjoint(gfc_r,size_c)
            counter=2
         end if
         current_temp = matmul_zgemm(Gamma_R,gfc_a)
         current_temp = matmul_zgemm(gfc_r,current_temp)
         current_temp = matmul_zgemm(Gamma_L,current_temp)
         do i=1,size_c
           transm(k) = transm(k) + real(current_temp(i,i))
         end do
      end do
    
      do j=1,energy_val
         write(*,*) 'Transm vs Energy curve',real(energy_list(j)),real(transm(j))
       end do

       outfile = trim(inputfile) // ".negf_dat"
       open(unit=7,file=outfile,action='write',iostat=ioerror)

       write(7,*) energy_val
       do j=1,energy_val
         write(7,*) real(energy_list(j)),real(transm(j))
       end do

      ! close(7)

     currentflag=.true.

    do j=1,volt_val
       temp = (j-1)*delta_volt + volt_start
       voltage(j) = temp
       current(j) = 0
       transm_curr=0
      do k=1,energy_val
        !  transm(k) = 0
          energy = energy_start + (k-1)*delta_en
!          current_temp = 0
!         if((hf_flag.eqv..true.).or.(dft_flag.eqv..true.)) then
!          gfc_r = 0
!          gfc_a = 0
!          gfc_r = energy*Smat_cen-H_Two_cen
!          gfc_r = gfc_r - Sigma_l - Sigma_r
!          gfc_r = inv(gfc_r)
!          gfc_a = adjoint(gfc_r,size_c)
!         else if((cisd_flag.eqv..true.).or.(rdm_flag.eqv..true.)) then
!            gfc_r = 0
!            gfc_a = 0
!            Call Build_G_SD_Invert(gfc_r,Sigma_l,Sigma_r,energy,size_l,size_c,size_lc,size_lcr,norb,inputfile,numocc,numvirt,counter,B1data,B2data,mo_ener,mo_coeff,mo_coeff2,doubles,currentflag,energy_val,k,G_S,corr_ener,numatomic,B0_coeff,use_b0,gamess,maple,numfcore,numfvirt,b0_type)
!            gfc_a = adjoint(gfc_r,size_c)
!            counter=2
!
!
!          end if


!          current_temp = matmul_zgemm(Gamma_R,gfc_a)
!          current_temp = matmul_zgemm(gfc_r,current_temp)
!          current_temp = matmul_zgemm(Gamma_L,current_temp)
          !do i=1,size_c
          ! transm(k) = transm(k) + real(current_temp(i,i))
          !end do
          fermi_l=fermi_function(energy-temp*0.5,fermi_enL,KT)
          fermi_r=fermi_function(energy+temp*0.5,fermi_enR,KT)
          !write(*,*) k,fermi_l,fermi_r
          transm_curr(k) = transm(k)*(fermi_l-fermi_r)
         end do
         !write(*,*) temp,KT,transm_curr(1)
        do k=1,energy_val
         current(j) = current(j) + delta_en*current_con*transm_curr(k)
        end do
        !write(*,*) 'Done with current at voltage:',real(voltage(j))
       end do
       
       do j=1,volt_val
         write(*,*) 'IV curve',real(voltage(j)),real(current(j))
       end do

       !outfile = trim(inputfile) // ".negf_dat"
       !open(unit=7,file=outfile,action='append',iostat=ioerror)

       write(7,*) volt_val
       do j=1,volt_val
         write(7,*) real(voltage(j)),real(current(j))
       end do

       close(7)

      elseif(trim(CalcType).eq."Transmission".or.trim(CalcType).eq."transmission") then

       write(*,*) "Starting Transmission Calculation at Energy"
       energy_val = abs(int((energy_start-energy_end)/delta_en))
       allocate(gfc_r(1:size_c,1:size_c))
       allocate(gfc_a(1:size_c,1:size_c))
       allocate(transm(1:energy_val))
       allocate(energy_list(1:energy_val))
       allocate(current_temp(1:size_c,1:size_c))


       energy = energy_start
       counter=1
      do k=1,energy_val
          transm(k) = 0
          current_temp = 0
          energy = energy_start + (k-1)*delta_en
          energy_list(k) = energy

         if((hf_flag.eqv..true.).or.(dft_flag.eqv..true.)) then
            gfc_r = 0
            gfc_a = 0
            gfc_r = energy*Smat_cen-H_Two_cen
            gfc_r = gfc_r - Sigma_l - Sigma_r
            gfc_r = inv(gfc_r)
            gfc_a = adjoint(gfc_r,size_c)

         else if((cisd_flag.eqv..true.).or.(rdm_flag.eqv..true.)) then
            gfc_r = 0
            gfc_a = 0
            Call Build_G_SD_Invert(gfc_r,Sigma_l,Sigma_r,energy,size_l,size_c,size_lc,size_lcr,norb,inputfile,numocc,numvirt,counter,B1data,B2data,mo_ener,mo_coeff,mo_coeff2,doubles,currentflag,energy_val,k,G_S,corr_ener,numatomic,B0_coeff,use_b0,gamess,maple,numfcore,numfvirt,b0_type)
            gfc_a = adjoint(gfc_r,size_c)
            counter=2


         end if

         current_temp = matmul_zgemm(Gamma_L,gfc_r)
         current_temp = matmul_zgemm(current_temp,Gamma_R)
         current_temp = matmul_zgemm(current_temp,gfc_a)
         do i=1,size_c
            transm(k) = transm(k) + real(current_temp(i,i))
          end do
          write(*,*) 'Done with transmission function at energy:',transm(k)
        end do

       do j=1,energy_val
         write(*,*) 'Transm vs Energy curve',real(energy_list(j)),real(transm(j))
       end do
 
       outfile = trim(inputfile) // ".negf_dat"
       open(unit=7,file=outfile,action='write',iostat=ioerror)
    
       write(7,*) energy_val
       do j=1,energy_val
         write(7,*) real(energy_list(j)),real(transm(j))
       end do

       close(7)
      end if

      else if(trim(ElectrodeType).eq."Molecule_WBL") then
       write(*,*) 'Using Molecule WBL Electrodes'
       write(*,*) "***This option only supports DFT/HF calcs and is outdated/buggy***"
       write(*,*) "***Slated for removal and likely to be incorrect***"
       write(*,*) "***Use at own risk***"
       allocate(Sigma_l(1:norb,1:norb))
       allocate(Sigma_r(1:norb,1:norb))
       Sigma_l = 0
       Sigma_r = 0
       write(*,*) 'Calculating Coupling'
       Call Calculate_Coupling_MoleculeWBL(Coupling_R,Coupling_L,localden_fermi)
       write(*,*) 'Calculating Electrodes'
       Call Electrodes_MoleculeWBL(Sigma_l,Sigma_r,Smat,Coupling_R,Coupling_L,size_c,size_lc,size_lcr)
       allocate(Gamma_L(1:size_c,1:size_c))
       allocate(Gamma_R(1:size_c,1:size_c))
       Gamma_L=DIMAG(Sigma_L)-DIMAG(adjoint(Sigma_L,norb))
       Gamma_R=DIMAG(Sigma_R)-DIMAG(adjoint(Sigma_R,norb))

      if(trim(CalcType).eq."current".or.trim(CalcType).eq."Current") then
       write(*,*) "Starting Current Calculation"
       volt_val = abs(int((volt_end - volt_start)/delta_volt))
       energy_val = abs(int((energy_start-energy_end)/delta_en))
       current_values = volt_val 
       allocate(gfc_r(1:norb,1:norb))
       allocate(gfc_a(1:norb,1:norb))
       allocate(current(1:current_values))
       allocate(voltage(1:current_values))
       allocate(transm(1:energy_val))
       allocate(current_temp(1:norb,1:norb))


       energy = energy_start
       current_con = 1.6021766E-19*4.135667E-15**(-1)
       counter = 1

        current_temp = 0

        gfc_r = 0
        gfc_a = 0
        gfc_r = energy*Smat-H_Two
        gfc_r = gfc_r - Sigma_l - Sigma_r
        gfc_r = inv(gfc_r)
        gfc_a = adjoint(gfc_r,norb)

        current_temp = matmul_zgemm(Gamma_R,gfc_a)
        current_temp = matmul_zgemm(gfc_r,current_temp)
        current_temp = matmul_zgemm(Gamma_L,current_temp)


    do j=1,volt_val
       temp = (j-1)*delta_volt + volt_start
       voltage(j) = temp
       current(j) = 0
      do k=1,energy_val
          transm(k) = 0
          energy = energy_start + (k-1)*delta_en
        !  current_temp = 0

        !  gfc_r = 0
        !  gfc_a = 0
        !  gfc_r = energy*Smat-H_Two
        !  gfc_r = gfc_r - Sigma_l - Sigma_r
        !  gfc_r = inv(gfc_r)
        !  gfc_a = adjoint(gfc_r,norb)

         ! current_temp = matmul_zgemm(Gamma_R,gfc_a)
         ! current_temp = matmul_zgemm(gfc_r,current_temp)
         ! current_temp = matmul_zgemm(Gamma_L,current_temp)

          do i=1,norb
           transm(k) = transm(k) + real(current_temp(i,i))
          end do
          transm(k) = transm(k)*(fermi_function(energy+temp*0.5,fermi_enL,KT)-fermi_function(energy-temp*0.5,fermi_enR,KT))
          
         end do
        do k=1,energy_val
         current(j) = current(j) + delta_en*current_con*transm(k)
        end do
       end do

       do j=1,volt_val
         write(*,*) 'IV curve',real(voltage(j)),real(current(j))
       end do

       outfile = trim(inputfile) // ".negf_dat"
       open(unit=7,file=outfile,action='write',iostat=ioerror)
       write(7,*) volt_val
       do j=1,volt_val
         write(7,*) real(voltage(j)),real(current(j))
       end do

       close(7)

      elseif(trim(CalcType).eq."Transmission".or.trim(CalcType).eq."transmission") then

       write(*,*) "Starting Transmission Calculation at Energy"
       energy_val = abs(int((energy_start-energy_end)/delta_en))
       allocate(gfc_r(1:norb,1:norb))
       allocate(gfc_a(1:norb,1:norb))
       allocate(transm(1:energy_val))
       allocate(energy_list(1:energy_val))
       allocate(current_temp(1:norb,1:norb))


       energy = energy_start
       temp = volt_start
      do k=1,energy_val
          transm(k) = 0
          energy = energy_start + (k-1)*delta_en
          energy_list(k) = energy

          gfc_r = 0
          gfc_a = 0
          gfc_r = energy*Smat-H_Two
          gfc_r = gfc_r - Sigma_l - Sigma_r
          gfc_r = inv(gfc_r)
          gfc_a = adjoint(gfc_r,norb)


          current_temp = matmul_zgemm(Gamma_R,gfc_a)
          current_temp = matmul_zgemm(gfc_r,current_temp)
          current_temp = matmul_zgemm(Gamma_L,current_temp)

          do i=1,norb
           transm(k) = transm(k) + real(current_temp(i,i))
          end do
         end do

       do j=1,energy_val
         write(*,*) 'Transm vs Energy curve',real(energy_list(j)),real(transm(j))
       end do

       outfile = trim(inputfile) // ".dat"
       open(unit=7,file=outfile,action='write',iostat=ioerror)

       do j=1,energy_val
         write(7,*) real(energy_list(j)),real(transm(j))
       end do

       close(7)
      end if
     end if

      call cpu_time(time_end)
      write(*,*) 'CPU Time:',time_end-time_start

      end program

      function adjoint(A,norb)
      implicit none
 
      complex(8),allocatable,dimension(:,:) :: A,adjoint
      integer :: i,j,norb

       allocate(adjoint(1:norb,1:norb))

       do i=1,norb
        do j=1,norb
         adjoint(j,i) =  DCONJG(A(i,j))
        end do
       end do
     
      end function


      subroutine Get_HF_QChem(inputfile,norb,H_two,Smat)
      !Use InterfaceMod 
      implicit none

      character(len=100) :: inputfile,datafile
      real(8),allocatable,dimension(:,:) :: H_Two,Smat
      integer :: norb,ioerror,i,j

      20 format(A)
      !write(*,*) inputfile
      allocate(H_two(1:norb,1:norb))
      allocate(Smat(1:norb,1:norb))
      H_two=0
      Smat=0

      datafile = trim(inputfile) // "_Smat"
      open(unit=2,file=datafile,action='READ', iostat = ioerror)

      do i=1,norb
       do j=1,norb
        read(2,*) Smat(i,j)
       end do
      end do
      close(2)

       write(*,*) 'Start Htwo'
      datafile = trim(inputfile) // "_Htwo"
      open(unit=3,file=datafile,action='READ', iostat = ioerror)

      do i=1,norb
       do j=1,norb
        read(3,*) H_Two(i,j)
       end do
      end do
      close(3)
      write(*,*) 'Done Getting HF Values'
      end subroutine

      subroutine Get_HF_Molcas(inputfile,norb,H_two,Smat,state_num)
      !Use InterfaceMod
      implicit none

      character(len=100) :: inputfile,datafile
      real(8),allocatable,dimension(:,:) :: H_Two,Smat
      integer :: norb,ioerror,i,j,x,y
      real(8) :: readtemp
      integer :: state_num
      character(len=100) :: state_char
      20 format(A)
      !write(*,*) inputfile
      allocate(H_two(1:norb,1:norb))
      allocate(Smat(1:norb,1:norb))
      H_two=0
      Smat=0

      datafile = "Overlap"
      open(unit=2,file=datafile,action='READ', iostat = ioerror)

      do i=1,norb
       do j=i,norb
        read(2,*) x,y,readtemp
         Smat(x,y)=readtemp
         Smat(y,x)=readtemp
       end do
      end do
      close(2)

       write(*,*) 'Start Fock Matrix'
      
      write(*,*) state_num
      if(state_num.lt.10) then
       write(state_char,"(I1)") state_num
      else
       write(state_char,"(I2)") state_num
      end if

      datafile = "FOCK_AO_"//trim(state_char)
      write(*,*) datafile
      open(unit=3,file=datafile,action='READ', iostat = ioerror)

      do i=1,norb
       do j=1,norb
        read(3,*) x,y,readtemp
        H_Two(x,y)=readtemp
       end do
      end do
      close(3)
      write(*,*) 'Done Getting HF Values'
      end subroutine

      subroutine Get_HF_Molcas2(datafile,norb,H_two,Smat,state_num)
      !Use InterfaceMod
      implicit none

      character(len=9) :: datafile
      real(8),allocatable,dimension(:,:) :: H_Two,Smat
      integer :: norb,ioerror,i,j,x,y
      real(8) :: readtemp
      integer :: state_num,num_states,norb_2,nelec_2,actorb,actel
      character(len=100) :: state_char,readtemp_str
      20 format(A)
      !write(*,*) inputfile
      allocate(H_two(1:norb,1:norb))
      allocate(Smat(1:norb,1:norb))
      H_two=0
      Smat=0

      open(unit=2,file=datafile,action='READ', iostat = ioerror)

      read(2,*) readtemp_str
      read(2,*) num_states,norb_2,nelec_2,actorb,actel
      read(2,*) readtemp_str

      if(norb_2.ne.norb) then
        write(*,*) "Your orbital count is incorrect. Please check your input and MolEl.dat files"
        stop
       end if

      do i=1,norb
       do j=i,norb
        read(2,*) x,y,readtemp
         Smat(x,y)=readtemp
         Smat(y,x)=readtemp
       end do
      end do

      do while (trim(readtemp_str).ne."Effective")
       read(2,*) readtemp_str
      end do
      
      x=0
      do while (trim(readtemp_str).ne."State".and.x.ne.state_num)
       read(2,*) readtemp_str,x
      end do

      write(*,*) "Reading Fock Matrix for ",state_num

      do i=1,norb
       do j=1,norb
        read(2,*) x,y,readtemp
        H_Two(x,y)=readtemp
       end do
      end do
      close(2)
      write(*,*) 'Done Getting HF Values'
      end subroutine

 
      
      Subroutine Get_HF_GAMESS(inputfile,numatomic,H_two,Smat,norb)
      use FunctionMod
      implicit none

      character(len=100) :: inputfile,datafile
      real(8),allocatable,dimension(:,:) :: H_Two,Smat,mo_data,mo_data2,ECP_m,ECP_a,ECP_temp
      integer :: numatomic,ioerror,i,j,norb

      20 format(A)
      allocate(H_two(1:numatomic,1:numatomic))
      allocate(Smat(1:numatomic,1:numatomic))
      H_two=0
      Smat=0

      datafile = trim(inputfile) // "_Smat"
      open(unit=2,file=datafile,action='READ', iostat = ioerror)

      do i=1,numatomic
       do j=1,numatomic
        read(2,*) Smat(i,j)
       end do
      end do
      close(2)

       write(*,*) 'Start Htwo'
      datafile = trim(inputfile) // "_Htwo"
      open(unit=3,file=datafile,action='READ', iostat = ioerror)

      do i=1,numatomic
       do j=1,numatomic
        read(3,*) H_Two(i,j)
       end do
      end do
      close(3)

      write(*,*) 'Get ECP'
      datafile = trim(inputfile) // "_ecp"
      open(unit=4,file=datafile,status='OLD',action='READ', iostat = ioerror)

     if(ioerror.eq.0) then
     write(*,*) 'ECP found'
     allocate(ECP_m(1:norb,1:norb))
     allocate(ECP_temp(1:numatomic,1:norb))
     allocate(ECP_a(1:numatomic,1:numatomic))
     allocate(mo_data(1:numatomic,1:norb))
     allocate(mo_data2(1:norb,1:numatomic))
 
      do i=1,norb
       do j=1,norb
        read(4,*) ECP_m(i,j)
       end do
      end do

      close(5)

      datafile = trim(inputfile) // ".mo_dat"
      open(unit=6,file=datafile,action='READ', iostat = ioerror)

      do i=1,numatomic
       do j=1,norb
        read(6,*) mo_data(i,j)
       end do
      end do
      close(6)
      mo_data2 = transpose(mo_data)

      ECP_temp=matmul_dgemm(mo_data,ECP_m)
      ECP_a=matmul_dgemm(ECP_temp,mo_data2)
      H_two=H_two-ECP_a
      else
       write(*,*) "No removal of ECP"
      end if
      write(*,*) 'H_Two',size(H_Two)
      write(*,*) 'Smat',size(Smat)
      !stop
      write(*,*) 'Done Getting HF Values from GAMESS'
      end subroutine

      Subroutine Get_HF_PySCF(inputfile,numatomic,H_two,Smat,norb)
      !Use InterfaceMod
      use FunctionMod
      implicit none

      character(len=100) :: inputfile,datafile,readtemp
      real(8),allocatable,dimension(:,:) :: H_Two,Smat,mo_data,mo_data2,ECP_m,ECP_a,ECP_temp
      integer :: numatomic,ioerror,i,j,norb,tempx,tempy
      real(8) :: tempval

      20 format(A)
      allocate(H_two(1:numatomic,1:numatomic))
      allocate(Smat(1:numatomic,1:numatomic))
      H_two=0
      Smat=0

      datafile = trim(inputfile) // ".scf_dat"
      open(unit=2,file=datafile,action='READ', iostat = ioerror)

      do while(trim(readtemp)/="Overlap Matrix")
       read(2,'(A)') readtemp
       end do
      do i=1,numatomic
       do j=1,numatomic
        read(2,*) tempx,tempy,tempval
        Smat(tempx,tempy)=tempval
       end do
      end do

      do while(trim(readtemp)/="Fock Matrix")
        read(2,'(A)') readtemp
       end do

      do i=1,numatomic
       do j=1,numatomic
        read(2,*) tempx,tempy,tempval
        H_Two(tempx,tempy)=tempval
       end do
      end do
      close(2)
     end subroutine


      subroutine Get_HF_libint(inputfile,norb,numact,H_one,Smat,mo_coeff,OneInts,TwoIntsCompact)
      use TypeMod
      use FunctionMod
      implicit none

      character(len=100) :: inputfile,datafile
      real(8),allocatable,dimension(:,:) :: H_one,H_Two,Smat,OneInts,mo_coeff
      real(8),allocatable,dimension(:) :: TwoIntsCompact
      integer :: norb,numact,ioerror,i,j,k,l
      integer(8) :: index1,index2,compindex1
      integer, dimension(1:4) :: orbind
      real(8) :: integral  

      allocate(H_one(1:norb,1:norb))
      allocate(Smat(1:norb,1:norb))
      allocate(mo_coeff(1:norb,1:norb))
      allocate(OneInts(1:numact,1:numact))
      allocate(TwoIntsCompact(1:numact*(numact+1)/2*(numact*(numact+1)/2+1)/2))

      20 format(A)
      datafile = inputfile // ".Hone"
      open(unit=1,file=datafile,action='READ', iostat = ioerror)

      do i=1,norb
       do j=1,norb
        read(1,20) H_one(i,j)
       end do
      end do
      close(1)

      datafile = inputfile // ".Smat"
      open(unit=2,file=datafile,action='READ', iostat = ioerror)

      do i=1,norb
       do j=1,norb
        read(2,20) Smat(i,j)
       end do
      end do
      close(2)


      open(unit=6,file=datafile,action='READ',iostat= ioerror)

      do i=1,norb
       do j=1,norb
        read(6,20) mo_coeff(i,j)
       end do
      end do

 
      datafile = inputfile // ".OneInts"
      open(unit=4,file=datafile,action='READ', iostat = ioerror)

      do i=1,norb
       do j=1,norb
        read(4,20) OneInts(i,j)
       end do
      end do
      close(4)

      datafile = inputfile // ".TwoInts"
      open(unit=5,file=datafile,action='READ', iostat = ioerror)

       do while(orbind(1).ne.0)

         read(5,*) orbind(1),orbind(2),orbind(3),orbind(4),integral
                     i = orbind(1)
                     k = orbind(2)
                     j = orbind(3)
                     l = orbind(4)
                     index1 = FirstIndex(i,k)
                     index2 = FirstIndex(j,l)
                     compindex1 = CompositeIndex(index1,index2)

                    TwoIntsCompact(compindex1) = integral

        end do


      end subroutine


      subroutine IntTransform(TwoIntsCompact,mo_coeff)
      use TypeMod
      implicit none

      real(8), allocatable,dimension(:,:) :: mo_coeff
      real(8), allocatable,dimension(:) :: TwoIntsCompact
      !There is nothing here yet!      

      end subroutine

      subroutine Calculate_Coupling_MoleculeWBL(Coupling_R,Coupling_L,localden_fermi)
      implicit none    

      real(8) :: Coupling_R,Coupling_L,localden_fermi

      Coupling_R =  0.2
      Coupling_L =  0.2
      end subroutine


      subroutine Electrodes_MoleculeWBL(Sigma_l,Sigma_r,Smat,Coupling_R,Coupling_L,size_c,size_lc,size_lcr)
      implicit none
      integer :: size_lc,size_c,size_lcr,i,j
      real(8) :: Coupling_R,Coupling_L,temp
      complex(8), allocatable, dimension(:,:) :: Sigma_L,Sigma_R
      real(8), allocatable, dimension(:,:) :: Smat_LE,Smat_RE,Smat

      do i = size_c+1,size_lc
       do j = 1,size_c
      Sigma_L(i,j) = CMPLX(0,-0.5*Coupling_L*Smat(i,j),8)
      Sigma_L(j,i) = CMPLX(0,-0.5*Coupling_L*Smat(j,i),8)
       end do
      end do
      do i=size_lc+1,size_lcr
       do j= 1,size_c
      Sigma_R(i,j) = CMPLX(0,-0.5*Coupling_R*Smat(i,j),8)
      Sigma_R(j,i) = CMPLX(0,-0.5*Coupling_R*Smat(j,i),8)
       end do
      end do
      end subroutine


      subroutine Electrodes_MetalWBL(Sigma_l,Sigma_r,Smat_re,Smat_le,H_Two_le,H_Two_re,localden_fermi_l,localden_fermi_r,size_c,size_lc,size_lcr,size_l,size_r,H_Two_le_trans,H_Two_re_trans,write_ruqt_data,inputfile)
      use FunctionMod
      implicit none
      integer :: size_lc,size_c,size_lcr,i,j,size_l,size_r,ioerror
      real(8) :: Coupling_R,Coupling_L,temp,localden_fermi_l,localden_fermi_r,pi
      complex(8), allocatable, dimension(:,:) :: Sigma_L,Sigma_R
      real(8), allocatable, dimension(:,:) :: Smat_LE,Smat_RE,Smat,H_Two_re,H_Two_le,Sigma_L_temp,Sigma_R_temp,Sigma_R_temp2,Sigma_L_temp2,sigma_r_temp3,sigma_l_temp3,H_Two_le_trans,H_Two_re_trans
      logical :: write_ruqt_data
      character(len=100) :: inputfile,outfile

      pi = 3.14159265359

      allocate(sigma_l_temp(1:size_l,1:size_l))
      allocate(sigma_l_temp2(1:size_c,1:size_l))
      allocate(sigma_l_temp3(1:size_c,1:size_c))

       write(*,*) 'alloc l done'
      Sigma_L_temp = -pi*localden_fermi_l*Smat_le

       call matmul_dgemm2(H_Two_le_trans,Sigma_L_temp,Sigma_L_temp2)
       call matmul_dgemm2(Sigma_L_temp2,H_Two_le,Sigma_L_temp3)
      Sigma_L = CMPLX(0,sigma_l_temp3,8)

      write(*,*) 'l done'
      deallocate(sigma_l_temp)
      deallocate(sigma_l_temp2)
      deallocate(sigma_l_temp3)

      allocate(sigma_r_temp(1:size_r,1:size_r))
      allocate(sigma_r_temp2(1:size_c,1:size_r))
      allocate(sigma_r_temp3(1:size_c,1:size_c))

      Sigma_R_temp = -pi*localden_fermi_r*Smat_re

       call matmul_dgemm2(H_Two_re_trans,Sigma_R_temp,Sigma_R_temp2)
       call matmul_dgemm2(Sigma_R_temp2,H_Two_re,Sigma_R_temp3)
      Sigma_R = CMPLX(0,Sigma_R_temp3,8)

      deallocate(Sigma_R_temp)
      deallocate(Sigma_R_temp2)
      deallocate(Sigma_R_temp3)

      if(write_ruqt_data) then
        outfile = trim(inputfile) // ".Sigma"
        open(unit=8,file=outfile,action='write',iostat=ioerror)

        write(8,*) size_l,size_c,size_r
        write(8,*) "RUQT Sigma Matrices"
        do j=1,size_c
          do i=1,size_c
             write(8,*) j,i,Sigma_R(j,i)
             end do
            end do
         do j=1,size_c
           do i=1,size_c
              write(8,*) j,i,Sigma_R(j,i)
             end do
            end do
       close(8)
      end if
      write(*,*) 'Sigma Calculated'
      end subroutine

              subroutine PartitionHS_MetalWBL(Smat,H_Two,size_l,size_r,size_c,size_lc,size_lcr,Smat_le,Smat_re,Smat_cen,H_Two_le,H_Two_re,H_Two_cen,H_Two_le_trans,H_Two_re_trans,write_ruqt_data,inputfile)
              implicit none

              real(8), allocatable, dimension(:,:) :: Smat,Smat_le,Smat_re,Smat_cen,H_Two,H_Two_cen,H_Two_re,H_Two_le,temp,H_Two_le_trans,H_Two_re_trans
              real(8), allocatable, dimension(:) :: eigen,work
              integer :: size_l,size_c,size_r,size_lc,size_lcr,info,lwork,i,j,ioerror
              logical :: write_ruqt_data
              character(len=100) :: inputfile,outfile

              allocate(Smat_le(1:size_l,1:size_l))
              allocate(Smat_re(1:size_r,1:size_r))
              allocate(Smat_cen(1:size_c,1:size_c))
              Smat_le(1:size_l,1:size_l)=Smat(1:size_l,1:size_l)
              Smat_re(1:size_r,1:size_r)=Smat(size_lc+1:size_lcr,size_lc+1:size_lcr)
              Smat_cen(1:size_c,1:size_c)=Smat(size_l+1:size_lc,size_l+1:size_lc)

               write(*,*) 'Done partitioning S Matrix'
              deallocate(Smat)

              allocate(H_Two_le(1:size_l,1:size_c))
              allocate(H_Two_re(1:size_r,1:size_c))
              allocate(H_Two_cen(1:size_c,1:size_c))
              allocate(H_Two_le_trans(1:size_c,1:size_l))
              allocate(H_Two_re_trans(1:size_c,1:size_r))

              H_Two_le(1:size_l,1:size_c)=27.2114*H_Two(1:size_l,size_l+1:size_lc)
              H_Two_re(1:size_r,1:size_c)=27.2114*H_Two(size_lc+1:size_lcr,size_l+1:size_lc)
              H_Two_cen(1:size_c,1:size_c)=27.2114*H_Two(size_l+1:size_lc,size_l+1:size_lc)
              H_Two_le_trans(1:size_c,1:size_l) = 27.2114*H_Two(size_l+1:size_lc,1:size_l) !transpose(H_Two_le)
              H_Two_re_trans(1:size_c,1:size_r) = 27.2114*H_Two(size_l+1:size_lc,size_lc+1:size_lcr)!transpose(H_Two_re)
              write(*,*) 'Done partitioning Fock Matrix'
              deallocate(H_Two) 
              if(write_ruqt_data) then
               outfile = trim(inputfile) // ".partdat"
               open(unit=8,file=outfile,action='write',iostat=ioerror)


               write(8,*) size_l,size_c,size_r
               write(8,*) "RUQT Overlap Matrices"
               do j=1,size_l
                 do i=1,size_l
                  write(8,*) j,i,Smat_le(j,i)
                end do
               end do
               do j=1,size_c
                do i=1,size_c
                 write(8,*) j,i,Smat_cen(j,i)
                end do
               end do
               do j=1,size_r
                do i=1,size_r
                 write(8,*) j,i, Smat_re(j,i)
                end do
               end do

               write(8,*) "RUQT Fock Matrices"
               do j=1,size_l
                do i=1,size_c
                  write(8,*) j,i,H_Two_le(j,i)/27.2114
                end do
               end do
               do j=1,size_r
                do i=1,size_c
                 write(8,*) j,i,H_Two_re(j,i)/27.2114
                end do
               end do
               do j=1,size_c
                do i=1,size_c
                 write(8,*) j,i,H_Two_cen(j,i)/27.2114
                end do
               end do

               close(8)
              end if
              end subroutine


              function fermi_function(energy,fermi_energy,KT)
              implicit none
              real(8) :: energydiff,KT,energy,fermi_energy,fermi_function

              energydiff = energy - fermi_energy
              fermi_function = 1.00/(exp(energydiff/KT)+1)

              end function

              function inv(A) result(Ainv)
               implicit none
               complex(8),allocatable, dimension(:,:), intent(in) :: A
               complex(8), allocatable, dimension(:,:) :: Ainv

               complex(8),allocatable, dimension(:) :: work  ! work array for LAPACK
               integer,allocatable, dimension(:) :: ipiv   ! pivot indices
               integer :: n, info

               allocate(Ainv(1:size(A,1),1:size(A,2)))
               allocate(work(1:size(A,1)))
               allocate(ipiv(1:size(A,1)))

               Ainv = A
               n = size(A,1)

               call ZGETRF(n, n, Ainv, n, ipiv, info)

               if (info /= 0) then
                 write(*,*) info!,A
                 stop 'Matrix is numerically singular!'
                end if

               call ZGETRI(n, Ainv, n, ipiv, work, n, info)

               if (info /= 0) then
                 write(*,*) info
                 stop 'Matrix inversion failed!'
                end if

                deallocate(work)
                deallocate(ipiv)
               end function inv     


              function inv_real(A) result(Ainv)
               implicit none
               real(8),allocatable, dimension(:,:), intent(in) :: A
               real(8), allocatable, dimension(:,:) :: Ainv

               real(8),allocatable, dimension(:) :: work  ! work array for LAPACK
               integer,allocatable, dimension(:) :: ipiv   ! pivot indices
               integer :: n, info


               allocate(Ainv(1:size(A,1),1:size(A,2)))
               allocate(work(1:size(A,1)))
               allocate(ipiv(1:size(A,1)))

               Ainv = A
               n = size(A,1)

               call DGETRF(n, n, Ainv, n, ipiv, info)

               if (info /= 0) then
                 write(*,*) info!,A
                 stop 'Matrix is numerically singular!'
                end if

               call DGETRI(n, Ainv, n, ipiv, work, n, info)

               if (info /= 0) then
                 write(*,*) info
                 stop 'Matrix inversion failed!'
                end if

                deallocate(work)
                deallocate(ipiv)
               end function inv_real

               function FirstIndex(i,k)

               Implicit None
               integer :: i,k
               integer(8) :: FirstIndex

               if(i.lt.k) then
                 FirstIndex = (k-1)*k/2 + i
               else
                 FirstIndex = (i-1)*i/2 + k
               end if

               end function FirstIndex


      function matmul_zgemm(leftmatrix,rightmatrix)
      Implicit NONE
      complex(8), allocatable, dimension(:,:) :: leftmatrix,rightmatrix
      complex(8), allocatable, dimension(:,:) :: matmul_zgemm
      integer :: lmr, lmc
      integer :: rmr, rmc

      lmr = size(leftmatrix,1)  !left matrix row size
      lmc = size(leftmatrix,2)  !left matrix col size
      rmr = size(rightmatrix,1) !right matrix row size
      rmc = size(rightmatrix,2) !right matrix col size

      allocate(matmul_zgemm(1:lmr,1:rmc))
      matmul_zgemm(1:lmr,1:rmc) = 0.0
      if(lmr.ne.0.and.lmc.ne.0.and.rmr.ne.0.and.rmc.ne.0)  then
      Call ZGEMM('N','N',lmr,rmc,lmc,1.d0,leftmatrix,lmr,rightmatrix,rmr,0.d0,matmul_zgemm,lmr)
      end if
      end function matmul_zgemm

      function matmul_dgemm(leftmatrix,rightmatrix)
      Implicit NONE
      real(8), allocatable, dimension(:,:) :: leftmatrix,rightmatrix
      real(8), allocatable, dimension(:,:) :: matmul_dgemm
      integer :: lmr, lmc
      integer :: rmr, rmc

      lmr = size(leftmatrix,1)  !left matrix row size
      lmc = size(leftmatrix,2)  !left matrix col size
      rmr = size(rightmatrix,1) !right matrix row size
      rmc = size(rightmatrix,2) !right matrix col size

      matmul_dgemm = 0.0
      if(lmr.ne.0.and.lmc.ne.0.and.rmr.ne.0.and.rmc.ne.0)  then
      Call DGEMM('N','N',lmr,rmc,lmc,1.d0,leftmatrix,lmr,rightmatrix,rmr,0.d0,matmul_dgemm,lmr)
      end if
      end function matmul_dgemm

      subroutine matmul_dgemm2(leftmatrix,rightmatrix,outmat)
      Implicit NONE
      real(8), allocatable, dimension(:,:) :: leftmatrix,rightmatrix
      real(8), allocatable, dimension(:,:) :: outmat
      integer :: lmr, lmc,cmc
      integer :: rmr, rmc,cmr

      lmr = size(leftmatrix,1)  !left matrix row size
      lmc = size(leftmatrix,2)  !left matrix col size
      rmr = size(rightmatrix,1) !right matrix row size
      rmc = size(rightmatrix,2) !right matrix col size
      cmr= size(outmat,1)
      cmc=size(outmat,2)
      outmat = 0.0
      if(lmr.ne.0.and.lmc.ne.0.and.rmr.ne.0.and.rmc.ne.0)  then
      Call DGEMM("N","N",lmr,rmc,lmc,1.d0,leftmatrix,lmr,rightmatrix,rmr,0.d0,outmat,cmr)
      end if
      end subroutine matmul_dgemm2


              function CompositeIndex(ik,jl)
              Implicit NONE

              integer(8) :: ik, jl
              integer(8) :: CompositeIndex

               CompositeIndex = 0

              if(ik.lt.jl) then
                CompositeIndex =  (jl-1)*jl/2 +ik
              else
                CompositeIndex =  (ik-1)*ik/2 + jl
              end if

              end function CompositeIndex

     subroutine ReadInput(inputfile,norb,numfcore,numfvirt,numocc,numvirt,size_l,size_r,size_c,energy_start,energy_end,delta_en,volt_start,volt_end,delta_volt,inputcode,KT,Electrode_Type,Fermi_enl,Fermi_enr,CalcType,localden_fermi_l,localden_fermi_r,doubles,numatomic,functional,num_threads,use_b0,b0_type,write_ruqt_data,state_num)
     implicit none
     character(len=100) :: inputfile
     character(len=40) :: inputcode,filename,Electrode_Type,CalcType,functional,b0_type
     integer :: norb, size_c,size_r,size_l,numfcore,numfvirt,numocc,numvirt,numatomic,num_threads
     real(8) :: energy_start,energy_end,delta_en,volt_start,volt_end,delta_volt,KT
     logical :: libint,doubles,use_b0,write_ruqt_data
     real(8) :: Fermi_enl,Fermi_enr,localden_fermi_l,localden_fermi_r
     integer :: state_num

              filename = trim(inputfile)
              open(unit=1,file=filename,action="read")
              20 format(A) 
              num_threads=0
              read(1,20) CalcType
              read(1,20) Electrode_Type
              read(1,*) Fermi_enl
              read(1,*) Fermi_enr
              read(1,*) localden_fermi_l
              read(1,*) localden_fermi_r
              read(1,*) norb
              read(1,*) numatomic
              read(1,*) numfcore
              read(1,*) numfvirt
              read(1,*) numocc
              read(1,*) numvirt
              read(1,*) size_c
              read(1,*) size_l
              read(1,*) size_r
              read(1,*) energy_start
              read(1,*) energy_end
              read(1,*) delta_en
              read(1,*) volt_start
              read(1,*) volt_end
              read(1,*) delta_volt
              read(1,*) KT
              read(1,*) inputcode
              read(1,*) doubles
              read(1,*) functional
              read(1,*) use_b0
              read(1,*) b0_type
              read(1,*) write_ruqt_data
              read(1,*) num_threads
              read(1,*) state_num
                      
              close(1)
              if(num_threads.eq.0) then
                 num_threads=1
                end if
              end subroutine 

                      subroutine Flag_set(inputcode,functional,cisd_flag,rdm_flag,hf_flag,dft_flag,qchem,gamess,pyscf,maple,molcas)
                      implicit none

                      character(len=40) :: inputcode,functional
                      logical :: rdm_flag,cisd_flag,dft_flag,hf_flag
                      logical :: qchem,gamess,pyscf,maple,molcas

                       rdm_flag=.false.
                       cisd_flag=.false.
                       dft_flag=.false.
                       hf_flag=.false.
                       qchem=.false.
                       gamess=.false.
                       pyscf=.false.
                       maple=.false.
                       molcas=.false.

                      if(inputcode.eq."qchem") then
                         qchem=.true.
                       if(functional.eq."dft") then
                         dft_flag=.true. 
                        elseif(functional.eq."hf") then
                         hf_flag=.true.
                        else
                         write(*,*) "Your method and QC code choice do not work together. Exiting"
                         stop
                        end if

                       elseif(inputcode.eq."molcas") then
                         molcas=.true.
                       if(functional.eq."dft") then
                         dft_flag=.true.
                        elseif(functional.eq."hf") then
                         hf_flag=.true.
                        else
                         write(*,*) "Your method and QC code choice do not work together. Exiting"
                         stop
                        end if


                       elseif(inputcode.eq."gamess") then
                         gamess=.true.
                       if(functional.eq."rdm") then
                         rdm_flag=.true.
                elseif(functional.eq."hf") then
                 hf_flag=.true.
                elseif(functional.eq."cisd") then
                 cisd_flag=.true.
                else
                 write(*,*) "Your method and QC code choice do not work together. Exiting"
                 stop
                end if

               elseif(inputcode.eq."pyscf") then
                 pyscf=.true.
                 if(functional.eq."hf") then
                  hf_flag=.true.
                 elseif(functional.eq."dft") then
                  dft_flag=.true.
                 else
                  write(*,*) "Your method and QC code choice do not work together. Exiting"
                  stop
                 end if


               elseif(inputcode.eq."maple") then
                 maple=.true.
               if(functional.eq."rdm") then
                 rdm_flag=.true.
                elseif(functional.eq."hf") then
                 hf_flag=.true.
                elseif(functional.eq."cisd") then
                 cisd_flag=.true.
                elseif(functional.eq."dft") then
                 dft_flag=.true.
                else
                 write(*,*) "Your method and QC code choice do not work together. Exiting"
                 stop
                end if


               else
                write(*,*) "Your QC code choice is not supported. Exiting"
                stop
               end if
              end subroutine
