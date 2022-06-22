       Module InterfaceMod
        Interface

       subroutine Flag_set(inputcode,functional,cisd_flag,rdm_flag,hf_flag,dft_flag,qchem,gamess,pyscf,maple,molcas)
       implicit none
       character(len=40) :: inputcode,functional
       logical :: rdm_flag,cisd_flag,dft_flag,hf_flag
       logical :: qchem,gamess,pyscf,maple,molcas
       end subroutine

       subroutine Get_HF_QChem(inputfile,norb,H_two,Smat)
       implicit none

       character(len=100) :: inputfile
       integer :: norb
       real(8), allocatable, dimension(:,:) :: H_Two,Smat
       end subroutine

       subroutine Get_HF_Molcas(inputfile,norb,H_two,Smat,state_num)
       implicit none

       character(len=100) :: inputfile
       integer :: norb
       real(8), allocatable, dimension(:,:) :: H_Two,Smat
       integer :: state_num
       end subroutine

       subroutine Get_HF_GAMESS(inputfile,numatomic,H_two,Smat,norb)
       implicit none

       character(len=100) :: inputfile
       integer :: numatomic,norb
       real(8), allocatable, dimension(:,:) :: H_Two,Smat
       end subroutine

       subroutine Get_HF_PySCF(inputfile,numatomic,H_two,Smat,norb)
       implicit none

       character(len=100) :: inputfile
       integer :: numatomic,norb
       real(8), allocatable, dimension(:,:) :: H_Two,Smat
       end subroutine


       subroutine Get_HF_libint(inputfile,norb,numact,H_one,Smat,mo_coeff,OneInts,TwoIntsCompact)
       implicit none

       character(len=100) :: inputfile
       real(8),allocatable,dimension(:,:) :: H_one,H_two,Smat,OneInts,mo_coeff
       real(8),allocatable,dimension(:) :: TwoIntsCompact
       integer :: norb,numact
       end subroutine

       subroutine Calculate_Coupling_MoleculeWBL(Coupling_R,Coupling_L,localden_fermi)
       implicit none

       real(8) :: Coupling_R,Coupling_L,localden_fermi
       end subroutine

      subroutine Electrodes_MoleculeWBL(Sigma_l,Sigma_r,Smat,Coupling_R,Coupling_L,size_c,size_lc,size_lcr)
      implicit none

      integer :: size_lc,size_c,size_lcr
      real(8) :: Coupling_R,Coupling_L
      complex(8), allocatable, dimension(:,:) :: Sigma_L,Sigma_R
      real(8), allocatable, dimension(:,:) :: Smat_le,Smat_re,Smat
      end subroutine

      subroutine PartitionHS_MoleculeWBL(Smat,size_l,size_r,size_c,size_lc,size_lcr,Smat_el,Smat_er)
      implicit none

      real(8), allocatable, dimension(:,:) :: Smat,Smat_el,Smat_er
      integer :: size_l,size_c,size_r,size_lc,size_lcr
      end subroutine

      subroutine Electrodes_MetalWBL(Sigma_l,Sigma_r,Smat_re,Smat_le,H_Two_le,H_Two_re,localden_fermi_l,localden_fermi_r,size_c,size_lc,size_lcr,size_l,size_r,H_Two_le_trans,H_Two_re_trans,write_ruqt_data,inputfile)
      implicit none

      integer :: size_lc,size_c,size_lcr,size_l,size_r
      real(8) :: localden_fermi_l,localden_fermi_r
      complex(8), allocatable, dimension(:,:) :: Sigma_L,Sigma_R
      real(8), allocatable, dimension(:,:) :: Smat_le,Smat_re,H_Two_re,H_Two_le,H_Two_re_trans,H_Two_le_trans
      logical :: write_ruqt_data
      character(len=100) :: inputfile
      end subroutine

      subroutine PartitionHS_MetalWBL(Smat,H_Two,size_l,size_r,size_c,size_lc,size_lcr,Smat_el,Smat_er,Smat_cen,H_Two_le,H_Two_re,H_Two_cen,H_Two_le_trans,H_Two_re_trans,write_ruqt_data,inputfile)
      implicit none

      real(8), allocatable, dimension(:,:) :: Smat,Smat_el,Smat_er,Smat_cen,H_Two,H_Two_le,H_Two_re,H_Two_cen,H_Two_le_trans,H_Two_re_trans
      integer :: size_l,size_c,size_r,size_lc,size_lcr
      logical :: write_ruqt_data
      character(len=100) :: inputfile
      end subroutine

      subroutine ReadInput(inputfile,norb,numfcore,numfvirt,numocc,numvirt,size_l,size_r,size_c,energy_start,energy_end,delta_en,volt_start,volt_end,delta_volt,inputcode,KT,Electrode_Type,Fermi_enl,Fermi_enr,CalcType,localden_fermi_l,localden_fermi_r,doubles,numatomic,functional,num_threads,use_b0,b0_type,write_ruqt_data,state_num)
      implicit none

      character(len=100) :: inputfile
      character(len=40) :: inputcode,filename,Electrode_Type,CalcType,functional,b0_type
      integer :: norb,size_c,size_r,size_l,numfvirt,numfcore,numocc,numvirt,numatomic,num_threads,state_num
      real(8) :: energy_start,energy_end,delta_en,volt_start,volt_end,delta_volt,KT,Fermi_enl,Fermi_enr,localden_fermi_l,localden_fermi_r
      logical :: libint,doubles,rdm,use_b0,write_ruqt_data
      end subroutine

      Subroutine Build_G_SD_Invert(G_C,Sigma_l,Sigma_r,energy,size_l,size_c,size_lc,size_lcr,norb,inputfile,numocc,numvirt,iter,B1data,B2data,mo_ener,mo_coeff,mo_coeff2,doubles,currentflag,energy_values,ener_val,G_S,corr_ener,numatomic,B0_coeff,use_b0,gamess,maple,numfcore,numfvirt,b0_type)
      use TypeMod
      implicit none
      character(len=100) :: inputfile,Bfile
      character(len=40) :: b0_type
      type(B1) :: B1data
      type(B2) :: B2data
      type(energr) :: G_S
      integer :: size_l,size_lc,size_c,size_lcr,i,j,ioerror,norb,numocc,numvirt,iter,energy_values,ener_val,numatomic
      integer :: numfcore,numfvirt
      real(8) :: temp,corr_ener,hf_ener,energy
      real(8), allocatable, dimension(:,:) :: mo_coeff,mo_coeff2
      real(8), allocatable, dimension(:) :: mo_ener,B0_coeff
      complex(8), allocatable, dimension(:,:) :: G_C,Sigma_L,Sigma_R
      logical :: doubles,currentflag,use_b0,gamess,maple
      end subroutine

      end interface
      end module

      module InterfaceMod2
      interface

      subroutine Build_B0_CISD(B0_coeff,norb,numfcore,numfvirt,numocc,numvirt,B1data,B2data)
      use TypeMod
      implicit none
      integer :: size_l,size_lc,size_c,size_lcr,i,j,ioerror,norb,iter,a,b
      integer :: numocc,numvirt,p,q,z,y,x,r,s,t,energy_values,numatomic,numfvirt,numfcore
      type(B2) :: B2data
      type(B1) :: B1data
      real(8) :: B0_coeff
      end subroutine


      end interface 
      end module

      module FunctionMod
      interface
       function matmul_dgemm(leftmatrix,rightmatrix)
       Implicit NONE
       double precision, allocatable, dimension(:,:) :: leftmatrix,rightmatrix,matmul_dgemm
       end function

       function matmul_zgemm(leftmatrix,rightmatrix)
       Implicit NONE
       complex(8), allocatable, dimension(:,:) :: leftmatrix,rightmatrix,matmul_zgemm
       end function

       subroutine matmul_dgemm2(leftmatrix,rightmatrix,matsol)
       Implicit NONE
       double precision, allocatable, dimension(:,:) :: leftmatrix,rightmatrix,matsol
       end subroutine
         
       function inv(A) result(Ainv)
       implicit none
       complex(8), allocatable, dimension(:,:), intent(in) :: A
       complex(8), allocatable, dimension(:,:) :: Ainv
       end function

       function inv_real(A) result(Ainv)
       implicit none
       real(8), allocatable, dimension(:,:), intent(in) :: A
       real(8), allocatable, dimension(:,:) :: Ainv
       end function

       function adjoint(A,norb)
       implicit none

       integer :: norb
       complex(8),allocatable,dimension(:,:) :: A,adjoint
       end function

      function fermi_function(energy,fermi_energy,KT)
      implicit none
      real(8) :: energydiff,KT,energy,fermi_energy,fermi_function
      end function

       function FirstIndex(i,k)
       Implicit None
       integer :: i,k
       integer(8) :: FirstIndex
       end function

       function CompositeIndex(ik,jl)
       Implicit NONE

       integer(8) :: ik, jl
       integer(8) :: CompositeIndex
       end function        

      end interface
      end module
