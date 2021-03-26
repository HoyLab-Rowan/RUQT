       Module TypeMod
       Implicit None

      !B2 type
      type :: inner
           real(8), allocatable, dimension(:,:) :: n
      end type inner

      type :: outer
           type(inner), allocatable, dimension(:,:) :: m
      end type outer

      type :: B2
           type(outer) :: aa
           type(outer) :: ab
           type(outer) :: bb
           type(outer) :: nospin
      end type B2

      !B1 type
      type :: virtorb
           real(8), allocatable, dimension(:,:) :: o
      end type virtorb

      type :: B1
           type(virtorb) :: a
           type(virtorb) :: b
      end type B1

      type :: ingreen
        real(8), allocatable,dimension(:,:) :: gf
      end type ingreen

      type :: energr
        type(ingreen),allocatable,dimension(:) :: en
      end type energr

      integer, allocatable, dimension(:) :: Virt_Index,Occ_Index

      end Module
