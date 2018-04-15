! ChemFlow - Computational Chemistry is Great Again
!
! @brief
! Reads an multimolecule IFP and reports the number of IFP matches for each
! interaction, by residue.
!
! So far only polar, and aromatic are implemented.
! 
! IFP by residue report - for IChem
!
! @Authors
! Priscila S.F.C. Gomes, Diego Enry B. Gomes*
! mailto: dgomes@pq.cnpbr
!
! Change stuff where necessary
! IFP*3000 will hold 9bit FPs for 51 residues (and a bit more).
! 
! IChem IFP positions
! 1st : hydrophobic
! 2nd : aromatic (face-to-face)
! 3rd : aromatic (edge-to-face)
! 4th : h-bond (protein donor)
! 5th : h-bond (ligand donor)
! 6th : ionic (protein charged +)
! 7th : ionic (ligand charged +)
! 8th : pi-cation
! 9th : metal coordination)
!
! Sun Apr 15 21:08:38 CEST 2018
program IFP_aromatic
  character :: name*50
  character :: IFP*3000
  character :: filename*256

  integer,allocatable   :: nIFP_hydrophobic(:)
  integer,allocatable   :: nIFP_polar(:) 
  integer,allocatable   :: nIFP_naromatic_f2f(:)
  integer,allocatable   :: nIFP_naromatic_e2f(:)

  integer*1,allocatable :: nhydrophobic(:)
  integer*1,allocatable :: npolar(:)
  integer*1,allocatable :: naromatic_f2f(:)
  integer*1,allocatable :: naromatic_e2f(:)
  integer :: nres,nbits

  write(*,*) "IFP file name: "
  read(*,*)  filename

  write(*,*) "Number of residues ? "
  read(*,*)  nres

  write(*,*) "Number molecules with IFP ? ( wc -l filename to discover) "
  read(*,*)  nlines

  write(*,*) "Number of bits in IFP (reduced = 5, normal = 7 or extended = 9)"
  read(*,*)  nbits

  open(1,file=trim(filename))  ! Open a clean IFP file.

  ! Manual configuration 
!  nres=51        ! number of residues
!  nbits=9        ! number of fingerprint bits 
!  nlines=232     ! number of lines in IFP file
!  open(1,file='sample.ifp')  ! Open a clean IFP file.
! End of manual configuration. 
  

allocate(nIFP_hydrophobic(nres))
allocate(nIFP_polar(nres))
allocate(nIFP_naromatic_f2f(nres))
allocate(nIFP_naromatic_e2f(nres))

allocate(nhydrophobic(nres))
allocate(npolar(nres))
allocate(naromatic_f2f(nres))
allocate(naromatic_e2f(nres))

! Zero out vectors 
nIFP_hydrophobic=0
nIFP_polar=0
nIFP_naromatic_f2f=0
nIFP_naromatic_e2f=0

do i=1,nlines
   
  ! Zero out temporary
  nhydrophobic=0
  npolar=0 
  naromatic_f2f=0
  naromatic_e2f=0

  read(1,*) name,IFP

  k=1 !counter for residue
  do j=1,nres*nbits,nbits

  ! Check if residue has any hydrophobic IFPs
    if (index(IFP(j:j),"1") /= 0) then
       nhydrophobic(k)=1
    endif

  ! Check if residue has any polar IFPs (pi-cation or metal)
    if (index(IFP(j+3:j+8),"1") /= 0) then
       npolar(k)=1
    endif

  ! Check if residue has aromatic face-to-face IFPs
    if (index(IFP(j+1:j+1),"1") /= 0) then
         naromatic_f2f(k)=1
      endif

  ! Check if residue has aromatic edge-to-face IFPs
    if (index(IFP(j+2:j+2),"1") /= 0) then
       naromatic_e2f(k)=1
    endif

  ! Counter for number of molecules with polar IFP to the current residue
    nIFP_hydrophobic(k)=nIFP_hydrophobic(k)+nhydrophobic(k)
    nIFP_polar(k)=nIFP_polar(k)+npolar(k)
    nIFP_naromatic_f2f(k)=nIFP_naromatic_f2f(k)+naromatic_f2f(k)
    nIFP_naromatic_e2f(k)=nIFP_naromatic_e2f(k)+naromatic_e2f(k)

    k=k+1
      
enddo
!     write(2,'(A,",",i2,",",i2)') trim(name),sum(naromatic_f2f),sum(naromatic_e2f)
!    write(*,'(A,51(",",i1))') trim(name), npolar
  enddo

  write(*,'(A,51(",",i5))') "          Hydrophobic",nIFP_hydrophobic
  write(*,'(A,51(",",i5))') "                Polar",nIFP_polar
  write(*,'(A,51(",",i5))') "Aromatic Face-to-Face",nIFP_naromatic_f2f
  write(*,'(A,51(",",i5))') "Aromatic Edge-to-Face",nIFP_naromatic_e2f

end
