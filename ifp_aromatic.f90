! ChemFlow - Computational Chemistry is Great Again
! 
! IFP report for IChem[1]
!
! @Authors
! Priscila S.F.C. Gomes, Diego Enry B. Gomes*
! mailto: dgomes@pq.cnpbr
!
! Change stuff where necessary
! IFP*256 will hold 5bit FPs for 51 residues. Do the math.
! 
! Sun Apr 15 12:19:22 CEST 2018
program IFP_aromatic
  character :: name*50,IFP*3000
  integer*1,allocatable :: naromatic_f2f(:)
  integer*1,allocatable :: naromatic_e2f(:)
  integer :: nres,nbits

  nres=242          ! number of residues
  nbits=9           ! number of fingerprint bits 
  nlines=26727      ! number of lines in IFP file

  allocate(naromatic_f2f(nres))
  allocate(naromatic_e2f(nres))


  open(1,file='TUDO.IFP')  ! Open a clean IFP file.
  open(2,file='tudo.csv')

  write(2,'(A)') "Name,Face-to-Face,Edge-to-Face" 


  do i=1,nlines
    
    naromatic_f2f=0
    naromatic_e2f=0

    read(1,*) name,IFP

    k=1 !counter for residue
    do j=1,nres*nbits,nbits

      ! Check if residue has aromatic face-to-face IFPs
      if (index(IFP(j+1:j+1),"1") /= 0) then
         naromatic_f2f(k)=1
      endif

      ! Check if residue has aromatic edge-to-face IFPs
      if (index(IFP(j+2:j+2),"1") /= 0) then
         naromatic_e2f(k)=1
      endif

      k=k+1

   enddo
     write(2,'(A,",",i2,",",i2)') trim(name),sum(naromatic_f2f),sum(naromatic_e2f)

  enddo

end
