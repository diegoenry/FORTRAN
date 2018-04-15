! ChemFlow - Computational Chemistry is Great Again
! 
! IFP counter
!
! @Authors
! Priscila Gomes, Diego Gomes
!
! Change stuff where necessary
! IFP*256 will hold 5bit FPs for 51 residues. Do the math.
! 
program IFP_count
  character :: name*50,IFP*256
  integer*1,allocatable :: npolar(:)
  integer :: nres,nbits

  nres=51           ! number of residues
  nbits=5           ! number of fingerprint bits 
  nlines=10718123   ! number of lines in IFP file

  allocate(npolar(nres))

  open(1,file='ifp_polar.csv')  ! Open a clean IFP file.
  open(3,file='ifp_polar_count.dat')


  do i=1,nlines
    
    npolar=0

    read(1,*) name,IFP
    !write(*,*) IFP

    k=1 !counter for residue
    do j=1,nres*nbits,nbits

      ! Check if residue has IFPs
      if (index(IFP(j:j+nbits-1),"1") /= 0) then
         npolar(k)=1
      endif
      k=k+1
   enddo
    write(3,'(A,",",51i1)') trim(name),npolar

  enddo

  open(2,file='ifp_count.dat')
  write(2,*) npolar
end
