!program hme_wrapper
!implicit none
!integer(kind=4)::natom,ngrid,qtot
!real(kind=8),allocatable::gridcrd(:,:),atomcrd(:,:)
!real(kind=8),allocatable::esp(:),mpinit(:),dpinit(:,:),qpinit(:,:)
!integer(kind=4)::level(3)
!character(len=30)::espfile
!character(len=10)::fitmethod
!integer(kind=4)::istrnt
!real(kind=8)::strength
!integer(kind=4)::nlgn
!integer(kind=4),allocatable::lgnatomidx(:,:)
!real(kind=8),allocatable::lgncons(:)
!integer(kind=4)::nlgnatom_in_this_cnstrnt
!integer(kind=4)::ilgn
!integer(kind=4)::nsymmconstraint,isymmconstraint
!integer(kind=4),allocatable::isymmconstraintpair(:,:)
!integer(kind=4),allocatable::atmidx(:)
!integer(kind=4),allocatable::coef(:)
!real(kind=8)::ssres,crosscorr
!integer(kind=4)::i,j,k,m,n
!espfile='aceticacid.esp'
!open(11,file=espfile)
!read(11,*)natom,ngrid
!if(allocated(gridcrd))deallocate(gridcrd)
!allocate(gridcrd(3,ngrid))
!if(allocated(atomcrd))deallocate(atomcrd)
!allocate(atomcrd(3,natom))
!if(allocated(esp))deallocate(esp)
!allocate(esp(ngrid))
!if(allocated(mpinit))deallocate(mpinit)
!allocate(mpinit(natom))
!if(allocated(dpinit))deallocate(dpinit)
!allocate(dpinit(3,natom))
!if(allocated(qpinit))deallocate(qpinit)
!allocate(qpinit(5,natom))
!if(allocated(atmidx))deallocate(atmidx)
!allocate(atmidx(natom))
!if(allocated(coef))deallocate(coef)
!allocate(coef(natom))
!read(11,'(16X,3E16.7)')((atomcrd(i,j),i=1,3),j=1,natom)
!do j = 1, natom
!!  write(*,'(A,I7,16X,3ES16.7)')'Atom ', j, atomcrd(:,j)
!end do
!read(11,'(4E16.7)')((esp(j),(gridcrd(i,j),i=1,3)),j=1,ngrid)
!do j = 1, ngrid
!!   write(*,'(A,I7,4ES16.7)')'Grid ', j, esp(j),gridcrd(:,j)
!end do
!close(11)
!qtot=0
!mpinit(1)=  -0.602012
!mpinit(2)=   0.224271
!mpinit(3)=   0.224271
!mpinit(4)=   0.224271
!mpinit(5)=   0.762176
!mpinit(6)=  -0.701809
!mpinit(7)=   0.459737
!mpinit(8)=  -0.533535
!mpinit=0.d0
!dpinit=0.d0
!qpinit=0.d0
!fitmethod='resp'
!istrnt=1
!strength=0.0005d0
!level(1)=1
!level(2)=1
!level(3)=1
!read*,nlgn
!allocate(lgnatomidx(natom,nlgn))
!allocate(lgncons(nlgn))
!do ilgn = 1, nlgn
! ! parse the input for Lagrange constraints
!  read(*,*)nlgnatom_in_this_cnstrnt,lgncons(ilgn)
!  lgnatomidx(:,ilgn)=0
!  coef=0
!  read(*,'(8I5)')((atmidx(j),Coef(j)),j=1,nlgnatom_in_this_cnstrnt)
!  do j = 1, nlgnatom_in_this_cnstrnt
!    lgnatomidx(atmidx(j),ilgn)=Coef(j)
!  end do
!end do
!read*,nsymmconstraint
!allocate(isymmconstraintpair(2,nsymmconstraint))
!do isymmconstraint = 1, nsymmconstraint
!  read(*,'(2I5)')isymmconstraintpair(:,isymmconstraint)
!end do
!call hme(natom,ngrid,gridcrd,atomcrd,level,esp,qtot,mpinit,dpinit,qpinit,fitmethod,istrnt,strength,nlgn,lgnatomidx,lgncons,nsymmconstraint,isymmconstraintpair,40,ssres,crosscorr)
!deallocate(gridcrd,atomcrd,esp,mpinit,dpinit,qpinit,atmidx,coef,lgnatomidx,lgncons,isymmconstraintpair)
!
!end program hme_wrapper



subroutine hme(natom,ngrid,gridcrd,atomcrd,level,espinit,qtot,mpinit,dpinit,qpinit,fitmethod0,  &
istrnt,strength,nlgn0,lgnatomidx0,lgncons0,nsymmconstraint,isymmconstraintpair,outid,ssres,crosscorr)
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              Restrained Electrostatic Potential Fitting of Atomic Charges               !!
!!                    Written by Ye Mei, East China Normal University                      !!
!!            During his stay at Dr. Bernie Brooks' group as Special Volunteer             !!
!!                             Ver. 1.0, Feb. 12, 2014                                     !!
!!                               All rights Reserved                                       !!
!!  References:                                                                            !!
!!  1. C. I. Bayly, P. Cieplak, W. Cornell, P. A. Kollman,                                 !!
!!     J. Phys. Chem., 1993, 97 (40), pp 10269–10280                                       !!
!!  2. J. Zeng, L. L., Duan, J. Z. H., Zhang, Y. Mei,                                      !!
!!     J. Comput. Chem., 34, 847-853 (2013)                                                !!
!!  Input variables                                                                        !!
!!    natom: Number of atoms (multipole sites)                                             !!
!!    ngrid: Number of ESP grid points                                                     !!
!!    gridcrd(3,ngrid): Coordinates of ESP grid points                                     !!
!!    atomcrd(3,natoms): Coordinates of atoms (multipole sites)                            !!
!!    level(3): if level(i)=1, multipole L=i-1 will be fitted. if 0, No.  =2? No Kidding!  !!
!!    espinit(ngrid): Standard (QM) ESP on grids                                           !!
!!    qtot : total charge                                                                  !!
!!    mpinit(natom): Initial guess of monopole on sites, the fitted monopole on return     !!
!!    dpinit(natom): Initial guess of dipole on sites, the fitted dipole on return         !!
!!    qpinit(natom): Initial guess of quadrupole on sites, the fitted quadrupole on return !!
!!    fitmethod: Fitting method (RESP, DRESP, to be added)                                 !!
!!    istrnt: restraint function (0: Harmonic; 1: hyperbolic                               !!
!!    strength: strength of restraint                                                      !!
!!    nlgn0: number of Lagrange constraints (excluding the constraint of total charge)     !!
!!    lgnatomidx0: coefficients for the atoms in each Lagrange constraints                 !!
!!    lgncons0: constraint for each Lagrange constraints                                   !!
!!    nsymmconstraint: number of symmetry constraints                                      !!
!!    isymmconstraintpair: site pair in symmetry                                           !!
!!    outid: file pipe for output                                                          !!
!!    ssres: residual sum of squares                                                       !!
!!    crosscorr: correlation between fitted and target ESP                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

interface
  subroutine getlgn(natom,nlgn0,lgnatomidx0,lgncons0,nlgn,lgnatomidx,lgncons,qresidual,qbase,outid)
  integer(kind=4), intent(in) :: natom
  integer(kind=4), intent(in) :: nlgn0
  integer(kind=4), intent(in) :: lgnatomidx0(natom,nlgn0)
  real(kind=8), intent(in) :: lgncons0(nlgn0)
  integer(kind=4), intent(out) :: nlgn
  integer(kind=4), allocatable, intent(out) :: lgnatomidx(:,:)
  real(kind=8), allocatable, intent(out) :: lgncons(:)
  real(kind=8), intent(in) :: qresidual
  real(kind=8), intent(in) :: qbase(natom)
  integer(kind=4), intent(in) :: outid
  end subroutine getlgn

  subroutine gendismat(natom,ngrid,atomcrd,gridcrd,r_atm_grid,r_atm_grid_vec,ibigmem,outid)
  integer(kind=4), intent(in) :: natom, ngrid
  real(kind=8), intent(in) :: atomcrd(3,natom), gridcrd(3,ngrid)
  real(kind=8), allocatable, intent(out) :: r_atm_grid(:,:), r_atm_grid_vec(:,:,:)
  integer(kind=4), intent(out) :: ibigmem
  integer(kind=4), intent(in) :: outid
  end subroutine gendismat
end interface

integer(kind=4), intent(in) :: natom, ngrid
real(kind=8), intent(in) :: gridcrd(3,ngrid), atomcrd(3,natom), espinit(ngrid)
integer(kind=4), intent(in) :: level(3)
integer(kind=4), intent(in) :: qtot
real(kind=8), intent(in out) :: mpinit(natom), dpinit(3,natom), qpinit(5,natom)
character(len=*), intent(in) :: fitmethod0
integer(kind=4), intent(in) :: istrnt
real(kind=8), intent(in) :: strength
integer(kind=4), intent(in) :: nlgn0
integer(kind=4), intent(in) :: lgnatomidx0(natom,nlgn0)
real(kind=8), intent(in) :: lgncons0(nlgn0)
integer(kind=4), intent(in) :: nsymmconstraint
integer(kind=4), intent(in) :: isymmconstraintpair(2,nsymmconstraint)
integer(kind=4), intent(in) :: outid
real(kind=8), intent(out) :: ssres, crosscorr


character(len=10)::fitmethod
real(kind=8) :: espfitted(ngrid), espresidual(ngrid)
integer(kind=4) :: nlgn
integer(kind=4), allocatable :: lgnatomidx(:,:)
real(kind=8), allocatable :: lgncons(:)
integer(kind=4) :: isymmconstraintpairlocal(2,nsymmconstraint)
real(kind=8), allocatable :: q(:), qcopy(:)
real(kind=8), allocatable :: awork(:,:), bwork(:)
real(kind=8), allocatable :: acopy(:,:), bcopy(:)
real(kind=8) :: q0(natom), qbase(natom), qresidual
real(kind=8) :: monopoleout(natom), dipoleout(3,natom), quadrupoleout(5,natom)
integer(kind=4) :: neffecsize
real(kind=8), allocatable :: qwt(:)
real(kind=8), allocatable :: r_atm_grid(:,:), r_atm_grid_vec(:,:,:)
real(kind=8) :: sstot, ssreg
integer(kind=4) :: ibigmem
logical :: isconverged
integer(kind=4) :: istep
integer(kind=4) :: icycle
real(kind=8), parameter :: criterion=1.d-4
integer(kind=4),allocatable :: iskilledsite(:)          ! some charge sites are not freely variable due to the existence of constraints
integer(kind=4) :: deallocatestatus
integer(kind=4) :: AllocateStatus
integer(kind=4) :: i, j, k, m, n

fitmethod=trim(adjustl(fitmethod0))
do i = 1, len(trim(fitmethod))
  if(ichar(fitmethod(i:i))<=ichar('z').and.ichar(fitmethod(i:i))>=ichar('a'))then
    fitmethod(i:i)=char(ichar(fitmethod(i:i))-ichar('a')+ichar('A'))
  end if
end do
if(fitmethod/='RESP'.and.fitmethod/='DRESP')then
  write(outid,'(A)')' Error: Only RESP and DRESP methods are allowed. '
  write(outid,'(A)')' Error: Please check your input for "fitmethod" and run me again.'
  stop
else 
  write(outid,'(A)')' Running '//trim(fitmethod)//' multipole fitting method'
end if

write(outid,'(A)')' Fit multipole moments using '//trim(fitmethod)//' fitting method'
write(outid,'(A,I7)')' Number of atoms ', natom
write(outid,'(A,I7)')' Number of grids ', ngrid
write(outid,*)

if(fitmethod/='RESP'.and.fitmethod/='DRESP')then
  write(outid,'(A)')' Initial multipole moments'
  write(outid,'(A)')' Monopoles'
  write(outid,'(8F10.6)')(mpinit(i),i=1,natom)
  write(outid,*)
  write(outid,'(A)')' Dipoles'
  do i = 1, natom
    write(outid,'(3F10.6)')dpinit(:,i)
  end do
  write(outid,*)
  write(outid,'(A)')' Quadrupoles'
  do i = 1, natom
    write(outid,'(5F10.6)')qpinit(:,i)
  end do
  write(outid,*)
  ! give a score to the initial multipole
  call multipoleESP(ngrid,gridcrd,natom,atomcrd,mpinit,dpinit,qpinit,espfitted)
  call qualitymetric(ngrid,espinit,espfitted,sstot,ssreg,ssres,crosscorr)
  write(outid,'(A)')'    Statistics of the initial multipoles:'
  write(outid,'(A,F10.3)')'  The initial sum of squares (ssvpot)      ', sstot
  write(outid,'(A,F10.3)')'  The residual sum of squares (chipot)     ', ssres
  write(outid,'(A,F10.5)')'  The std err of estimate (sqrt(chipot/N)) ', sqrt(ssres/ngrid)
  write(outid,'(A,F10.5)')'  ESP relative RMS (SQRT(chipot/ssvpot))   ', sqrt(ssres/sstot)
  write(outid,'(A,F10.5)')'  Correlation coefficient                  ', 1.d0-ssres/sstot
  write(outid,'(A,F10.5)')'  Cross correlation                        ', crosscorr
  write(outid,*)
end if

! generate distance matrix between atoms and grids if memory allows
call gendismat(natom,ngrid,atomcrd,gridcrd,r_atm_grid,r_atm_grid_vec,ibigmem,outid)

! initialize 
if(allocated(qwt))then
  deallocate(qwt,stat=deallocatestatus)
  if(deallocatestatus/=0)then
    write(outid,'(A)')'deallocating qwt failed. deallocatestatus=',deallocatestatus
    stop
  end if
end if
allocate(qwt(natom),stat=AllocateStatus)
if(AllocateStatus/=0)then
  write(outid,'(A)')'Failed when allocating qwt'
  stop
end if


call initialize(natom,ngrid,gridcrd,atomcrd,espinit,espresidual,mpinit,qbase,qwt,r_atm_grid,ibigmem,fitmethod,outid)

write(outid,'(A,I8)')' Total charge of this molecule:', qtot
write(outid,'(A,F16.5)')' Sum of initial charge:', sum(mpinit)
qresidual=qtot-sum(qbase)
if(fitmethod=='DRESP')then
  write(outid,'(A,F14.5)')' Residueal charge to fit:', qresidual
end if
write(outid,*)

if(allocated(iskilledsite))then
  deallocate(iskilledsite,stat=deallocatestatus)
  if(deallocatestatus/=0)then
    write(outid,'(A)')'deallocating iskilledsite failed. deallocatestatus=',deallocatestatus
    stop
  end if
end if
allocate(iskilledsite(natom),stat=AllocateStatus)
if(AllocateStatus/=0)then
  write(outid,'(A)')'Failed when allocating iskilledsite'
  stop
end if
iskilledsite=0

espfitted = 0.d0
istep = 0
if(level(1)<=0)then
  if(abs(qresidual)>1.D-4)then
    write(outid,'(A)')'Error: Monopole fitting is disabled. However, I noticed that the sum of the initial charge does not equal to qtot.'
    write(outid,'(A)')'Error: This is abnormal.'
    stop
  end if
else
  istep = istep + 1
  write(outid,'(A,I2,A)')' Step ',istep,': To fit monopoles'
  ! get info for Lagrange multiplier
  call getlgn(natom,nlgn0,lgnatomidx0,lgncons0,nlgn,lgnatomidx,lgncons,qresidual,qbase,outid)
  
  isymmconstraintpairlocal=isymmconstraintpair
  write(outid,'(A,I5)')' Total number of symmetry constraints:',nsymmconstraint
  do i = 1, nsymmconstraint
    write(outid,'(2I3)')isymmconstraintpairlocal(:,i)
  end do
  write(outid,*)

  if(allocated(awork))then
    deallocate(awork,stat=deallocatestatus)
    if(deallocatestatus/=0)then
      write(outid,'(A)')'deallocating awork failed. deallocatestatus=',deallocatestatus
      stop
    end if
  end if
  allocate(awork(natom+nlgn,natom+nlgn),stat=AllocateStatus)
  if(AllocateStatus/=0)then
    write(outid,'(A)')'Failed when allocating awork'
    stop
  end if
  
  if(allocated(bwork))then
    deallocate(bwork,stat=deallocatestatus)
    if(deallocatestatus/=0)then
      write(outid,'(A)')'deallocating bwork failed. deallocatestatus=',deallocatestatus
      stop
    end if
  end if
  allocate(bwork(natom+nlgn),stat=AllocateStatus)
  if(AllocateStatus/=0)then
    write(outid,'(A)')'Failed when allocating bwork'
    stop
  end if

  if(allocated(q))then
    deallocate(q,stat=deallocatestatus)
    if(deallocatestatus/=0)then
      write(outid,'(A)')'deallocating q failed. deallocatestatus=',deallocatestatus
      stop
    end if
  end if
  allocate(q(natom),stat=AllocateStatus)
  if(AllocateStatus/=0)then
    write(outid,'(A)')'Failed when allocating q'
    stop
  end if

  neffecsize=natom+nlgn-nsymmconstraint
  if(allocated(acopy))then
    deallocate(acopy,stat=deallocatestatus)
    if(deallocatestatus/=0)then
      write(outid,'(A)')'deallocating acopy failed. deallocatestatus=',deallocatestatus
      stop
    end if
  end if
  allocate(acopy(neffecsize,neffecsize),stat=AllocateStatus)
  if(AllocateStatus/=0)then
    write(outid,'(A)')'Failed when allocating acopy'
    stop
  end if

  if(allocated(bcopy))then
    deallocate(bcopy,stat=deallocatestatus)
    if(deallocatestatus/=0)then
      write(outid,'(A)')'deallocating bcopy failed. deallocatestatus=',deallocatestatus
      stop
    end if
  end if
  allocate(bcopy(neffecsize),stat=AllocateStatus)
  if(AllocateStatus/=0)then
    write(outid,'(A)')'Failed when allocating bcopy'
    stop
  end if

  if(allocated(qcopy))then
    deallocate(qcopy,stat=deallocatestatus)
    if(deallocatestatus/=0)then
      write(outid,'(A)')'deallocating bcopy failed. deallocatestatus=',deallocatestatus
      stop
    end if
  end if
  allocate(qcopy(neffecsize),stat=AllocateStatus)
  if(AllocateStatus/=0)then
    write(outid,'(A)')'Failed when allocating qcopy'
    stop
  end if

  monopoleout=0.d0
  dipoleout=0.d0
  quadrupoleout=0.d0
  q0=mpinit-qbase
  isconverged=.false.
  icycle=0
  do while(.not.isconverged)
    icycle = icycle + 1
    write(outid,'(A,I3)')' Cycle:', icycle
    ! generate least-square matrix 
    call bldmatmonopole(natom,nlgn,atomcrd,ngrid,gridcrd,espresidual,r_atm_grid,ibigmem,awork,bwork)
!    write(outid,'(<natom>F10.5,2X,F10.4)')((awork(i,j),j=1,natom),bwork(i),i=1,natom)
!    write(outid,*)
    ! modify working matrices with charge restraints
    call getrestraints(natom,nlgn,awork,bwork,istrnt,strength,q0,qwt)
!    write(outid,'(<natom>F10.5,2X,F10.4)')((awork(i,j),j=1,natom),bwork(i),i=1,natom)
!    write(outid,*)
    ! add Lagrange condition to least-square matrix
    call appendlgn2mat(natom,nlgn,awork,bwork,lgnatomidx,lgncons)
!  write(outid,'(<natom+nlgn>F10.5,2X,F10.4)')((awork(i,j),j=1,natom+nlgn),bwork(i),i=1,natom+nlgn)
!  write(outid,*)

    if(nsymmconstraint>0)then
      call sortsymmconstraint(nsymmconstraint,isymmconstraintpairlocal)
      call applysymmetry(natom,nlgn,awork,bwork,nsymmconstraint,isymmconstraintpair,iskilledsite)
   end if

!  write(outid,'(<natom+nlgn>F10.5,2X,F10.4)')((awork(i,j),j=1,natom+nlgn),bwork(i),i=1,natom+nlgn)
!  write(outid,*)

  ! backups for working matrices. Essential for iterations
    m = 0
    do i = 1, natom+nlgn
      if(i<=natom)then
        if(iskilledsite(i)/=0)cycle
      end if
      m=m+1
      n=0
      do j = 1, natom+nlgn
        if(j<=natom)then
          if(iskilledsite(j)/=0)cycle
        end if
        n=n+1
        acopy(m,n)=awork(i,j)
      end do
      bcopy(m)=bwork(i)
    end do
!    write(outid,'(<neffecsize>F10.5,2X,F10.4)')((acopy(i,j),j=1,neffecsize),bcopy(i),i=1,neffecsize)
!    write(outid,*)
    do i = 1, neffecsize
      if(abs(acopy(i,i))<1.D-10)acopy(i,i)=1.D-10
    end do   
    
!    if(allocated(awork))deallocate(awork)
!    allocate(awork(neffecsize,neffecsize))
!    if(allocated(bwork))deallocate(bwork)
!    allocate(bwork(neffecsize))
!    awork=acopy
!    bwork=bcopy
    ! SVD fit of charges
    call LESVD(neffecsize,neffecsize,acopy,bcopy,qcopy,outid)
!    call LELUD(neffecsize,neffecsize,acopy,bcopy,qcopy,outid)
!    print*,'killed charge:'
!    do i = 1, natom
!      if(iskilledsite(i)>0)print*,i
!    end do

    do i = 1, min(natom,neffecsize)
      q(i) = qcopy(i)
    end do

    call rebuildchargearray(natom,q,nsymmconstraint,isymmconstraintpairlocal,iskilledsite)

    forall(i=1:natom)
      monopoleout(i)=qbase(i)+q(i)
    end forall
    write(outid,'(A)')' Fitted monopoles'
    write(outid,'(8F10.6)')(q(i),i=1,natom)
    write(outid,*)
    if(fitmethod=='DRESP')then
      write(outid,'(A)')' Fitted + initial monopoles'
      write(outid,'(8F10.6)')(monopoleout(i),i=1,natom)
      write(outid,*)
    end if
  
    ! check if converged
    if(istrnt==1)then
      call checkconvergence(natom,q,q0,criterion,isconverged)
      forall(i=1:natom)
        q0(i)=q(i)
      end forall
    else
      isconverged=.true.
    endif
  end do
! 
! call multipoleESP(ngrid,gridcrd,natom,atomcrd,monopoleout,dipoleout,quadrupoleout,espfitted)
! call qualitymetric(ngrid,espinit,espfitted,sstot,ssreg,ssres,crosscorr)
! write(outid,'(A)')'  Statistics of the fitted multipoles (up to L=0):'
! write(outid,'(A,F10.3)')'  The initial sum of squares (ssvpot)      ', sstot
! write(outid,'(A,F10.3)')'  The residual sum of squares (chipot)     ', ssres
! write(outid,'(A,F10.5)')'  The std err of estimate (sqrt(chipot/N)) ', sqrt(ssres/ngrid)
! write(outid,'(A,F10.5)')'  ESP relative RMS (SQRT(chipot/ssvpot))   ', sqrt(ssres/sstot)
! write(outid,'(A,F10.5)')'  Correlation coefficient (1-chipot/ssvpot)', 1.d0-ssres/sstot
! write(outid,'(A,F10.5)')'  Cross correlation                        ', crosscorr
! write(outid,*)
end if




if(level(2)>0)then
  istep = istep + 1
  write(outid,'(A,I2,A)')' Step ',istep,': To fit dipoles'
  if(allocated(awork))deallocate(awork)
  allocate(awork(3*natom,3*natom))
  if(allocated(bwork))deallocate(bwork)
  allocate(bwork(3*natom))
  
  espresidual=espinit-espfitted
  call bldmatdipole(natom,atomcrd,ngrid,gridcrd,espresidual,r_atm_grid,r_atm_grid_vec,ibigmem,awork,bwork)
  call LESVD(3*natom,3*natom,awork,bwork,dipoleout,outid)
  write(outid,'(A)')' Fitted dipoles'
  do i = 1, natom
    write(outid,'(3F10.6)')dipoleout(:,i)
  end do
  write(outid,*)
  
  call multipoleESP(ngrid,gridcrd,natom,atomcrd,monopoleout,dipoleout,quadrupoleout,espfitted)
  call qualitymetric(ngrid,espinit,espfitted,sstot,ssreg,ssres,crosscorr)
  write(outid,'(A)')'  Statistics of the fitted multipoles (up to L=1):'
  write(outid,'(A,F10.3)')'  The initial sum of squares (ssvpot)      ', sstot
  write(outid,'(A,F10.3)')'  The residual sum of squares (chipot)     ', ssres
  write(outid,'(A,F10.5)')'  The std err of estimate (sqrt(chipot/N)) ', sqrt(ssres/ngrid)
  write(outid,'(A,F10.5)')'  ESP relative RMS (SQRT(chipot/ssvpot))   ', sqrt(ssres/sstot)
  write(outid,'(A,F10.5)')'  Correlation coefficient (1-chipot/ssvpot)', 1.d0-ssres/sstot
  write(outid,'(A,F10.5)')'  Cross correlation                        ', crosscorr
  write(outid,*)
end if
  
if(level(3)>0)then
  istep = istep + 1
  write(outid,'(A,I2,A)')' Step ',istep,': To fit quadrupoles'
  if(allocated(awork))deallocate(awork)
  allocate(awork(5*natom,5*natom))
  if(allocated(bwork))deallocate(bwork)
  allocate(bwork(5*natom))
  
  espresidual=espinit-espfitted
  call bldmatquadrupole(natom,atomcrd,ngrid,gridcrd,espresidual,r_atm_grid,r_atm_grid_vec,ibigmem,awork,bwork)
  call LESVD(5*natom,5*natom,awork,bwork,quadrupoleout,outid)
  write(outid,'(A)')' Fitted quadrupoles'
  do i = 1, natom
    write(outid,'(5F10.6)')quadrupoleout(:,i)
  end do
  write(outid,*)
  
  call multipoleESP(ngrid,gridcrd,natom,atomcrd,monopoleout,dipoleout,quadrupoleout,espfitted)
  call qualitymetric(ngrid,espinit,espfitted,sstot,ssreg,ssres,crosscorr)
  write(outid,'(A)')'  Statistics of the fitted multipoles (up to L=2):'
  write(outid,'(A,F10.3)')'  The initial sum of squares (ssvpot)      ', sstot
  write(outid,'(A,F10.3)')'  The residual sum of squares (chipot)     ', ssres
  write(outid,'(A,F10.5)')'  The std err of estimate (sqrt(chipot/N)) ', sqrt(ssres/ngrid)
  write(outid,'(A,F10.5)')'  ESP relative RMS (SQRT(chipot/ssvpot))   ', sqrt(ssres/sstot)
  write(outid,'(A,F10.5)')'  Correlation coefficient (1-chipot/ssvpot)', 1.d0-ssres/sstot
  write(outid,'(A,F10.5)')'  Cross correlation                        ', crosscorr
  write(outid,*)
end if

write(outid,'(A)')' Summary of multipoles'
write(outid,'(A)')' Fitted monopoles'
write(outid,'(8F10.6)')(monopoleout(i),i=1,natom)
write(outid,*)
write(outid,'(A)')' Fitted dipoles'
do i = 1, natom
  write(outid,'(3F10.6)')dipoleout(:,i)
end do
write(outid,*)
write(outid,'(A)')' Fitted quadrupoles'
do i = 1, natom
  write(outid,'(5F10.6)')quadrupoleout(:,i)
end do
write(outid,*)

if(level(1)==1)then
  mpinit=monopoleout
else if (level(2)==1) then
  dpinit=dipoleout
else if (level(3)==1)then
  qpinit=quadrupoleout
end if

if(allocated(awork))deallocate(awork)
if(allocated(bwork))deallocate(bwork)
if(allocated(acopy))deallocate(acopy,STAT=deallocatestatus)
if(allocated(bcopy))deallocate(bcopy)
if(allocated(bcopy))deallocate(qcopy)
if(allocated(q))deallocate(q)
if(allocated(qwt))deallocate(qwt)
if(allocated(lgnatomidx))deallocate(lgnatomidx)
if(allocated(lgncons))deallocate(lgncons)
if(allocated(r_atm_grid))deallocate(r_atm_grid)
if(allocated(r_atm_grid_vec))deallocate(r_atm_grid_vec)
end subroutine hme

subroutine bldmatquadrupole(natom,atomcrd,ngrid,gridcrd,esp,rnorm,r,ibigmem,awork,bwork)
implicit none
integer(kind=4), intent(in) :: natom, ngrid
real(kind=8), intent(in) :: atomcrd(3,natom), gridcrd(3,ngrid)
real(kind=8), intent(in) :: esp(ngrid)
real(kind=8), intent(in) :: rnorm(ngrid,natom),r(3,ngrid,natom)
integer(kind=4), intent(in) :: ibigmem
real(kind=8), intent(out) :: awork(5*natom,5*natom), bwork(5*natom)
integer(kind=4) :: idx1, idx2
integer(kind=4) :: i1, j1, k1
integer(kind=4) :: i2, j2, k2
real(kind=8) :: rji1, rji2
real(kind=8) :: rji1_v(3), rji2_v(3)
real(kind=8) :: factor
integer(kind=4) :: i, j, k, m, n
if(ibigmem<4)then
  !do i = 1, natom
  !  do j = 1, ngrid
  !    do k = 1, 3
  !      r(k,j,i)=gridcrd(k,j)-atomcrd(k,i)
  !    end do
  !    rnorm(j,i)=sqrt(r(1,j,i)**2+r(2,j,i)**2+r(3,j,i)**2)
  !  end do
  !end do
  
  awork=0.d0
  bwork=0.d0
  do i1 = 1, natom
    idx1=5*(i1-1)+1
    do j = 1, ngrid
      do k = 1, 3
        rji1_v(k)=gridcrd(k,j)-atomcrd(k,i1)
      end do
      rji1=sqrt(rji1_v(1)**2+rji1_v(2)**2+rji1_v(3)**2)
      factor=(3.d0*rji1_v(3)**2-rji1**2)/(2.d0*rji1**5)
      do i2 = 1, natom
        do k = 1, 3
          rji2_v(k)=gridcrd(k,j)-atomcrd(k,i2)
        end do
        rji2=sqrt(rji2_v(1)**2+rji2_v(2)**2+rji2_v(3)**2)
        idx2=5*(i2-1)+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*(3.d0*rji2_v(3)**2-rji2**2)/(2.d0*rji2**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(1)*rji2_v(3)/rji2**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(2)*rji2_v(3)/rji2**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*0.5d0*(rji2_v(1)**2-rji2_v(2)**2)/(2.d0*rji2**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(1)*rji2_v(2)/rji2**5
      end do
      bwork(idx1)=bwork(idx1)+esp(j)*factor
    end do
  
    idx1=idx1+1
    do j = 1, ngrid
      do k = 1, 3
        rji1_v(k)=gridcrd(k,j)-atomcrd(k,i1)
      end do
      rji1=sqrt(rji1_v(1)**2+rji1_v(2)**2+rji1_v(3)**2)
      factor=sqrt(3.d0)*rji1_v(1)*rji1_v(3)/rji1**5
      do i2 = 1, natom
        do k = 1, 3
          rji2_v(k)=gridcrd(k,j)-atomcrd(k,i2)
        end do
        rji2=sqrt(rji2_v(1)**2+rji2_v(2)**2+rji2_v(3)**2)
        idx2=5*(i2-1)+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*(3.d0*rji2_v(3)**2-rji2**2)/(2.d0*rji2**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(1)*rji2_v(3)/rji2**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(2)*rji2_v(3)/rji2**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*0.5d0*(rji2_v(1)**2-rji2_v(2)**2)/(2.d0*rji2**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(1)*rji2_v(2)/rji2**5
      end do
      bwork(idx1)=bwork(idx1)+esp(j)*factor
    end do
  
    idx1=idx1+1
    do j = 1, ngrid
      do k = 1, 3
        rji1_v(k)=gridcrd(k,j)-atomcrd(k,i1)
      end do
      rji1=sqrt(rji1_v(1)**2+rji1_v(2)**2+rji1_v(3)**2)
      factor=sqrt(3.d0)*rji1_v(2)*rji1_v(3)/rji1**5
      do i2 = 1, natom
        do k = 1, 3
          rji2_v(k)=gridcrd(k,j)-atomcrd(k,i2)
        end do
        rji2=sqrt(rji2_v(1)**2+rji2_v(2)**2+rji2_v(3)**2)
        idx2=5*(i2-1)+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*(3.d0*rji2_v(3)**2-rji2**2)/(2.d0*rji2**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(1)*rji2_v(3)/rji2**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(2)*rji2_v(3)/rji2**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*0.5d0*(rji2_v(1)**2-rji2_v(2)**2)/(2.d0*rji2**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(1)*rji2_v(2)/rji2**5
      end do
      bwork(idx1)=bwork(idx1)+esp(j)*factor
    end do
  
    idx1=idx1+1
    do j = 1, ngrid
      do k = 1, 3
        rji1_v(k)=gridcrd(k,j)-atomcrd(k,i1)
      end do
      rji1=sqrt(rji1_v(1)**2+rji1_v(2)**2+rji1_v(3)**2)
      factor=sqrt(3.d0)*0.5d0*(rji1_v(1)**2-rji1_v(2)**2)/(2.d0*rji1**5)
      do i2 = 1, natom
        do k = 1, 3
          rji2_v(k)=gridcrd(k,j)-atomcrd(k,i2)
        end do
        rji2=sqrt(rji2_v(1)**2+rji2_v(2)**2+rji2_v(3)**2)
        idx2=5*(i2-1)+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*(3.d0*rji2_v(3)**2-rji2**2)/(2.d0*rji2**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(1)*rji2_v(3)/rji2**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(2)*rji2_v(3)/rji2**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*0.5d0*(rji2_v(1)**2-rji2_v(2)**2)/(2.d0*rji2**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(1)*rji2_v(2)/rji2**5
      end do
      bwork(idx1)=bwork(idx1)+esp(j)*factor
    end do
  
    idx1=idx1+1
    do j = 1, ngrid
      do k = 1, 3
        rji1_v(k)=gridcrd(k,j)-atomcrd(k,i1)
      end do
      rji1=sqrt(rji1_v(1)**2+rji1_v(2)**2+rji1_v(3)**2)
      factor=sqrt(3.d0)*rji1_v(1)*rji1_v(2)/rji1**5
      do i2 = 1, natom
        do k = 1, 3
          rji2_v(k)=gridcrd(k,j)-atomcrd(k,i2)
        end do
        rji2=sqrt(rji2_v(1)**2+rji2_v(2)**2+rji2_v(3)**2)
        idx2=5*(i2-1)+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*(3.d0*rji2_v(3)**2-rji2**2)/(2.d0*rji2**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(1)*rji2_v(3)/rji2**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(2)*rji2_v(3)/rji2**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*0.5d0*(rji2_v(1)**2-rji2_v(2)**2)/(2.d0*rji2**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*rji2_v(1)*rji2_v(2)/rji2**5
      end do
      bwork(idx1)=bwork(idx1)+esp(j)*factor
    end do
  end do
else
  awork=0.d0
  bwork=0.d0
  do i1 = 1, natom
    idx1=5*(i1-1)+1
    do j = 1, ngrid
      factor=(3.d0*r(3,j,i1)**2-rnorm(j,i1)**2)/(2.d0*rnorm(j,i1)**5)
      do i2 = 1, natom
        idx2=5*(i2-1)+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*(3.d0*r(3,j,i2)**2-rnorm(j,i2)**2)/(2.d0*rnorm(j,i2)**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(1,j,i2)*r(3,j,i2)/rnorm(j,i2)**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(2,j,i2)*r(3,j,i2)/rnorm(j,i2)**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*0.5d0*(r(1,j,i2)**2-r(2,j,i2)**2)/(2.d0*rnorm(j,i2)**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(1,j,i2)*r(2,j,i2)/rnorm(j,i2)**5
      end do
      bwork(idx1)=bwork(idx1)+esp(j)*factor
    end do
  
    idx1=idx1+1
    do j = 1, ngrid
      factor=sqrt(3.d0)*r(1,j,i1)*r(3,j,i1)/rnorm(j,i1)**5
      do i2 = 1, natom
        idx2=5*(i2-1)+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*(3.d0*r(3,j,i2)**2-rnorm(j,i2)**2)/(2.d0*rnorm(j,i2)**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(1,j,i2)*r(3,j,i2)/rnorm(j,i2)**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(2,j,i2)*r(3,j,i2)/rnorm(j,i2)**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*0.5d0*(r(1,j,i2)**2-r(2,j,i2)**2)/(2.d0*rnorm(j,i2)**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(1,j,i2)*r(2,j,i2)/rnorm(j,i2)**5
      end do
      bwork(idx1)=bwork(idx1)+esp(j)*factor
    end do
  
    idx1=idx1+1
    do j = 1, ngrid
      factor=sqrt(3.d0)*r(2,j,i1)*r(3,j,i1)/rnorm(j,i1)**5
      do i2 = 1, natom
        idx2=5*(i2-1)+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*(3.d0*r(3,j,i2)**2-rnorm(j,i2)**2)/(2.d0*rnorm(j,i2)**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(1,j,i2)*r(3,j,i2)/rnorm(j,i2)**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(2,j,i2)*r(3,j,i2)/rnorm(j,i2)**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*0.5d0*(r(1,j,i2)**2-r(2,j,i2)**2)/(2.d0*rnorm(j,i2)**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(1,j,i2)*r(2,j,i2)/rnorm(j,i2)**5
      end do
      bwork(idx1)=bwork(idx1)+esp(j)*factor
    end do
  
    idx1=idx1+1
    do j = 1, ngrid
      factor=sqrt(3.d0)*0.5d0*(r(1,j,i1)**2-r(2,j,i1)**2)/(2.d0*rnorm(j,i1)**5)
      do i2 = 1, natom
        idx2=5*(i2-1)+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*(3.d0*r(3,j,i2)**2-rnorm(j,i2)**2)/(2.d0*rnorm(j,i2)**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(1,j,i2)*r(3,j,i2)/rnorm(j,i2)**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(2,j,i2)*r(3,j,i2)/rnorm(j,i2)**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*0.5d0*(r(1,j,i2)**2-r(2,j,i2)**2)/(2.d0*rnorm(j,i2)**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(1,j,i2)*r(2,j,i2)/rnorm(j,i2)**5
      end do
      bwork(idx1)=bwork(idx1)+esp(j)*factor
    end do
  
    idx1=idx1+1
    do j = 1, ngrid
      factor=sqrt(3.d0)*r(1,j,i1)*r(2,j,i1)/rnorm(j,i1)**5
      do i2 = 1, natom
        idx2=5*(i2-1)+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*(3.d0*r(3,j,i2)**2-rnorm(j,i2)**2)/(2.d0*rnorm(j,i2)**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(1,j,i2)*r(3,j,i2)/rnorm(j,i2)**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(2,j,i2)*r(3,j,i2)/rnorm(j,i2)**5
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*0.5d0*(r(1,j,i2)**2-r(2,j,i2)**2)/(2.d0*rnorm(j,i2)**5)
        idx2=idx2+1
        awork(idx1,idx2)=awork(idx1,idx2)+factor*sqrt(3.d0)*r(1,j,i2)*r(2,j,i2)/rnorm(j,i2)**5
      end do
      bwork(idx1)=bwork(idx1)+esp(j)*factor
    end do
  end do

end if 
end subroutine bldmatquadrupole




subroutine bldmatdipole(natom,atomcrd,ngrid,gridcrd,esp,rnorm,r,ibigmem,awork,bwork)
implicit none
integer(kind=4), intent(in) :: natom, ngrid
real(kind=8), intent(in) :: atomcrd(3,natom), gridcrd(3,ngrid)
real(kind=8), intent(in) :: esp(ngrid)
real(kind=8), intent(in) :: rnorm(ngrid,natom),r(3,ngrid,natom)
integer(kind=4), intent(in) :: ibigmem
real(kind=8), intent(out) :: awork(3*natom,3*natom), bwork(3*natom)
integer(kind=4) :: idx1, idx2
integer(kind=4) :: i1, j1, k1
integer(kind=4) :: i2, j2, k2
real(kind=8) :: rji1, rji2
integer(kind=4) :: i, j, k, m, n
if(ibigmem==4)then
  awork=0.d0
  bwork=0.d0
  idx1=0
  do i1 = 1, natom
    do k1 = 1, 3
       idx1=idx1+1
       idx2=0
       do i2 = 1, natom
         do k2 = 1, 3
           idx2=idx2+1
           do j = 1, ngrid
             awork(idx1,idx2)=awork(idx1,idx2)+r(k1,j,i1)*r(k2,j,i2)/(rnorm(j,i1)*rnorm(j,i2))**3
           end do
         end do
       end do
       do j = 1, ngrid
         bwork(idx1)=bwork(idx1)+esp(j)*r(k1,j,i1)/rnorm(j,i1)**3
       end do
    end do
  end do
else
  awork=0.d0
  bwork=0.d0
  idx1=0
  do i1 = 1, natom
    do k1 = 1, 3
       idx1=idx1+1
       idx2=0
       do i2 = 1, natom
         do k2 = 1, 3
           idx2=idx2+1
           do j = 1, ngrid
             rji1=sqrt((gridcrd(1,j)-atomcrd(1,i1))**2+(gridcrd(2,j)-atomcrd(2,i1))**2+(gridcrd(3,j)-atomcrd(3,i1))**2)
             rji2=sqrt((gridcrd(1,j)-atomcrd(1,i2))**2+(gridcrd(2,j)-atomcrd(2,i2))**2+(gridcrd(3,j)-atomcrd(3,i2))**2)
             awork(idx1,idx2)=awork(idx1,idx2)+(gridcrd(k1,j)-atomcrd(k1,i1))*(gridcrd(k2,j)-atomcrd(k2,i2))/(rji1*rji2)**3
           end do
         end do
       end do
       do j = 1, ngrid
         rji1=sqrt((gridcrd(1,j)-atomcrd(1,i1))**2+(gridcrd(2,j)-atomcrd(2,i1))**2+(gridcrd(3,j)-atomcrd(3,i1))**2)
         bwork(idx1)=bwork(idx1)+esp(j)*(gridcrd(k1,j)-atomcrd(k1,i1))/rji1**3
       end do
    end do
  end do
end if
end subroutine bldmatdipole


subroutine qualitymetric(ngrid,esp0,esp1,sstot,ssreg,ssres,crosscorr)
implicit none
integer(kind=4), intent(in) :: ngrid
real(kind=8), intent(in) :: esp0(ngrid), esp1(ngrid)
real(kind=8), intent(out) :: sstot, ssreg, ssres, crosscorr
real(kind=8) :: esp0mean, esp1mean
real(kind=8) :: esp0width, esp1width
integer(kind=4) :: i, j, k, m, n
esp0mean = sum(esp0)/max(1,ngrid)
esp1mean = sum(esp1)/max(1,ngrid)
sstot = 0.d0
ssreg = 0.d0
ssres = 0.d0
do i = 1, ngrid
  sstot = sstot + (esp0(i)-esp0mean)**2
  ssreg = ssreg + (esp1(i)-esp0mean)**2
  ssres = ssres + (esp1(i)-esp0(i))**2
end do
crosscorr=0.d0
esp0width=0.d0
esp1width=0.d0
do i = 1, ngrid
  crosscorr=crosscorr+(esp0(i)-esp0mean)*(esp1(i)-esp1mean)
  esp0width=esp0width+(esp0(i)-esp0mean)**2
  esp1width=esp1width+(esp1(i)-esp1mean)**2
end do
crosscorr = crosscorr/sqrt(esp0width*esp1width+1.D-10)
end subroutine qualitymetric



subroutine multipoleESP(ngrid,gridcrd,natom,atomcrd,monopole,dipole,quadrupole,esp)
implicit none
integer(kind=4), intent(in) :: ngrid, natom
real(kind=8), intent(in) :: gridcrd(3,ngrid), atomcrd(3,natom), monopole(natom), dipole(3,natom), quadrupole(5,natom)
real(kind=8), intent(out) :: esp(ngrid)
real(kind=8) :: thisdipole(3), thisquadrupole(5)
real(kind=8) :: dr(3), rnorm
real(kind=8) :: dipoleesp, quadrupoleesp
integer(kind=4) :: i, j, k, m, n
esp=0.d0
do i =1, ngrid
  do j = 1, natom
    do k =1, 3
      dr(k) = gridcrd(k,i)-atomcrd(k,j)
    end do
    rnorm=sqrt(dot_product(dr,dr))
    thisdipole(:)=dipole(:,j)
    thisquadrupole(:)=quadrupole(:,j)
    esp(i)=esp(i)+monopole(j)/rnorm
    esp(i)=esp(i)+dipoleesp(thisdipole,dr)
    esp(i)=esp(i)+quadrupoleesp(thisquadrupole,dr)
  end do
end do
end subroutine multipoleESP


function dipoleesp(dipole,r)
implicit none
real(kind=8), intent(in) :: dipole(3),r(3)
real(kind=8) :: dipoleesp
real(kind=8) :: rnorm
rnorm=sqrt(dot_product(r,r))
dipoleesp=dot_product(dipole,r)/rnorm**3
end function dipoleesp



function quadrupoleesp(quadrupole,r)
implicit none
real(kind=8),intent(in) :: quadrupole(5),r(3)
real(kind=8) :: quadrupoleesp
real(kind=8) :: rnorm
rnorm=sqrt(dot_product(r,r))
quadrupoleesp=(quadrupole(1)*(1.5d0*r(3)**2-0.5d0*rnorm**2) &
              +quadrupole(2)*sqrt(3.d0)*r(1)*r(3) &
              +quadrupole(3)*sqrt(3.d0)*r(2)*r(3) &
              +quadrupole(4)*0.5d0*sqrt(3.d0)*(r(1)**2-r(2)**2) &
              +quadrupole(5)*sqrt(3.d0)*r(1)*r(2) &
              )/rnorm**5
end function quadrupoleesp


subroutine checkconvergence(natom,q,q0,criterion,isconverged)
implicit none
integer(kind=4), intent(in) :: natom
real(kind=8), intent(in) :: q(natom), q0(natom)
real(kind=8), intent(in) :: criterion
logical, intent(out) :: isconverged
integer(kind=4) :: i, j, k, m, n
isconverged=.true.
do i =1, natom
  if(abs(q(i)-q0(i))>criterion)then
    isconverged=.false.
    return
  end if
end do
end subroutine checkconvergence

subroutine getrestraints(natom,nlgn,awork,bwork,istrnt,strength,q0,qwt)
implicit none
integer(kind=4), intent(in)::natom,nlgn
real(kind=8), intent(in out)::awork(natom+nlgn,natom+nlgn),bwork(natom+nlgn)
integer(kind=4), intent(in)::istrnt
real(kind=8), intent(in)::strength
real(kind=8), intent(in)::q0(natom),qwt(natom)
integer(kind=4)istrntlocal
integer(kind=4)::i,j,k,m,n
istrntlocal=istrnt
if(sum(abs(q0))<1.D-8)istrntlocal=0
if(istrntlocal==0)then
  do i = 1, natom
    awork(i,i)=awork(i,i)+strength*qwt(i)
    bwork(i)=bwork(i)+strength*qwt(i)*q0(i)
  end do
else
  do i = 1, natom
    awork(i,i)=awork(i,i)+strength*qwt(i)/sqrt(q0(i)**2+0.01d0) 
  end do
end if
end subroutine getrestraints

subroutine LELUD(m,n,a,b,x,outid)
  implicit none
  integer(kind=4), intent(in) :: m, n
  real(kind=8), intent(in) :: a(m,n), b(m)
  real(kind=8), intent(out) :: x(n)
  integer(kind=4), intent(in) :: outid
  real(kind=8) :: acopy(m,n), bcopy(m)
  integer(kind=4) :: ipiv(n)
  integer(kind=4) :: info
  integer(kind=4) :: i
  if(m<n) then
    write(outid,'(A)')'The lower dimension (m) of a must be no smaller than the order of A (n)'
    stop
  end if
  acopy = a
  call dgetrf( m, n, acopy, m, ipiv, info )
  if(info/=0)then
    write(outid,'(A)') 'dgetrf failed'
    stop
  end if
  bcopy = b
  call dgetrs( 'N', n, 1, acopy, m, ipiv, bcopy, m, info )
  if(info/=0)then
    write(outid,'(A)') 'dgetrs failed'
    stop
  end if
  forall(i=1:n)
    x(i)=bcopy(i)
  end forall
end subroutine LELUD

subroutine LESVD(m,n,a,b,x,outid)
implicit none
integer(kind=4), intent(in) :: m, n
real(kind=8), intent(in) :: a(m,n), b(m)
real(kind=8), intent(out) :: x(n)
integer(kind=4), intent(in) :: outid
integer(kind=4) :: maxmn, minmn
real(kind=8), allocatable :: sigma(:), sigma_inv(:)
real(kind=8) :: u(m,m), vt(n,n)
integer(kind=4) :: ldu, ldvt
integer(kind=4) :: lwork
real(kind=8), allocatable :: work(:)
integer(kind=4) :: info
real(kind=8) :: utu(m,m), vvt(n,n), vs(n,m)
real(kind=8), allocatable :: vsut(:,:)
real(kind=8) :: acopy(m,n)
integer(kind=4) :: i, j, k
maxmn=m
if(n>m)maxmn=n
minmn=m
if(n<m)minmn=n
lwork=10*maxmn
if(allocated(sigma))deallocate(sigma)
allocate(sigma(maxmn))
if(allocated(work))deallocate(work)
allocate(work(lwork))
ldu=m
ldvt=n
acopy=a
call dgesvd('A','A',m,n,acopy,m,sigma,u,ldu,vt,ldvt,work,lwork,info)
if(info<0)then
  write(outid,*)'In LESVD dgesvd failed with info=',info
  stop
else if(info>0) then
  write(outid,*)'In LESVD dgesvd failed with info=',info
  stop
else
!  write(outid,*)'Checking SVD '
!  utu=0.d0
!  do i = 1, m
!    do j = 1, m
!       do k = 1, m
!          utu(i,j)=utu(i,j)+u(k,i)*u(k,j)
!       end do
!    end do
!  end do
!  vvt=0.d0
!  do i = 1, n
!    do j = 1, n
!       do k = 1, n
!          vvt(i,j)=vvt(i,j)+vt(k,i)*vt(k,j)
!       end do
!    end do
!  end do
!  do i = 1, m
!    write(outid,'(<m>(F12.3,1X))')utu(i,:)
!  end do
!  write(outid,*)
!  do i = 1, n
!    write(outid,'(<n>(F12.3,1X))')vvt(i,:)
!  end do
  write(outid,*)'Singular value of least-square matrix:'
  write(outid,'(8F10.4)')sigma
  write(outid,*)
  if(allocated(sigma_inv))deallocate(sigma_inv)
  allocate(sigma_inv(maxmn))
  do i = 1, minmn
    sigma_inv(i)=1.d0/sigma(i)
  end do
  do i = minmn+1, maxmn
    sigma_inv(i)=0.d0
  end do
  if(allocated(vsut))deallocate(vsut)
  allocate(vsut(n,m))
  do i = 1, n
    do j =1, n
      vs(i,j)=vt(j,i)*sigma_inv(j)
    end do
  end do
  vsut=0.d0
  do i = 1, n
    do j = 1, m
      do k = 1, m 
        vsut(i,j)=vsut(i,j)+vs(i,k)*u(j,k)
      end do
    end do
  end do
  x=0.d0
  do i = 1, n
    do j = 1, m
      x(i)=x(i)+vsut(i,j)*b(j)
    end do
  end do 
end if
if(allocated(sigma))deallocate(sigma)
if(allocated(work))deallocate(work)
if(allocated(sigma_inv))deallocate(sigma_inv)
if(allocated(vsut))deallocate(vsut)
end subroutine LESVD


subroutine sortsymmconstraint(nsymmconstraint,isymmconstraintpair)
implicit none
integer(kind=4), intent(in)::nsymmconstraint
integer(kind=4), intent(in out)::isymmconstraintpair(2,nsymmconstraint)
real(kind=8)::isymmconstraintpair_sorted(2,nsymmconstraint)
integer(kind=4)::lastsite(nsymmconstraint)
integer(kind=4)::location(nsymmconstraint)
integer(kind=4)::loctmp(1)
integer(kind=4)::iswaptmp
integer(kind=4)::i,j,k,m,n
lastsite=0
do i = 1, nsymmconstraint
  if(isymmconstraintpair(2,i)<isymmconstraintpair(1,i))then
    iswaptmp=isymmconstraintpair(2,i)
    isymmconstraintpair(2,i)=isymmconstraintpair(1,i)
    isymmconstraintpair(1,i)=iswaptmp
  end if
  lastsite(i)=isymmconstraintpair(2,i)
end do
isymmconstraintpair_sorted=0
do i = 1, nsymmconstraint
  loctmp=maxloc(lastsite)
  location(i)=loctmp(1)
  lastsite(location(i))=-1
end do

do i = 1, nsymmconstraint
  isymmconstraintpair_sorted(:,i)=isymmconstraintpair(:,location(i))
end do
isymmconstraintpair=isymmconstraintpair_sorted
end subroutine sortsymmconstraint

subroutine rebuildchargearray(natom,q,nsymmconstraint,isymmconstraintpair,iskilledsite)
implicit none
integer(kind=4),intent(in)::natom
real(kind=8),intent(in out)::q(natom)
integer(kind=4),intent(in)::nsymmconstraint
integer(kind=4),intent(in)::isymmconstraintpair(2,nsymmconstraint)
integer(kind=4),intent(in)::iskilledsite(natom)
real(kind=8)::qfull(natom)
integer(kind=4)::i,j,k,m,n
j=0
do i = 1, natom
  if(iskilledsite(i)==0)then
    j=j+1
    qfull(i)=q(j)
  else
    qfull(i)=qfull(isymmconstraintpair(1,iskilledsite(i)))
  end if
end do
q(1:natom)=qfull(1:natom)
end subroutine rebuildchargearray

subroutine applysymmetry(natom,nlgn,awork,bwork,nsymmconstraint,isymmconstraintpair,iskilledsite)
implicit none
integer(kind=4),intent(in)::natom,nlgn
real(kind=8),intent(inout)::awork(natom+nlgn,natom+nlgn),bwork(natom+nlgn)
integer(kind=4),intent(in)::nsymmconstraint
integer(kind=4),intent(in out)::isymmconstraintpair(2,nsymmconstraint)
integer(kind=4),intent(out)::iskilledsite(natom)
integer(kind=4)::iatom2bekilled
real(kind=8)::column2bekilled(natom+nlgn),row2bekilled(natom+nlgn)
integer(kind=4)::i,j,k,m,n
do i = 1, nsymmconstraint
  iatom2bekilled=isymmconstraintpair(2,i)
  row2bekilled(:)=awork(iatom2bekilled,:)
  awork(isymmconstraintpair(1,i),:)=awork(isymmconstraintpair(1,i),:)+row2bekilled(:)
  column2bekilled(:)=awork(:,iatom2bekilled)
  awork(:,isymmconstraintpair(1,i))=awork(:,isymmconstraintpair(1,i))+column2bekilled(:)
  bwork(isymmconstraintpair(1,i))=bwork(isymmconstraintpair(1,i))+bwork(iatom2bekilled)
  awork(:,iatom2bekilled)=0.d0
  awork(iatom2bekilled,:)=0.d0
  iskilledsite(iatom2bekilled)=i
end do 
end subroutine applysymmetry


subroutine appendlgn2mat(natom,nlgn,awork,bwork,lgnatomidx,lgncons)
implicit none
integer(kind=4),intent(in)::natom,nlgn
real(kind=8),intent(inout)::awork(natom+nlgn,natom+nlgn),bwork(natom+nlgn)
integer(kind=4),intent(in)::lgnatomidx(natom,nlgn)
real(kind=8),intent(in)::lgncons(nlgn)
integer(kind=4)::aworkrowidx
integer(kind=4)::i,j,k,m,n
do i = 1, nlgn
  aworkrowidx=i+natom
  do j = 1, natom
    awork(aworkrowidx,j)=lgnatomidx(j,i)
    awork(j,aworkrowidx)=lgnatomidx(j,i)
  end do
  bwork(aworkrowidx)=lgncons(i)
end do
end subroutine appendlgn2mat





subroutine gendismat(natom,ngrid,atomcrd,gridcrd,r_atm_grid,r_atm_grid_vec,ibigmem,outid)
implicit none
integer(kind=4),intent(in)::natom,ngrid
real(kind=8),intent(in)::atomcrd(3,natom),gridcrd(3,ngrid)
real(kind=8),allocatable,intent(out)::r_atm_grid(:,:), r_atm_grid_vec(:,:,:)
integer(kind=4),intent(out)::ibigmem
integer(kind=4),intent(in)::outid
integer(kind=4)::AllocateStatus
integer(kind=4)::i,j,k,m,n
ibigmem=0
write(outid,*)'Allocate memory for distance matrix between ATOMS and GRIDS'
write(outid,*)'It takes ',natom*ngrid*4*8, ' bytes'
if(allocated(r_atm_grid))deallocate(r_atm_grid)
allocate(r_atm_grid(ngrid,natom), STAT = AllocateStatus)
if(AllocateStatus /= 0)then
  write(outid,*)'No enough memory for distance matrix'
  write(outid,*)'Switching to small memory mode'
  write(outid,*)
  ibigmem=0
  return
else
  write(outid,*)'Distance matrix allocated'
  write(outid,*)
  ibigmem=1
  if(allocated(r_atm_grid_vec))deallocate(r_atm_grid_vec)
  allocate(r_atm_grid_vec(3,ngrid,natom), STAT = AllocateStatus)
  if(AllocateStatus /= 0)then
    write(outid,*)'No enough memory for distance vector matrix'
    write(outid,*)'Switching to median memory mode'
    do i =1, ngrid
      do j =1, natom
        r_atm_grid(i,j)=(atomcrd(1,j)-gridcrd(1,i))**2 &
                       +(atomcrd(2,j)-gridcrd(2,i))**2 &
                       +(atomcrd(3,j)-gridcrd(3,i))**2 
      end do
    end do
    r_atm_grid=sqrt(r_atm_grid)
  else
    write(outid,*)'Distance vector matrix allocated'
    write(outid,*)
    ibigmem=4
    do i =1, ngrid
      do j =1, natom
        do m = 1, 3
          r_atm_grid_vec(m,i,j)=gridcrd(m,i)-atomcrd(m,j)
        end do
        r_atm_grid(i,j)=(r_atm_grid_vec(1,i,j)**2 &
                       + r_atm_grid_vec(2,i,j)**2 &
                       + r_atm_grid_vec(3,i,j)**2)
      end do
    end do
    r_atm_grid=sqrt(r_atm_grid)
  end if
end if
end subroutine gendismat




subroutine bldmatmonopole(natom,nlgn,atomcrd,ngrid,gridcrd,esp,r_atm_grid,ibigmem,awork,bwork)
implicit none
integer(kind=4), intent(in) :: natom, nlgn, ngrid
real(kind=8), intent(in) :: atomcrd(3,natom), gridcrd(3,ngrid), esp(ngrid)
real(kind=8), intent(in) :: r_atm_grid(ngrid,natom)
integer(kind=4), intent(in) :: ibigmem
real(kind=8), intent(out) :: awork(natom+nlgn,natom+nlgn), bwork(natom+nlgn)
real(kind=8)::rik,rjk
integer(kind=4)::i,j,k,m,n
awork=0.d0
bwork=0.d0
if(ibigmem>=1)then
  do i = 1, natom
    do j = i+1, natom
      do k = 1, ngrid
        awork(i,j)=awork(i,j)+1.d0/(r_atm_grid(k,i)*r_atm_grid(k,j))
      end do
      awork(j,i)=awork(i,j)
    end do
    do k = 1, ngrid
      awork(i,i)=awork(i,i)+1.d0/(r_atm_grid(k,i)*r_atm_grid(k,i))
      bwork(i)=bwork(i)+esp(k)/r_atm_grid(k,i)
    end do
  end do
else
  do i = 1, natom
    do j = i+1, natom
      do k = 1, ngrid
        rik=sqrt((atomcrd(1,i)-gridcrd(1,k))**2+(atomcrd(2,i)-gridcrd(2,k))**2+(atomcrd(3,i)-gridcrd(3,k))**2)
        rjk=sqrt((atomcrd(1,j)-gridcrd(1,k))**2+(atomcrd(2,j)-gridcrd(2,k))**2+(atomcrd(3,j)-gridcrd(3,k))**2)
        awork(i,j)=awork(i,j)+1.d0/(rik*rjk)
      end do
      awork(j,i)=awork(i,j)
    end do
    do k = 1, ngrid
      rik=sqrt((atomcrd(1,i)-gridcrd(1,k))**2+(atomcrd(2,i)-gridcrd(2,k))**2+(atomcrd(3,i)-gridcrd(3,k))**2)
      awork(i,i)=awork(i,i)+1.d0/rik**2
      bwork(i)=bwork(i)+esp(k)/rik
    end do
  end do
end if
end subroutine bldmatmonopole

subroutine getlgn(natom,nlgn0,lgnatomidx0,lgncons0,nlgn,lgnatomidx,lgncons,qresidual,qbase,outid)
implicit none
integer(kind=4),intent(in)::natom
integer(kind=4),intent(in)::nlgn0
integer(kind=4),intent(in)::lgnatomidx0(natom,nlgn0)
real(kind=8),intent(in)::lgncons0(nlgn0)
integer(kind=4),intent(out)::nlgn
integer(kind=4),allocatable,intent(out)::lgnatomidx(:,:)
real(kind=8),allocatable,intent(out)::lgncons(:)
real(kind=8),intent(in)::qresidual
real(kind=8),intent(in)::qbase(natom)
integer(kind=4),intent(in)::outid
real(kind=8)::sumq
integer(kind=4)::i,j,k,m,n
nlgn=nlgn0+1
allocate(lgnatomidx(natom,nlgn))
allocate(lgncons(nlgn))
do i = 1, nlgn0
  do j = 1, natom
    lgnatomidx(j,i)=lgnatomidx0(j,i)
  end do
  lgncons(i)=lgncons0(i)
end do
lgnatomidx(:,nlgn)=1
lgncons(nlgn)=qresidual
do i = 1, nlgn0
  sumq=0.d0
  do j = 1, natom
    sumq=sumq+lgnatomidx(j,i)*qbase(j)
  end do
  lgncons(i)=lgncons(i)-sumq
end do
write(outid,'(A,I5)')' Total number of Lagrange conditions:', nlgn
do i = 1, nlgn
  write(outid,'(<natom>I3,A,F8.4)')lgnatomidx(:,i), ' & ', lgncons(i)
end do
write(outid,*)
end subroutine getlgn

subroutine initialize(natom,ngrid,gridcrd,atomcrd,espinit,espresidual,mpinit,qbase,qwt,r_atm_grid,ibigmem,fitmethod,outid)
implicit none
integer(kind=4),intent(in)::natom,ngrid
real(kind=8),intent(in)::gridcrd(3,ngrid),atomcrd(3,natom)
real(kind=8),intent(in)::espinit(ngrid)
real(kind=8),intent(out)::espresidual(ngrid)
real(kind=8),intent(in)::mpinit(natom)
real(kind=8),intent(out)::qbase(natom)
real(kind=8),intent(out)::qwt(natom)
real(kind=8),intent(in)::r_atm_grid(ngrid,natom)
integer(kind=4),intent(in)::ibigmem
character(len=*),intent(in)::fitmethod
integer(kind=4),intent(in)::outid
real(kind=8)::totalwt
real(kind=8)::r
real(kind=8),parameter::machinezero=1.D-4
integer(kind=4)::i,j,k,m,n
espresidual=espinit  
qwt=1.0d0
qbase=0.d0
if(fitmethod=='DRESP')then
  qbase=mpinit
  do i = 1, natom
    qwt(i)=1.d0/(mpinit(i)**2+machinezero)
  end do
  totalwt=1.d0/((sum(abs(mpinit(1:natom)))/natom)**2+machinezero)
  qwt=qwt/totalwt
!  write(40,*)'Totalwt=',totalwt
  if(ibigmem>=1)then
    ! r_atm_grid is available
    do i = 1, natom
      do j = 1, ngrid
        espresidual(j)=espresidual(j)-qbase(i)/r_atm_grid(j,i)
      end do
    end do
  else
    ! r_atm_grid is NOT available
    do i = 1, natom
      do j = 1, ngrid
        r=sqrt((atomcrd(1,i)-gridcrd(1,j))**2+(atomcrd(2,i)-gridcrd(2,j))**2+(atomcrd(3,i)-gridcrd(3,j))**2)
        espresidual(j)=espresidual(j)-mpinit(i)/r
      end do
    end do
  end if
  write(outid,'(A)')' Weighting factor for each site:'
  write(outid,'(6E12.3)')qwt
  write(outid,*)
end if
end subroutine initialize
