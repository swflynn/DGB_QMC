!=============================================================================80
!        1D-Distributed Gaussian Basis (Coordinate Space Formulation)
!==============================================================================!
!    Discussion:
!DGB analysis for 1D potential: single atom (x coordinate only, mass:=1)
!Equations are formulated assuming Cartesian Coordinates
!Gaussian Centers selected with Uniform Grid (same widths;alpha parameter)
!==============================================================================!
!    Modified:
!       17 September 2018
!    Author:
!       Shane Flynn 
!==============================================================================!
module dgb_groundstate
implicit none
!==============================================================================!
!                           Global Paramaters
!==============================================================================!
double precision,parameter::Hmass=1d0
double precision,parameter::pi=4.*atan(1d0)
!==============================================================================!
!                           Global Variables 
!==============================================================================!
!       Discussion:
!Natoms             ==> Number of atoms 
!Dimen              ==> System Dimensionality
!==============================================================================!
integer::Natoms,Dimen
character(len=2),allocatable::atom_type(:)
double precision,allocatable::mass(:),sqrt_mass(:)
!==============================================================================!
contains
!==============================================================================!
function Atom_Mass(atom)
!==============================================================================!
implicit none
double precision::Atom_Mass
character(len=2)::atom
if(atom=='H'.or.atom=='h')then 
    Atom_mass=Hmass
else 
    write(*,*) 'atom ', atom, ' is not recognized'
    stop 'Check Atom_Mass Function' 
endif
end function Atom_Mass
!==============================================================================!
subroutine Toy_Potential(x,energies)
!==============================================================================!
!       Discussion:
!Hard-Coded Potential Energy (1-Dimensional)
!V:=1/2*(x)^2 
!==============================================================================!
implicit none
double precision::x(Dimen),energies
energies=0d0
energies=0.5*(x(1))**2 
end subroutine Toy_Potential
!==============================================================================!
subroutine Toy_Force(x,forces)
!==============================================================================!
!       Discussion:
!Returns the Forces associated with Toy_Potential Subroutine
!Forces are hard-coded based on Toy_Potential
!==============================================================================!
implicit none
integer::i
double precision::x(Dimen),forces(Dimen)
forces(1)=-x(1) 
end subroutine Toy_Force
!==============================================================================!
subroutine Toy_Hessian(x,Hess_Mat)
!==============================================================================!
!       Discussion:
!Numerically computed Hessian using forces from Toy_Force Subroutine
!s          ==> Perturbation Parameter 
!Hess_Mat   ==> (Dimen,Dimen); Symmetrized Mass-Scaled Hessian
!x          ==> (Dimen); XYZ Configuration at minimum
!==============================================================================!
implicit none 
integer::i,j
double precision::Hess_Mat(Dimen,Dimen),x(Dimen),r(Dimen),force(Dimen)
double precision::force0(Dimen)
double precision,parameter::s=1d-6
r=x
call Toy_Force(r, force0)
do i=1,Dimen
    r(i)=x(i)+s
    call Toy_Force(r,force)
    r(i)=x(i)
    do j=1,Dimen
        Hess_Mat(i,j)=(force0(j)-force(j))/s
    enddo
enddo
!==============================================================================!
!                   Symmetrize and Mass-Scale the Hessian
!==============================================================================!
do i=1,Dimen
    do j=1,i
        if(i.ne.j) Hess_Mat(i,j)=(Hess_Mat(i,j)+Hess_Mat(j,i))/2
        Hess_Mat(i,j)=Hess_Mat(i,j)/(sqrt_mass(i)*sqrt_mass(j))
        if(i.ne.j) Hess_Mat(j,i)=Hess_Mat(i,j)
    enddo
enddo
end subroutine Toy_Hessian
!==============================================================================!
subroutine Frequencies_From_Hess(Dimen,Hess,omega,U)
!==============================================================================!
!       Discussion:
!Compute Eigenvalues and Eigenvectors of Hessian
!Uses the LLAPACK real symmetric eigen-solver (dsygev)
!Hess   ==> (Dimen,Dimen); Hessian Matrix
!omega  ==> (Dimen); Eigenvalues of the Hessian
!U      ==> (Dimen,Dimen); Eigenvectors of the Hessian
!       LLAPACK (dsyev):
!v      ==> Compute both Eigenvalues and Eigenvectors
!u      ==> Use Upper-Triangle of matrix
!==============================================================================!
implicit none
integer::i,info,lwork,Dimen
double precision::Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen)
double precision,allocatable::work(:) 
lwork=max(1,3*Dimen-1)        !suggested by LAPACK Developers
allocate(work(max(1,lwork)))  !suggested by LAPACK Developers 
U=Hess  
call dsyev('v','u',Dimen,U,Dimen,omega,work,lwork,info) 
write(*,*) 'Frequencies from the Hessian:'
do i=Dimen,1,-1
    omega(i)=sign(sqrt(abs(omega(i))),omega(i))
    write(*,*) omega(i), 'normalized = 1?',sum(U(:,i)**2)
enddo
end subroutine frequencies_from_Hess
!==============================================================================!
end module dgb_groundstate
!==============================================================================!
!==============================================================================!
program DGB_ground_state
use dgb_groundstate
!==============================================================================!
!==============================================================================!
!               Discussion:
!coord_in       ==> Input Geometry File-Name (xyz file)
!NG             ==> Number of Gaussians
!Nsobol         ==> Number of Sobol Points for numerical integration
!ii,skip        ==> Integer(kind=8): necessary for using Sobol Module
!alpha          ==> (NG): Width of Gaussian
!q0             ==> (Dimen): Input Geometry coordinates
!q              ==> (Dimen,NG): Gaussian Centers
!q2             ==> (Dimen,Nsobol): Sequence for Potential Integration
!force          ==> (Dimen): Forces. See "Toy_Force"
!Hess           ==> (Dimen,Dimen): Mass-Scaled Hessian. See "Toy_Hessian"
!omega          ==> (Dimen): Frequencies (Eigenvalues). See "Freq_Hess"
!U              ==> (Dimen,Dimen): Hess Eigenvectors. See "Freq_Hess"
!z              ==> (Dimen,NG): Quasi-Random Sequence for Gaussian Centers
!z2             ==> (Dimen,Nsobol): Quasi-Random Sequence for Integration
!Smat           ==> (NG,NG): Overlap Matrix
!Tmat           ==> (NG,NG): Kinetic Energy Matrix
!Vmat           ==> (NG,NG): Potential Energy Matrix
!Hmat           ==> (NG,NG): Hamiltonian Matrix (V+T) for eigenvalue problem
!eigenvalues    ==> (NG): Eigenvalues of the Hamiltonian Matrix
!lambda         ==> (NG): Eigenvalues of the Overlap Matrix
!E0             ==> Energy Evaluation at the minimum configuration
!pot_ene        ==> Potential Energy evaluated in q space
!lmat           ==> (NG,NG): Lambda Matrix 
!==============================================================================!
!==============================================================================!
implicit none
character(len=50)::coord_in
integer::NG,Nsobol
integer::i,j,k,n
integer*8::ii,skip,skip2
double precision::E0,lsum,prod_omega,Tsum,pot_ene,alpha_par,Tpre,low_bound
!==============================================================================!
double precision,allocatable::q0(:),force(:),q(:,:),q2(:,:),Hess(:,:)
double precision,allocatable::omega(:),U(:,:),z(:,:),z2(:,:),alpha(:)
double precision,allocatable::lmat(:,:),Smat(:,:),S1mat(:,:),eigenvalues(:)
double precision,allocatable::Tmat(:,:),Vmat(:,:),Hmat(:,:)
!==============================================================================!
!                           dsygv variables                                    !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable::work(:)
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) coord_in
read(*,*) NG
read(*,*) Nsobol
read(*,*) alpha_par
read(*,*) low_bound
skip=NG
skip2=Nsobol
write(*,*) 'Test 1; Successfully Read Input Data File'
!==============================================================================!
!                         Input Water Geometry 
!==============================================================================!
open(16,File=coord_in)
read(16,*) Natoms
read(16,*) 
Dimen=1*Natoms 
write(*,*) 'Dimensionality ==> ', Dimen
!==============================================================================!
!                   Input Configuration Energy
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(Dimen),q0(Dimen),force(Dimen))
allocate(Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen),q(Dimen,NG),z(Dimen,NG))
allocate(z2(Dimen,Nsobol),alpha(NG),lmat(NG,NG),Smat(NG,NG),S1mat(NG,NG))
allocate(eigenvalues(NG),Tmat(NG,NG),Vmat(NG,NG),Hmat(NG,NG),q2(Dimen,Nsobol))
do i=1,Natoms
    read(16,*) atom_type(i), q0(1)   
    mass(i)=Atom_mass(atom_type(i))
    sqrt_mass(1)=SQRT(mass(i))
enddo
close(16)
call toy_potential(q0,E0)
write(*,*) 'E0 ==> ', E0
!==============================================================================!
! 			Compute Hessian and Frequencies
!==============================================================================!
call Toy_Hessian(q0, Hess)
call Frequencies_From_Hess(Dimen,Hess,omega,U)
write(*,*) 'Test 2; Successfully Computed Hessian'
!==============================================================================!
!           Scale Eigenvectors for std normal distribution
!==============================================================================!
do i=1,Dimen
    U(:,i)=U(:,i)/sqrt_mass(:)/sqrt(omega(i))
enddo
!==============================================================================!
!               Generate Gaussian Centers with Uniform Grid
!write to file (y:=1 for convenient plotting)
!==============================================================================!
z=0d0
open(unit=17,file='centers.dat')
q(:,1)=low_bound
do ii=2,NG
    q(:,ii)=q(:,ii-1)+(-low_bound-low_bound)/NG
    write(17,*) q(:,ii), 1
enddo
close(17)
write(*,*) 'Test 3; Successfully Generated Gaussian Centers' 
!==============================================================================!
!                           Generate Alpha Scaling 
!Single paramater for each Gaussian
!==============================================================================!
alpha=alpha_par
!==============================================================================!
!                               Lambda Matrix
!==============================================================================!
do i=1,NG
    do j=i,NG
        lsum=0d0
        do k=1,dimen
            lsum=lsum+omega(k)*(q(k,i)-q(k,j))**2
        enddo
        lmat(i,j)=exp(-lsum*(alpha(i)*alpha(j)/(2.*(alpha(i)+alpha(j)))))
        lmat(j,i)=lmat(i,j)
    enddo
enddo
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
Smat=1d0
prod_omega=1d0
do k=1,Dimen
    prod_omega=prod_omega*omega(k)
enddo
do i=1,NG
    do j=i,NG
        Smat(i,j)=sqrt(2.*pi/(prod_omega*(alpha(i)+alpha(j))))
        Smat(j,i)=Smat(i,j)
    enddo
enddo
!==============================================================================!
!                   Check to see if S is positive definite
!If this is removed, you need to allocate llapack arrays before Hamiltonian 
!==============================================================================!
S1mat=Smat
S1mat=S1mat*Lmat
eigenvalues=0d0
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
CALL dsyev('n','u',NG,S1mat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'info ==>', info
open(unit=17,file='overlap_eigenvalues.dat')
do i=1,NG
    write(17,*) eigenvalues(i)
end do
close(17)
write(*,*) 'Test 4; Successfully Computed Overlap Matrix'
!==============================================================================!
!                           Kinetic Matrix (T)
!==============================================================================!
Tmat=0d0
Tpre=0d0
Tsum=0d0
do i=1,NG
    do j=i,NG
        Tpre=alpha(i)*alpha(j)*(2*pi)**(dimen/2.)/(2*(alpha(i)+alpha(j))&
            **(1+(dimen/2.)))
        Tsum=0d0
        do k=1,Dimen
            Tsum=Tsum+(omega(k)-(alpha(i)*alpha(j)*omega(k)**2&
                *(q(k,j)-q(k,i))**2/(alpha(i)+alpha(j))))
        enddo
        Tmat(i,j)=Tpre*(prod_omega)**(-0.5)*Tsum
        Tmat(j,i)=Tmat(i,j)
    enddo 
enddo
write(*,*) 'Test 5; Successfully Computed Kinetic Matrix'
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
!Scale coordinates by 1/sqrt(2) to account for std norm dist.
!==============================================================================!
z2=0d0
q2=0d0
do ii=1,Nsobol
    call sobol_stdnormal(Dimen,skip2,z2(:,ii))
    do i=1,Dimen
        q2(:,ii)=q2(:,ii)+z2(i,ii)*U(:,i)*(1/sqrt(2.))
    enddo
enddo
!==============================================================================!
!                              Evaluate Potential 
!==============================================================================!
Vmat=0d0
do i=1,NG
    do j=i,NG
        pot_ene=0d0
        do n=1,Nsobol
            call Toy_Potential(q2(:,n)&
                +((alpha(i)*q(:,i)+alpha(j)*q(:,j))/(alpha(i)+alpha(j))),pot_ene)
            Vmat(i,j)=Vmat(i,j)+pot_ene
        enddo
        Vmat(i,j)=Vmat(i,j)/sqrt(alpha(i)+alpha(j))
        Vmat(j,i)=Vmat(i,j)
    enddo
enddo
Vmat=Vmat*(2*pi)**(dimen/2.)*(prod_omega)**(-0.5)/Nsobol
write(*,*) 'Test 6; Successfully Computed Potential Matrix'
!==============================================================================!
Hmat=Vmat
Hmat=Hmat+Tmat
Hmat=Hmat*Lmat
S1mat=Smat
S1mat=S1mat*Lmat
!==============================================================================!
!                       Generalized Eigenvalue Problem
!==============================================================================!
itype=1
eigenvalues=0d0
!==============================================================================!
! Allocations needed if overlap is not diagonalized above
!==============================================================================!
!lwork = max(1,3*NG-1)                                                         !
!allocate(work(max(1,lwork)))                                                  !
!==============================================================================!
CALL dsygv(itype,'n','u',NG,Hmat,NG,S1mat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'info ==>', info
open(unit=18,file='eigenvalues.dat')
!write(*,*) 'Computed,                   True Energy,                    %Error'
do i=1,NG
    write(18,*) eigenvalues(i), ((i-1)+.5)*omega(1),&
        ((eigenvalues(i) - ((i-1)+.5)*omega(1))/(((i-1)+.5)*omega(1)))  
enddo
close(18)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
open(19,file='simulation.dat')
write(19,*) 'Natoms ==>', Natoms
write(19,*) 'Dimensionality ==>', Dimen 
write(19,*) 'N_gauss ==>', NG
write(19,*) 'N_Sobol ==>', Nsobol
write(19,*) 'Alpha Parameter ==>', alpha_par
write(19,*) 'Omega ==>', omega 
close(19)
end program DGB_ground_state
