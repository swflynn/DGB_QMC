!=============================================================================80
!                 Distributed Gaussian Basis (Ground State Energy)
!==============================================================================!
!       Discussion:
!DGB analysis for 2D seperable potential (hard-coded) 
!Normal mode coordinates for analysis, Integration with quasi-random grid
!==============================================================================!
!       Modified:
!   7 January 2018
!       Author:
!   Shane Flynn 
!==============================================================================!
module dgb_groundstate
implicit none
!==============================================================================!
!                            Global Paramaters
!==============================================================================!
double precision,parameter::Hmass=1d0
double precision,parameter::c1=1d0
double precision,parameter::c2=2d0
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!       Discussion:
!Natoms             ==> Number of Atoms 
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
!Hard-coded Potential Energy 
!V:=c1*(x+y)^2+c2*(x-y)^2
!==============================================================================!
implicit none
double precision::x(Dimen),energies
energies=c1*(x(1)+x(2))**2+c2*(x(1)-x(2))**2
!write(*,*) 'Energy from Toy_Potential Subroutine ==> ', energies
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
forces(1)=-(2.*x(1)*(c1+c2) + 2.*(x(2)*(c1-c2)))
forces(2)=-(2.*x(1)*(c1-c2) + 2.*(x(2)*(c1+c2)))
!write(*,*) 'Forces from Toy_Force Subroutine ==> ', forces
end subroutine Toy_Force
!==============================================================================!
subroutine Toy_Hessian(x,Hess_Mat)
!==============================================================================!
!       Discussion:
!Numerically computed Hessian using forces from "Toy_Force Subroutine"
!       Variables:
!s          ==> Perturbation Parameter 
!Hess_Mat   ==> (Dimen,Dimen); Symmetrized Mass-Scaled Hessian
!x          ==> (Dimen); XYZ Configuration 
!==============================================================================!
implicit none 
integer::i,j
double precision::Hess_Mat(Dimen,Dimen),x(Dimen),r(Dimen),force(Dimen)
double precision::force0(Dimen)
double precision,parameter::s=1d-6
r=x
call Toy_Force(r,force0)
do i=1,Dimen
    r(i)=x(i)+s
    call Toy_Force(r, force)
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
!write(*,*) 'Hessian from Toy_Hessian Subroutine ==> ', Hess_Mat
end subroutine Toy_Hessian
!==============================================================================!
subroutine Freq_Hess(Dimen,Hess,omega,U)
!==============================================================================!
!       Discussion:
!Compute Eigenvalues and Eigenvectors of Hessian
!Uses the LLAPACK real symmetric eigen-solver (dsygev)
!       Variables:
!Hess   ==> (Dimen,Dimen); Hessian Matrix
!omega  ==> (Dimen); Hessian Eigenvalues
!U      ==> (Dimen,Dimen); Hessian Eigenvectors
!       LLAPACK (dsyev):
!v      ==> Compute both Eigenvalues and Eigenvectors
!u      ==> Use Upper-Triangle of matrix
!==============================================================================!
implicit none
integer::i,info,lwork,Dimen
double precision::Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen)
double precision,allocatable::work(:) 
lwork=max(1,3*Dimen-1)          !suggested by LAPACK Developers
allocate(work(max(1,lwork)))    !suggested by LAPACK Developers 
U=Hess  
call dsyev('v','u',Dimen,U,Dimen,omega,work,lwork,info) 
write(*,*) 'Frequencies from the Hessian:'
do i=Dimen,1,-1
    omega(i)=sign(sqrt(abs(omega(i))),omega(i))
    write(*,*) omega(i), 'normalized = 1?', sum(U(:,i)**2)
enddo
end subroutine Freq_Hess
!==============================================================================!
end module dgb_groundstate
!==============================================================================!
!==============================================================================!
program DGB_2D
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
!r              ==> (Dimen,NG): Gaussian Centers 
!r2             ==> (Dimen): ith Gaussian coordinate for integration
!force          ==> (Dimen): Forces. See "Toy_Force" 
!Hess           ==> (Dimen,Dimen): Mass-Scaled Hessian. See "Toy_Hessian"
!omega          ==> (Dimen): Frequencies (Eigenvalues). See "Freq_Hess"
!U              ==> (Dimen,Dimen): Hess Eigenvectors. See "Freq_Hess"
!z              ==> (Dimen,NG): Quasi-Random Sequence for Gaussian Centers
!z2             ==> (Dimen,Nsobol): Quasi-Random Sequence for Integration 
!Smat           ==> (NG,NG): Overlap Matrix
!kinetic        ==> Kinetic Energy Matrix (computed on the fly for i,j element)
!Vmat           ==> (NG,NG): Potential Energy Matrix
!Hmat           ==> (NG,NG): Hamiltonian Matrix (V+T) for eigenvalue problem
!eigenvalues    ==> (NG): Eigenvalues of the Hamiltonian Matrix
!lambda         ==> (NG): Eigenvalues of the Overlap Matrix
!r_ij           ==> center of the i,j matrix element 
!E0             ==> Energy Evaluation at the minimum configuration
!pot_ene        ==> Potential Energy evaluated in q space
!==============================================================================!
implicit none
character(len=50)::coord_in
integer::NG,Nsobol
integer::i,j,k,l,m,n,counter,data_freq
integer*8::ii,skip,skip2
double precision::E0,pot_ene,s_sum,kinetic,alpha0
!==============================================================================!
double precision,allocatable::q0(:),q1(:),force(:),r(:,:),r2(:),eigenvalues(:)
double precision,allocatable::Hess(:,:),omega(:),U(:,:),z(:,:),z2(:,:),Hmat(:,:)
double precision,allocatable::Smat(:,:),alpha(:),Vmat(:,:),lambda(:),r_ij(:)
!==============================================================================!
!                       LLAPACK dsygv variables                                !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable::work(:)
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) coord_in
read(*,*) alpha0
read(*,*) NG
read(*,*) Nsobol
read(*,*) data_freq
skip=NG
skip2=Nsobol
write(*,*) 'Test 1; Successfully Read Input Data File!'
!==============================================================================!
!                         Set Input Water Geometry 
!==============================================================================!
open(16,File=coord_in)
read(16,*) Natoms
read(16,*) 
Dimen=2*Natoms 
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(Dimen),q0(Dimen),q1(Dimen))
allocate(force(Dimen),Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen),z(Dimen,NG))
allocate(z2(Dimen,Nsobol),r(Dimen,NG),alpha(NG),Smat(NG,NG),eigenvalues(NG))
allocate(Vmat(NG,NG),Hmat(NG,NG),r2(Dimen),lambda(NG),r_ij(Dimen))
!==============================================================================!
!                         Input Configuration Energy
!==============================================================================!
do i=1,Natoms
    read(16,*) atom_type(i),q0(2*i-1:2*i)   
    mass(i)=Atom_mass(atom_type(i))
    sqrt_mass(2*i-1:2*i)=sqrt(mass(i))
enddo
close(16)
write(*,*) 'Dimen ==> ', Dimen
call toy_potential(q0,E0)
write(*,*) 'E0 ==> ', E0
!==============================================================================!
! 			Compute Hessian and Frequencies
!==============================================================================!
call Toy_Hessian(q0,Hess)
call Freq_Hess(Dimen,Hess,omega,U)
write(*,*) 'Test 2; Successfully Computed Hessian'
!==============================================================================!
!                 Generate Gaussian Centers with a quasi grid
! Assume centers are in normal-mode space r (not coordinate space q)
!==============================================================================!
r=0d0
do ii=1,NG
    call sobol_stdnormal(Dimen,skip,z(:,ii))
    r(:,ii)=z(:,ii)/sqrt(abs(omega(:)))
enddo
write(*,*) 'Test 3; Successfully Generated Gaussian Centers'
!==============================================================================!
!                       Generate Alpha Scaling 
!==============================================================================!
do i=1,NG
    alpha(i)=alpha0*NG**(2./dimen)*exp(-1./dimen*sum(r(:,i)**2*omega(:)))
enddo
open(unit=17,file='centers.dat')
open(unit=18,file='width.dat')
do i=1,NG
    write(17,*) alpha(i), r(1,i), r(2,i)
    write(18,*) 1./alpha(i), r(1,i), r(2,i)
enddo
close(17)
close(18)
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
do i=1,NG
    do j=i,NG
        s_sum=sum(omega(:)*(r(:,i)-r(:,j))**2)
        Smat(i,j)=sqrt(alpha(i)+alpha(j))**(-dimen)&           
                 *exp(-0.5*alpha(i)*alpha(j)/(alpha(i)+alpha(j))*s_sum)
        Smat(j,i)=Smat(i,j)
    enddo
enddo
!==============================================================================!
!                   Check to see if S is positive definite
! If this is removed, you need to allocate llapack arrays before Hamiltonian 
!==============================================================================!
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
call dsyev('v','u',NG,Smat,NG,lambda,work,Lwork,info)
write(*,*) 'Info (Overlap Matrix) ==> ', info
open(unit=19,file='overlap_eigenvalues.dat')
do i=1,NG
    write(19,*) lambda(i)
enddo
close(19)
write(*,*) 'Test 4; Successfully Computed Overlap Matrix'
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
!==============================================================================!
do ii=1,Nsobol
    call sobol_stdnormal(Dimen,skip2,z2(:,ii))
enddo
write(*,*) 'Test 5; Successfully Generated Integration Sequence'
!==============================================================================!
!                           Evaluate Potential 
!==============================================================================!
Vmat=0d0
counter=0
open(unit=20,file='eigenvalues.dat')
do while(counter.lt.Nsobol/data_freq)
    do i=1,NG
        do j=i,NG
            r_ij(:)=(alpha(i)*r(:,i)+alpha(j)*r(:,j))/(alpha(i)+alpha(j))
            do l=1+(counter*data_freq),(counter+1)*data_freq
                r2(:)=r_ij(:)+z2(:,l)/sqrt(omega(:)*(alpha(i)+alpha(j)))
!==============================================================================!
!                           Transformation to Cartesian
!==============================================================================!
                q1=q0
                q1(1:Dimen)=q1(1:Dimen)+matmul(U,r2)/sqrt_mass(:)
                call Toy_Potential(q1,pot_ene)
                Vmat(i,j)=Vmat(i,j)+pot_ene
            enddo
        enddo
    enddo
!==============================================================================!
!                   Solve Generalized Eigenvalue Problem
!==============================================================================!
    write(*,*) 'Iteration = ', (counter+1)*data_freq
    Hmat=Vmat
    do m=1,NG
        do n=m,NG
            s_sum=sum(omega(:)*(r(:,m)-r(:,n))**2)
            Smat(m,n)=sqrt(alpha(m)+alpha(n))**(-dimen)&           
                     *exp(-0.5*alpha(m)*alpha(n)/(alpha(m)+alpha(n))*s_sum)
            Smat(n,m)=Smat(m,n)
            kinetic=Smat(m,n)*0.5*alpha(m)*alpha(n)/(alpha(m)+alpha(n))&
            *sum(omega(:)-(alpha(m)*alpha(n)*(omega(:)**2*(r(:,m)-r(:,n))**2)&
            /(alpha(m)+alpha(n))))
            Hmat(m,n)=Hmat(m,n)*(Smat(m,n)/((counter+1)*data_freq))+kinetic
            Hmat(n,m)=Hmat(m,n)
        enddo
    enddo
    itype=1
    eigenvalues=0d0
    call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
    write(20,*) eigenvalues(:)
    write(*,*) 'info ==> ', info
    counter=counter + 1
enddo
close(20)
!==============================================================================!
!                           Theorectical Eigenvalues                           !
!==============================================================================!
open(unit=21,file='theory.dat')
do i=0,NG
    do j=0,NG
        write(21,*) (0.5+i)*omega(1)+(0.5+j)*omega(2)
    enddo
enddo
close(21)
write(*,*) 'Hello Universe!'
!==============================================================================!
!                               Output File                                    !
!==============================================================================!
open(90,file='simulation.dat')
write(90,*) 'Natoms ==>', Natoms
write(90,*) 'Dimensionality Input System ==>', Dimen 
write(90,*) 'N_gauss ==>', NG
write(90,*) 'N_Sobol ==>', Nsobol
write(90,*) 'omega ==>', omega(1), omega(2)
write(90,*) 'E0 ==>', E0
close(90)
end program DGB_2D
