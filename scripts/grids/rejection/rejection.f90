!==============================================================================!
!                             Rejection Method 
!==============================================================================!
!Implementation of the Rejection Method 
!See Numerical Recipies in Fortran 77; Section 7.3
!==============================================================================!
!       Modified:
!   1 September 2018
!       Author:
!   Shane Flynn 
!==============================================================================!
module rej_mod
implicit none
!==============================================================================!
!                               Global Variables 
!==============================================================================!
integer:: Dimen
!==============================================================================!
!                           Begin Module 
!==============================================================================!
contains
!==============================================================================!
subroutine pdf_eval(x,p)
!==============================================================================!
!       Discussion:
!2D potential (hard-coded)
!pdf:=(-(x-4)^2)/2 + (-(y-4)^2)/2
!==============================================================================!
implicit none
double precision::x(Dimen),p,N
double precision,parameter::pi=4*atan(1d0)
N=1
p=0d0
p=N*exp(-(((x(1)-4)**2)/2+((x(2)-4)**2)/2))
!write(*,*) 'evaluate PDF', p
end subroutine pdf_eval
!==============================================================================!
end module rej_mod
!==============================================================================!
!==============================================================================!
program reject
use rej_mod
!==============================================================================!
implicit none
integer::Nsobol,i, count_acc, count_rej
integer*8::ii,skip
double precision::accept,pxy
double precision,allocatable::z(:),x(:,:)
!==============================================================================!
!                          Read Input Data File
!==============================================================================!
read(*,*) Nsobol
skip=Nsobol
Dimen=2
allocate(z(Dimen),x(Dimen,Nsobol))
!==============================================================================!
count_acc=0
count_rej=0
z=0d0
x=0d0
open(unit=66,file='all_data.dat')
!==============================================================================!
!                     Generate Points Until total=Nsobol
!==============================================================================!
do while (count_acc < Nsobol)
    call sobol_stdnormal(Dimen,skip,z(:))
    write(66,*) z(:)
    call random_number(accept)
    call pdf_eval(z(:),pxy)
!==============================================================================!
!                    Use PDF to Accept/Reject Coordinates
!==============================================================================!
    if(pxy>accept) then
        count_acc=count_acc +1
        x(:,count_acc)=z(:)
    else
        count_rej=count_rej + 1
    endif
end do
do i=1,Nsobol
    write(*,*) x(:,i)
enddo
close(66)
open(unit=67,file='sim.dat')
write(67,*) 'Nsobol ', Nsobol
write(67,*) 'Dimensionality', dimen
write(67,*) 'Total Rejections', count_rej
close(67)
end program reject
