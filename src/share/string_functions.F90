!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! MODULE: string_functions --- 
!
! !DESCRIPTION:
!   String functions to manupulate with arrays of strings.
!
! COPYING
!
!   Copyright (C) 2023 BFM System Team (bfm_st@cmcc.it)
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation.
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! INCLUDE
#include"cppdefs.h"
!
! INTERFACE
   MODULE string_functions
!
! USES
   IMPLICIT NONE

   public empty, index_trim, getseq_number, replace_char

   contains

     logical function empty(string)
!
! DESCRIPTION:
!   to test on empty string

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     character(len=*)         ::string
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     empty=lle(string,' ')
     return

     end function empty

!-----------------------------------------------------------------------

     integer function index_trim(string,text)
! DESCRIPTION
!   Trim string

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     character(len=*)             :: string*(*),text*(*)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     if (.not.empty(text)) then
        index_trim=index(string,text(1:len_trim(text)))
     else
        index_trim=0
     endif
     return

     end function index_trim

!-----------------------------------------------------------------------

     integer function getseq_number(string,name,n,exact,start)
!
! DESCRIPTION
!
     implicit none

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     integer,intent(IN)                       :: n
     character(len=*),dimension(n),intent(IN) :: name
     character(len=*),intent(IN)              :: string
     logical,intent(IN)                       :: exact
     integer,intent(IN),optional              :: start

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     integer                 :: i,j,k
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     k=1
     if ( present(start) ) k=start
     do i=k,n
       j=index_trim(name(i),string)
       if (j == 1) then
         if (exact ) THEN
           if (name(i) /= string) j=0
         endif
         if (j == 1) then
           getseq_number=i
           return
         endif
       endif
     enddo
     getseq_number=0
     return
     end function getseq_number

!-----------------------------------------------------------------------

     subroutine replace_char(str,tar,rep)
!
! DESCRIPTION:
!   replace character
!
    implicit none

  ! INPUT/OUTPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    character(LEN=*), intent(INOUT) :: str
    character(LEN=*), intent(IN)    :: tar, rep

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     integer                 :: times
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     times = scan(str, tar)
     do while ( times .ne. 0 )
        str(times:times) = rep
        times = scan(str, tar)
     end do

   end subroutine replace_char

end module string_functions

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
