      subroutine hpsort2 (ra, nmode, n, m,list)
c
c heapsort of multi-indices according to 
c co-lexicographical (anti-lexicographical) order
c
      implicit none

      integer :: n,m, nmode
      integer :: ra(nmode,m)
      integer :: i,ir,j,l
      integer :: rra(nmode),list(m),t
      logical :: colexsor
      
      
      list=0
      do i=1,n
       list(i)=i
      end do
      !print*, list(:)
      
      
      if (n.lt.2) return
      
      l=n/2+1
      ir=n

10    continue
        if(l.gt.1)then
        
          l=l-1
          rra(:)=ra(:,l)
          t=list(l)
          
        else
        
          rra(:)=ra(:,ir)
          t=list(ir)
          
          ra(:,ir)=ra(:,1)          
          list(ir)=list(1)   
                    
          ir=ir-1
          
          if(ir.eq.1)then
          
            ra(:,1)=rra(:)
            list(1)=t
            
            return
            
          endif
          
        endif
        
        i=l
        j=l+l
        
20      if(j.le.ir)then

          if(j.lt.ir)then
          
            !print*, ra(1,j)
            if ( colexsor (ra(1,j), ra(1,j+1), nmode) ) j = j + 1
            
          endif
          
          if ( colexsor (rra, ra(1,j), nmode) ) then  
          
            ra(:,i)=ra(:,j)
            list(i)=list(j)        
            
            i=j
            j=j+j
          else
          
            j=ir+1
          endif
          
        goto 20
        
        endif
        
        ra(:,i)=rra(:)
        list(i)=t
        
      goto 10
      
      return
      end
  
c
c
c
c
c
c
      logical function colexsor (r, p, nmode)
c
c determines if r < p according to co-lexicographical ordering
c
      implicit none
      
      integer, intent(in) :: nmode
      integer, dimension(*), intent(in) :: r, p
c      
      integer :: i
      
      colexsor = .false.
      
      

      do i = nmode, 1, -1
         if ( r(nmode-i+1) .eq. p(nmode-i+1)) cycle
         if ( r(nmode-i+1) .lt. p(nmode-i+1) ) colexsor = .true.
         return
      enddo

      return
      end
c
c
c
c
      logical function lexsor (r, p, nmode)
c
c determines if r > p according to lexicographical ordering?
c
      implicit none
      
      integer, intent(in) :: nmode
      integer, dimension(*), intent(in) :: r, p
c      
      integer :: i
     
      
      lexsor = .false.

      do i = nmode, 1, -1
         if ( r(i) .eq. p(i)) cycle
         if ( r(i) .lt. p(i) ) lexsor = .true.
         return
      enddo

      return
      end

