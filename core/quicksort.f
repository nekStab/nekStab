c-----------------------------------------------------------------------
      subroutine quicksort2(n, arr, idx)
      implicit none
      integer n, m, nstack
      real arr(n)
      integer idx(n)
      parameter (m=7, nstack=50)
      integer i, ir, j, jstack, k, l, istack(nstack)
      real a, temp
      integer b, temp_int

      jstack = 0
      l = 1
      ir = n

 1    if (ir-l .lt. m) then
         do 12 j = l+1, ir
            a = arr(j)
            b = idx(j)
            do 11 i = j-1, l, -1
               if (arr(i) .le. a) goto 2
               arr(i+1) = arr(i)
               idx(i+1) = idx(i)
 11         continue
            i = l-1
 2          arr(i+1) = a
            idx(i+1) = b
 12      continue

         if (jstack .eq. 0) return

         ir = istack(jstack)
         l = istack(jstack - 1)
         jstack = jstack - 2
      else
         k = (l + ir) / 2
         temp = arr(k)
         arr(k) = arr(l+1)
         arr(l+1) = temp
         temp_int = idx(k)
         idx(k) = idx(l+1)
         idx(l+1) = temp_int

         if (arr(l) .gt. arr(ir)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp

            temp_int = idx(l)
            idx(l) = idx(ir)
            idx(ir) = temp_int
         endif

         if (arr(l+1) .gt. arr(ir)) then
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp

            temp_int = idx(l+1)
            idx(l+1) = idx(ir)
            idx(ir) = temp_int
         endif

         if (arr(l) .gt. arr(l+1)) then
            temp = arr(l)
            arr(l) = arr(l+1)
            arr(l+1) = temp
            temp_int = idx(l)
            idx(l) = idx(l+1)
            idx(l+1) = temp_int
         endif

         i = l + 1
         j = ir
         a = arr(l+1)
         b = idx(l+1)

 3       continue

         i = i+1

         if (arr(i) .lt. a) goto 3

 4       continue

         j = j -1

         if (arr(j) .gt. a) goto 4
         if (j .lt. i) goto 5

         temp = arr(i)
         arr(i) = arr(j)
         arr(j) = temp

         temp_int = idx(i)
         idx(i) = idx(j)
         idx(j) = temp_int
         goto 3

 5       arr(l+1) = arr(j)
         arr(j) = a
         idx(l+1) = idx(j)
         idx(j) = b
         jstack = jstack + 2
         if (jstack .gt. nstack) pause "..." !'NSTACK too small in quicksort2'

         if (ir-i+1 .ge. j-1) then
            istack(jstack) = ir
            istack(jstack-1) = i
            ir = j-1
         else
            istack(jstack) = j-1
            istack(jstack-1) = l
            l = i
         endif
      endif
      goto 1
      end subroutine quicksort2
c-----------------------------------------------------------------------