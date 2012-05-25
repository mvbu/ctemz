
      integer fixedRandFlag ! if non-zero, randproto() will "generate" numbers from a file
      integer fixedRandFileOpened, fixedRandCounter
      dimension fixedRandData(100000)
      real fixedRandData
      common/cfixedrand/fixedRandFlag, fixedRandFileOpened, fixedRandData, fixedRandCounter
      fixedRandFileOpened = 0
      fixedRandFlag = 1
      fixedRandCounter = 0
      do i=1,100020
         x = randproto(0)
         print *, x
      end do
      end


      ! Function to call rand() most of the time, but when fixedRandFlag is non-zero, 
      ! get the random values from a file. When all the numbers are used up, start
      ! back at the beginning of the number list
      function randproto(seed)
      common/cfixedrand/fixedRandFlag, fixedRandFileOpened, fixedRandData, fixedRandCounter
      dimension fixedRandData(100000)
      real fixedRandData
      integer fixedRandFlag, fixedRandFileOpened, fixedRandCounter
      integer*4 seed

      if(fixedRandFileOpened.eq.0) then
         ! need to open the file with all the pre-generated random nums
         ! and load into a shared array
         ! At some point, I want to have the ability to set an environment variable for the filename:
         !character(len=512) :: filename
         !call get_environment_variable('BLZFIXEDRAND', filename, 512, status)
         !if(status.eq.1) then
         !   filename = './fixedrand.dat'
         !endif
         open(99, FILE = 'fixedrand.dat') ! At the moment, as you can see, this file must be in running/current directory
         do i=1,100000
            read(99, *) fixedRandData(i)
         end do
         fixedRandFileOpened = 1
         fixedRandCounter = 1
         close(99)
      endif

      if(fixedRandFlag.eq.0) then
         randproto = rand(seed)
      else
         randproto = fixedRandData(fixedRandCounter)
         if(fixedRandCounter.ge.100000) then
            fixedRandCounter = 1
         else
            fixedRandCounter = fixedRandCounter + 1
         endif
      endif
      return
      end
