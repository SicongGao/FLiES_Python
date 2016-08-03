
      character*20 dum1,dum2,dum3,dum4,dum5,dum6,dum7

      open(12,file="gas_TR_lo")
      open(13,file="gas_TR_lo_conv")

      do i=1,53
         read(12,*) dum1
         read(12,*) num,w1,w2
         read(12,*) dum1
         write(13,*) dum1
         write(13,1000) num,w1,w2
         write(13,*) dum1
         do j=1,12
            read(12,*) iz,ext,abs1
            write(13,*) iz,ext,abs1,"  0.0  0.0"
         end do
      end do

      close(12)
      close(13)
 1000 format(I4,2(F6.2))
      stop 
      end
