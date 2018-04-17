      CALL OPFILE(FNID,10,LFKEY)
!
!     read in initial background damping coef., b.c. damping coef.
!     and newtonian damping
!
      READ(10,*)
#ifdef CAA_PML
      READ(10,*) ADAMP_RE, ADAMPSB, COENT, ADAMPPML, BETA, AX, AY
#else
      READ(10,*) ADAMP_RE, ADAMPSB, COENT
#endif
!
!     read in damping coef. along every b.c.
!
      READ(10,*)
      DO IMB = 1 , MB
        READ(10,*) ADAMPBC_RE(IMB,1),ADAMPBC_RE(IMB,2),
     #             ADAMPBC_RE(IMB,3),ADAMPBC_RE(IMB,4)
      END DO
!
!     read in  newtonian cooling type damping coef. along every b.c.
!
      READ(10,*)
      DO IMB = 1 , MB
        READ(10,*) COENTBC(IMB,1),COENTBC(IMB,2),
     #             COENTBC(IMB,3),COENTBC(IMB,4)
      END DO
!
!     read in extra added damping along those lines.
!
      READ(10,*) 
      READ(10,*) NOUTBC
      READ(10,*) 
!
      ALLOCATE(IOUTBCCON(NOUTBC))
      ALLOCATE(XOUTBC(NOUTBC), YOUTBC(NOUTBC), 
     1         WTHOUTBC(NOUTBC), ADAMPOUTBC(NOUTBC), COENTOUTBC(NOUTBC))
!
      DO IOUTBC = 1 , NOUTBC
        READ(10,*) IOUTBCCON(IOUTBC),XOUTBC(IOUTBC),YOUTBC(IOUTBC),
     #             WTHOUTBC(IOUTBC),ADAMPOUTBC(IOUTBC),
     #             COENTOUTBC(IOUTBC)
      END DO
!
!     read in extra added damping in points.
!
      READ(10,*)
      READ(10,*) NSPPT 
      READ(10,*) 
!
      ALLOCATE(XSPPT(NSPPT),YSPPT(NSPPT),
     1         WTHSPPT(NSPPT),ADAMPSPPT(NSPPT),COENTSPPT(NSPPT))
!
      DO ISPPT = 1 , NSPPT
        READ(10,*) XSPPT(ISPPT),YSPPT(ISPPT),
     #             WTHSPPT(ISPPT),ADAMPSPPT(ISPPT),COENTSPPT(ISPPT)
      END DO
!
!
      CLOSE(10)
