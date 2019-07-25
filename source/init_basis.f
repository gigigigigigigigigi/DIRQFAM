c======================================================================c

      subroutine init_basis( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /baspar/ hom, hb0, b0;
      common /defbas/ beta0, q, bp, bz;
      common /gaussh/ xh(0:NGH), wh(0:NGH), zb(0:NGH);
      common /gaussl/ xl(0:NGL), wl(0:NGL), sxl(0:NGL), rb(0:NGL);

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      common /wavefunc/ qh ( -NGH:NGH , NTX ),
     &                  qh1( -NGH:NGH , NTX ),
     &                  ql (    1:NGL , NTX ),
     &                  ql1(    1:NGL , NTX );

      common /basis/ qhql ( -NGH:NGH , 1:NGL , NTX ),
     &               qh1ql( -NGH:NGH , 1:NGL , NTX ),
     &               qhql1( -NGH:NGH , 1:NGL , NTX ),
     &               wqhql(    1:NGH , 1:NGL , NTX );

      common /k_PHI/ k_PHI, k_PHIr, k_dzPHI, k_drPHI;
      common /PHI/ PHI_U    (    NTX , KTRUNC ),
     &             PHI_SVt  ( KTRUNC , NCOORD ),
     &
     &             PHIr_U   (    NTX , KTRUNC ),
     &             PHIr_SVt ( KTRUNC , NCOORD ),
     &
     &             dzPHI_U  (    NTX , KTRUNC ),
     &             dzPHI_SVt( KTRUNC , NCOORD ),
     &
     &             drPHI_U  (    NTX , KTRUNC ),
     &             drPHI_SVt( KTRUNC , NCOORD );



      REAL*8 PHI  ( NTX , NCOORD );
      REAL*8 PHIr ( NTX , NCOORD );
      REAL*8 dzPHI( NTX , NCOORD );
      REAL*8 drPHI( NTX , NCOORD );

      parameter( THRESHOLD = 1.D-12 );
      REAL*8    A( NTX , NTX );
      REAL*8 EVal( NTX );
      REAL*8 EVec( NTX , KTRUNC );

      parameter( NWORK = 10000*NTX );
      REAL*8     WORK( NWORK );
      INTEGER*4 IWORK( 5*NTX );
      INTEGER*4 IFAIL(   NTX );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_basis() **************************';
      write(6,*) '';
      endif






c-----Calculation of qh
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nz = nz_spx(i,ib);
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  z = DBLE(isign(1,ih)) * zb(abs(ih));
                  qh(ih,i-1+ia_spx(ib)) = phi_nz(nz,z);
              enddo
          enddo
      enddo






c-----Calculation of ql
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nr = nr_spx(i,ib);
              ml = ml_spx(i,ib);
              do il = 1 , NGL
                  r = rb(il);
                  ql(il,i-1+ia_spx(ib)) = phi_nr_ml(nr,abs(ml),r);
              enddo
          enddo
      enddo






c-----Calculation of qh1
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nz = nz_spx(i,ib);
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  z = DBLE(isign(1,ih)) * zb(abs(ih));
                  qh1(ih,i-1+ia_spx(ib)) = d_phi_nz(nz,z);
              enddo
          enddo
      enddo






c-----Calculation of ql1
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nr = nr_spx(i,ib);
              ml = ml_spx(i,ib);
              do il = 1 , NGL
                  r = rb(il);
                  ql1(il,i-1+ia_spx(ib)) = d_phi_nr_ml(nr,abs(ml),r);
              enddo
          enddo
      enddo






c-----Calculation of qhql, qh1ql and qhql1
      do i = 1 , N_total
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  qhql (ih,il,i) = qh (ih,i) * ql (il,i);
                  qh1ql(ih,il,i) = qh1(ih,i) * ql (il,i);
                  qhql1(ih,il,i) = qh (ih,i) * ql1(il,i);
              enddo
          enddo
      enddo






c-----Calculation of wqhql
      do i = 1 , N_total
          do il = 1 , NGL
              do ih = 1 , NGH
                  ! Notice that in Ground State code, Gauss-Hermite
                  ! nodes are multiplied by a factor 2
                  !                                  |
                  !                                  ˇ
                  w = (0.5D0*b0*b0*b0*bp*bp*bz)*(0.5D0*wh(ih))*wl(il);
                  wqhql(ih,il,i) = DSQRT(w) * qhql(ih,il,i);
              enddo
          enddo
      enddo






c-----Calculation of PHI_U, PHI_SVt
      do i = 1 , N_total
          j = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  j = j+1;
                  PHI(i,j) = qh(ih,i)*ql(il,i);
              enddo
          enddo
      enddo

      ! A = PHI * PHI^T
      call dsyrk(     'U' ,    'N' ,
     &            N_total , NCOORD ,
     &               1.D0 ,
     &                PHI ,    NTX ,
     &               0.D0 ,
     &                 A  ,    NTX  );
      ! A = U * Sigma*Sigma^T * U^T, where PHI = U * Sigma * V^T
      ABSTOL = 2.D0*DLAMCH('S');
      IR     = N_total;
      IL     = max( IR-KTRUNC , 1 );
      call dsyevx(      'V' ,    'I' ,    'U' ,
     &              N_total ,      A ,    NTX ,
     &                 0.D0 ,   0.D0 ,
     &                   IL ,     IR ,
     &               ABSTOL ,
     &               NFOUND ,
     &                 EVal ,
     &                 EVec ,    NTX ,
     &                 WORK ,  NWORK ,   IWORK ,
     &                IFAIL ,  INFO   );
      if( INFO .ne. 0 ) then
          stop 'Error: dsyevx()!';
      endif

      if( EVal(1) .gt. THRESHOLD ) then
          stop 'Error: KTRUNC too small! Please, increase it!';
      endif

      k_PHI = 0;
      do j = NFOUND , 1 , -1
          if( EVal(j) .gt. THRESHOLD ) then
              k_PHI = k_PHI + 1;
              do i = 1 , N_total
                  PHI_U( i , k_PHI ) = EVec(i,j);
              enddo
          endif
      enddo


      ! Sigma*V^T = U^T * PHI
      call dgemm(   'T' ,    'N' ,
     &            k_PHI , NCOORD , N_total ,
     &             1.D0 ,
     &            PHI_U ,    NTX ,
     &              PHI ,    NTX ,
     &             0.D0 ,
     &          PHI_SVt , KTRUNC   );

      if( lpr ) then
          write(6,'(a,i3)') 'k_PHI = ' , k_PHI;

          acc1 = 0.D0;
          acc2 = 0.D0;
          do i = 1 , N_total
              do j = 1 , NCOORD

                  acc1 = acc1 + PHI(i,j)*PHI(i,j);

                  x = 0.D0;
                  do k = 1 , k_PHI
                      x = x + PHI_U(i,k)*PHI_SVt(k,j);
                  enddo
                  acc2 = acc2 + (x-PHI(i,j))*(x-PHI(i,j));

              enddo
          enddo

          write(6,'(a,E20.10)')
     &    '||PHI-PHI(k_PHI)||_F / ||PHI||_F = ',
     &    DSQRT(acc2) / DSQRT(acc1);

      endif






c-----Calculation of PHIr_U, PHIr_SVt
      do i = 1 , N_total
          j = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  j = j+1;
                  PHIr(i,j) = PHI(i,j) / rb(il);
              enddo
          enddo
      enddo

      ! A = PHIr * PHIr^T
      call dsyrk(     'U' ,    'N' ,
     &            N_total , NCOORD ,
     &               1.D0 ,
     &               PHIr ,    NTX ,
     &               0.D0 ,
     &                 A  ,    NTX  );
      ! A = U * Sigma*Sigma^T * U^T, where PHIr = U * Sigma * V^T
      ABSTOL = 2.D0*DLAMCH('S');
      IR     = N_total;
      IL     = max( IR-KTRUNC , 1 );
      call dsyevx(      'V' ,    'I' ,    'U' ,
     &              N_total ,      A ,    NTX ,
     &                 0.D0 ,   0.D0 ,
     &                   IL ,     IR ,
     &               ABSTOL ,
     &               NFOUND ,
     &                 EVal ,
     &                 EVec ,    NTX ,
     &                 WORK ,  NWORK ,   IWORK ,
     &                IFAIL ,  INFO   );
      if( INFO .ne. 0 ) then
          stop 'Error: dsyevx()!';
      endif

      if( EVal(1) .gt. THRESHOLD ) then
          stop 'Error: KTRUNC too small! Please, increase it!';
      endif

      k_PHIr = 0;
      do j = NFOUND , 1 , -1
          if( EVal(j) .gt. THRESHOLD ) then
              k_PHIr = k_PHIr + 1;
              do i = 1 , N_total
                  PHIr_U( i , k_PHIr ) = EVec(i,j);
              enddo
          endif
      enddo


      ! Sigma*V^T = U^T * PHIr
      call dgemm(   'T' ,    'N' ,
     &           k_PHIr , NCOORD , N_total ,
     &             1.D0 ,
     &           PHIr_U ,    NTX ,
     &             PHIr ,    NTX ,
     &             0.D0 ,
     &         PHIr_SVt , KTRUNC   );

      if( lpr ) then
          write(6,'(a,i3)') 'k_PHIr = ' , k_PHIr;

          acc1 = 0.D0;
          acc2 = 0.D0;
          do i = 1 , N_total
              do j = 1 , NCOORD

                  acc1 = acc1 + PHIr(i,j)*PHIr(i,j);

                  x = 0.D0;
                  do k = 1 , k_PHIr
                      x = x + PHIr_U(i,k)*PHIr_SVt(k,j);
                  enddo
                  acc2 = acc2 + (x-PHIr(i,j))*(x-PHIr(i,j));

              enddo
          enddo

          write(6,'(a,E20.10)')
     &    '||PHIr-PHIr(k_PHIr)||_F / ||PHIr||_F = ',
     &    DSQRT(acc2) / DSQRT(acc1);

      endif






c-----Calculation of dzPHI_U, dzPHI_SVt
      do i = 1 , N_total
          j = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  j = j+1;
                  dzPHI(i,j) = qh1(ih,i)*ql(il,i);
              enddo
          enddo
      enddo

      ! A = dzPHI * dzPHI^T
      call dsyrk(     'U' ,    'N' ,
     &            N_total , NCOORD ,
     &               1.D0 ,
     &              dzPHI ,    NTX ,
     &               0.D0 ,
     &                 A  ,    NTX  );
      ! A = U * Sigma*Sigma^T * U^T, where dzPHI = U * Sigma * V^T
      ABSTOL = 2.D0*DLAMCH('S');
      IR     = N_total;
      IL     = max( IR-KTRUNC , 1 );
      call dsyevx(      'V' ,    'I' ,    'U' ,
     &              N_total ,      A ,    NTX ,
     &                 0.D0 ,   0.D0 ,
     &                   IL ,     IR ,
     &               ABSTOL ,
     &               NFOUND ,
     &                 EVal ,
     &                 EVec ,    NTX ,
     &                 WORK ,  NWORK ,   IWORK ,
     &                IFAIL ,  INFO   );
      if( INFO .ne. 0 ) then
          stop 'Error: dsyevx()!';
      endif

      if( EVal(1) .gt. THRESHOLD ) then
          stop 'Error: KTRUNC too small! Please, increase it!';
      endif

      k_dzPHI = 0;
      do j = NFOUND , 1 , -1
          if( EVal(j) .gt. THRESHOLD ) then
              k_dzPHI = k_dzPHI + 1;
              do i = 1 , N_total
                  dzPHI_U( i , k_dzPHI ) = EVec(i,j);
              enddo
          endif
      enddo


      ! Sigma*V^T = U^T * dzPHI
      call dgemm(   'T' ,    'N' ,
     &          k_dzPHI , NCOORD , N_total ,
     &             1.D0 ,
     &          dzPHI_U ,    NTX ,
     &            dzPHI ,    NTX ,
     &             0.D0 ,
     &        dzPHI_SVt , KTRUNC   );

      if( lpr ) then
          write(6,'(a,i3)') 'k_dzPHI = ' , k_dzPHI;

          acc1 = 0.D0;
          acc2 = 0.D0;
          do i = 1 , N_total
              do j = 1 , NCOORD

                  acc1 = acc1 + dzPHI(i,j)*dzPHI(i,j);

                  x = 0.D0;
                  do k = 1 , k_dzPHI
                      x = x + dzPHI_U(i,k)*dzPHI_SVt(k,j);
                  enddo
                  acc2 = acc2 + (x-dzPHI(i,j))*(x-dzPHI(i,j));

              enddo
          enddo

          write(6,'(a,E20.10)')
     &    '||dzPHI-dzPHI(k_dzPHI)||_F / ||dzPHI||_F = ',
     &    DSQRT(acc2) / DSQRT(acc1);

      endif






c-----Calculation of drPHI_U, drPHI_SVt
      do i = 1 , N_total
          j = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  j = j+1;
                  drPHI(i,j) = qh(ih,i)*ql1(il,i);
              enddo
          enddo
      enddo

      ! A = drPHI * drPHI^T
      call dsyrk(     'U' ,    'N' ,
     &            N_total , NCOORD ,
     &               1.D0 ,
     &              drPHI ,    NTX ,
     &               0.D0 ,
     &                 A  ,    NTX  );
      ! A = U * Sigma*Sigma^T * U^T, where drPHI = U * Sigma * V^T
      ABSTOL = 2.D0*DLAMCH('S');
      IR     = N_total;
      IL     = max( IR-KTRUNC , 1 );
      call dsyevx(      'V' ,    'I' ,    'U' ,
     &              N_total ,      A ,    NTX ,
     &                 0.D0 ,   0.D0 ,
     &                   IL ,     IR ,
     &               ABSTOL ,
     &               NFOUND ,
     &                 EVal ,
     &                 EVec ,    NTX ,
     &                 WORK ,  NWORK ,   IWORK ,
     &                IFAIL ,  INFO   );
      if( INFO .ne. 0 ) then
          stop 'Error: dsyevx()!';
      endif

      if( EVal(1) .gt. THRESHOLD ) then
          stop 'Error: KTRUNC too small! Please, increase it!';
      endif

      k_drPHI = 0;
      do j = NFOUND , 1 , -1
          if( EVal(j) .gt. THRESHOLD ) then
              k_drPHI = k_drPHI + 1;
              do i = 1 , N_total
                  drPHI_U( i , k_drPHI ) = EVec(i,j);
              enddo
          endif
      enddo


      ! Sigma*V^T = U^T * drPHI
      call dgemm(   'T' ,    'N' ,
     &          k_drPHI , NCOORD , N_total ,
     &             1.D0 ,
     &          drPHI_U ,    NTX ,
     &            drPHI ,    NTX ,
     &             0.D0 ,
     &        drPHI_SVt , KTRUNC   );

      if( lpr ) then
          write(6,'(a,i3)') 'k_drPHI = ' , k_drPHI;

          acc1 = 0.D0;
          acc2 = 0.D0;
          do i = 1 , N_total
              do j = 1 , NCOORD

                  acc1 = acc1 + drPHI(i,j)*drPHI(i,j);

                  x = 0.D0;
                  do k = 1 , k_drPHI
                      x = x + drPHI_U(i,k)*drPHI_SVt(k,j);
                  enddo
                  acc2 = acc2 + (x-drPHI(i,j))*(x-drPHI(i,j));

              enddo
          enddo

          write(6,'(a,E20.10)')
     &    '||drPHI-drPHI(k_drPHI)||_F / ||drPHI||_F = ',
     &    DSQRT(acc2) / DSQRT(acc1);

      endif






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_basis() ****************************';
      write(6,*) '';
      endif

      return;
      end;
