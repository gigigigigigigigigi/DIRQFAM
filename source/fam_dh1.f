c======================================================================c

      subroutine fam_dh1( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      common /basis/ qhql ( -NGH:NGH , 1:NGL , NTX ),
     &               qh1ql( -NGH:NGH , 1:NGL , NTX ),
     &               qhql1( -NGH:NGH , 1:NGL , NTX ),
     &               wqhql(    1:NGH , 1:NGL , NTX );

      COMPLEX*16 dVpS, dVmS, dSig_z, dSig_r, dSig_p;
      common /fam_pot/ dVpS  ( -NGH:NGH , 1:NGL , 2 ),
     &                 dVmS  ( -NGH:NGH , 1:NGL , 2 ),
     &                 dSig_z( -NGH:NGH , 1:NGL , 2 ),
     &                 dSig_r( -NGH:NGH , 1:NGL , 2 ),
     &                 dSig_p( -NGH:NGH , 1:NGL , 2 );

      COMPLEX*16 dh_1, dh_2;
      common /delta_h/ dh_1( NTX , NTX , 2 ),
     &                 dh_2( NTX , NTX , 2 );



      CHARACTER fg1, fg2;
      COMPLEX*16 z;
      COMPLEX*16 dVpS_odd  ( 1:NGH , 1:NGL , 2 ),
     &           dVpS_evn  ( 1:NGH , 1:NGL , 2 ),
     &           dVmS_odd  ( 1:NGH , 1:NGL , 2 ),
     &           dVmS_evn  ( 1:NGH , 1:NGL , 2 ),
     &           dSig_z_odd( 1:NGH , 1:NGL , 2 ),
     &           dSig_z_evn( 1:NGH , 1:NGL , 2 ),
     &           dSig_r_odd( 1:NGH , 1:NGL , 2 ),
     &           dSig_r_evn( 1:NGH , 1:NGL , 2 ),
     &           dSig_p_odd( 1:NGH , 1:NGL , 2 ),
     &           dSig_p_evn( 1:NGH , 1:NGL , 2 );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_dh1() ******************************';
      write(6,*) '';
      endif






c-----Initializing odd/even parts of induced potentials
      do it = 1 , 2
        do il = 1 , NGL
          do ih = 1 , NGH

            dVpS_odd(ih,il,it)   = dVpS(ih,il,it)   - dVps(-ih,il,it);
            dVpS_evn(ih,il,it)   = dVpS(ih,il,it)   + dVps(-ih,il,it);

            dVmS_odd(ih,il,it)   = dVmS(ih,il,it)   - dVms(-ih,il,it);
            dVmS_evn(ih,il,it)   = dVmS(ih,il,it)   + dVms(-ih,il,it);

            dSig_z_odd(ih,il,it) = dSig_z(ih,il,it) - dSig_z(-ih,il,it);
            dSig_z_evn(ih,il,it) = dSig_z(ih,il,it) + dSig_z(-ih,il,it);

            dSig_r_odd(ih,il,it) = dSig_r(ih,il,it) - dSig_r(-ih,il,it);
            dSig_r_evn(ih,il,it) = dSig_r(ih,il,it) + dSig_r(-ih,il,it);

            dSig_p_odd(ih,il,it) = dSig_p(ih,il,it) - dSig_p(-ih,il,it);
            dSig_p_evn(ih,il,it) = dSig_p(ih,il,it) + dSig_p(-ih,il,it);

          enddo
        enddo
      enddo






c-----Calculation of dh_1 matrix
      dh_1 = COMPLEX( 0.D0 , 0.D0 );
      K    = K_multipole;
      do it = 1 , 2
        do ib2 = 1 , N_blocks
          do ib1 = 1 , ib2

            j0 = ia_spx(ib2);
            j1 = ia_spx(ib2)+id_spx(ib2)-1;
            do j = j0 , j1
              fg2 = fg_spx(j-j0+1,ib2);
              nz2 = nz_spx(j-j0+1,ib2);
              ml2 = ml_spx(j-j0+1,ib2);


              i0 = ia_spx(ib1);
              i1 = ia_spx(ib1)+id_spx(ib1)-1;
              if( ib1 .eq. ib2 ) i1 = j;
              do i = i0 , i1
                fg1 = fg_spx(i-i0+1,ib1);
                nz1 = nz_spx(i-i0+1,ib1);
                ml1 = ml_spx(i-i0+1,ib1);




                if( fg1.eq.'f' .and. fg2.eq.'f' ) then

                    if( abs(ml1-ml2) .ne. K ) CYCLE;

                    z = COMPLEX( 0.D0 , 0.D0 );
                    if( MOD(nz1+nz2,2) .eq. 0 ) then
                        do il = 1 , NGL
                            do ih = 1 , NGH
                                z = z + dVpS_evn(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                            enddo
                        enddo
                    else
                        do il = 1 , NGL
                            do ih = 1 , NGH
                                z = z + dVpS_odd(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                            enddo
                        enddo
                    endif

                    if( K .ne. 0 ) then
                        z = 0.5D0 * z;
                    endif

                    dh_1(i,j,it) = z;
                    dh_1(j,i,it) = z;

                endif


                if( fg1.eq.'g' .and. fg2.eq.'g' ) then

                    if( abs(ml1-ml2) .ne. K ) CYCLE;

                    z = COMPLEX( 0.D0 , 0.D0 );
                    if( MOD(nz1+nz2,2) .eq. 0 ) then
                        do il = 1 , NGL
                            do ih = 1 , NGH
                                z = z + dVmS_evn(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                            enddo
                        enddo
                    else
                        do il = 1 , NGL
                            do ih = 1 , NGH
                                z = z + dVmS_odd(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                            enddo
                        enddo
                    endif

                    if( K .ne. 0 ) then
                        z = 0.5D0 * z;
                    endif

                    dh_1(i,j,it) = z;
                    dh_1(j,i,it) = z;

                endif


                if( fg1.eq.'g' .and. fg2.eq.'f' ) then

                    if( abs(ml1+ml2+1) .eq. K ) then

                        z = COMPLEX( 0.D0 , 0.D0 );

                        if( K .eq. 0 ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_r_evn(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_r_odd(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif

                        if( K.ne.0 .and. (ml1+ml2+1).eq.(+K) ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_evn(ih,il,it)
     &                                    - dSig_p_evn(ih,il,it) )
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_odd(ih,il,it)
     &                                    - dSig_p_odd(ih,il,it) )
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif

                        if( K.ne.0 .and. (ml1+ml2+1).eq.(-K) ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_evn(ih,il,it)
     &                                    + dSig_p_evn(ih,il,it) )
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_odd(ih,il,it)
     &                                    + dSig_p_odd(ih,il,it) )
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif


                        z = z * COMPLEX( 0.D0 , 1.D0 );
                        if( K .ne. 0 ) then
                            z = 0.5D0 * z;
                        endif

                        dh_1(i,j,it) = + z;
                        dh_1(j,i,it) = - z;

                    endif

                    if( abs(ml1-ml2) .eq. K ) then

                        z = COMPLEX( 0.D0 , 0.D0 );
                        if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_z_evn(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                        else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_z_odd(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                        endif

                        if( K .ne. 0 ) then
                            z = 0.5D0 * z;
                        endif

                        dh_1(i,j,it) = - z;
                        dh_1(j,i,it) = - z;

                    endif

                endif


                if( fg1.eq.'f' .and. fg2.eq.'g' ) then

                    if( abs(ml1+ml2+1) .eq. K ) then

                        z = COMPLEX( 0.D0 , 0.D0 );

                        if( K .eq. 0 ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_r_evn(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_r_odd(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif

                        if( K.ne.0 .and. (ml1+ml2+1).eq.(+K) ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_evn(ih,il,it)
     &                                    - dSig_p_evn(ih,il,it) )
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_odd(ih,il,it)
     &                                    - dSig_p_odd(ih,il,it) )
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif

                        if( K.ne.0 .and. (ml1+ml2+1).eq.(-K) ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_evn(ih,il,it)
     &                                    + dSig_p_evn(ih,il,it) )
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_odd(ih,il,it)
     &                                    + dSig_p_odd(ih,il,it) )
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif


                        z = z * COMPLEX( 0.D0 , 1.D0 );
                        if( K .ne. 0 ) then
                            z = 0.5D0 * z;
                        endif

                        dh_1(i,j,it) = - z;
                        dh_1(j,i,it) = + z;

                    endif

                    if( abs(ml1-ml2) .eq. K ) then

                        z = COMPLEX( 0.D0 , 0.D0 );
                        if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_z_evn(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                        else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_z_odd(ih,il,it)
     &                                * wqhql(ih,il,i)*wqhql(ih,il,j);
                              enddo
                            enddo
                        endif

                        if( K .ne. 0 ) then
                            z = 0.5D0 * z;
                        endif

                        dh_1(i,j,it) = - z;
                        dh_1(j,i,it) = - z;

                    endif

                endif


              enddo
            enddo

          enddo
        enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_dh1() ********************************';
      write(6,*) '';
      endif

      return;
      end;
