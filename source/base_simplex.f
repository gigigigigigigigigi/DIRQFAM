c======================================================================c

      subroutine base_simplex( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      CHARACTER tb*6;
      CHARACTER tt*8;
      common /basnnn/ n0f, n0b;
      common /blokap/ nb, kb(NBX), mb(NBX), tb(NBX);
      common /bloosc/ ia(NBX,2), id(NBX,2);
      common /quaosc/ nt, nz(NTX), nr(NTX), ml(NTX),
     &                    ms(NTX), np(NTX), tt(NTX);

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);



      CHARACTER fg*1;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN base_simplex() ************************';
      endif






      N_total = 0;
      do ib = 1 , nb
          nf = id(ib,1);
          ng = id(ib,2);
          N_total = N_total + nf+ng;
      enddo

      N_blocks = 2*n0f + 3;

      if( N_blocks.ne.NBX .or. N_total.ne.NTX ) then
          stop 'Error: N_blocks =/= NBX or N_total =/= NTX!';
      endif






c-----Construction of simplex-y quantum numbers
      do ib = 1 , nb
          kap = kb(ib);
          nf  = id(ib,1);
          ng  = id(ib,2);

          i0f = ia(ib,1);
          i0g = ia(ib,2);


          il = 1;
          do i = i0f + 1 , i0f + nf
              mx = 2*( iabs(kap) - ml(i) ) - 1;
              nz_spx(il,ib) = nz(i);
              nr_spx(il,ib) = nr(i);
              ml_spx(il,ib) = ml(i) * mx;
              fg_spx(il,ib) = 'f';

              il = il + 1;
          enddo
          do i = i0g + 1 , i0g + ng
              mx = 2*( iabs(kap) - ml(i) ) - 1;
              nz_spx(il,ib) = nz(i);
              nr_spx(il,ib) = nr(i);
              ml_spx(il,ib) = ml(i) * mx;
              fg_spx(il,ib) = 'g';

              il = il + 1;
          enddo

          id_spx(ib) = il-1;
          if( id_spx(ib) .gt. NBSX ) then
              stop 'Error: NBSX (Block Size Max) too small!';
          endif

      enddo

      if( nb .ne. N_blocks ) then
          stop 'Error: nb =/= N_blocks!'
      endif






c-----Calculation of ia_spx
      ia_spx(1) = 1;
      do ib = 2 , N_blocks
          ia_spx(ib) = ia_spx(ib-1) + id_spx(ib-1);
      enddo






c-----Calculation of nf_size and ng_size
      do ib = 1 , N_blocks
          nf_size(ib) = 0;
          ng_size(ib) = 0;
          do i = 1 , id_spx(ib)
               if( fg_spx(i,ib) .eq. 'f' ) then
                   nf_size(ib) = nf_size(ib) + 1;
               else
                   ng_size(ib) = ng_size(ib) + 1;
               endif
          enddo
      enddo






c-----Sort basis with respect to nz parity
c-----This sort makes U, V matrices appear with pure real/imag blocks
c-----The code doesn't use this info at all and it is not necessary
      do ib = 1 , N_blocks

          do i = 1 , nf_size(ib)
              nz1 = nz_spx(i,ib);
              do j = i+1 , nf_size(ib)
                  nz2 = nz_spx(j,ib);
                  if( mod(nz1,2).eq.1 .and. mod(nz2,2).eq.0 ) then
                      nz_tmp       = nz_spx(i,ib);
                      nz_spx(i,ib) = nz_spx(j,ib);
                      nz_spx(j,ib) = nz_tmp;

                      nr_tmp       = nr_spx(i,ib);
                      nr_spx(i,ib) = nr_spx(j,ib);
                      nr_spx(j,ib) = nr_tmp;

                      ml_tmp       = ml_spx(i,ib);
                      ml_spx(i,ib) = ml_spx(j,ib);
                      ml_spx(j,ib) = ml_tmp;
                  endif
              enddo
          enddo

          do i = nf_size(ib)+1 , nf_size(ib)+ng_size(ib)
              nz1 = nz_spx(i,ib);
              do j = i+1 , nf_size(ib)+ng_size(ib)
                  nz2 = nz_spx(j,ib);
                  if( mod(nz1,2).eq.1 .and. mod(nz2,2).eq.0 ) then
                      nz_tmp       = nz_spx(i,ib);
                      nz_spx(i,ib) = nz_spx(j,ib);
                      nz_spx(j,ib) = nz_tmp;

                      nr_tmp       = nr_spx(i,ib);
                      nr_spx(i,ib) = nr_spx(j,ib);
                      nr_spx(j,ib) = nr_tmp;

                      ml_tmp       = ml_spx(i,ib);
                      ml_spx(i,ib) = ml_spx(j,ib);
                      ml_spx(j,ib) = ml_tmp;
                  endif
              enddo
          enddo

      enddo






c-----Test wheather the basis has correctly been constructed
      if( .true. ) then

          if( lpr ) then

              do ib = 1 , N_blocks
                  write(6,*)'-----------------------------------------';
                  do i = 1 , id_spx(ib)

                      write(6,200) i, '.  fg = ', fg_spx(i,ib),
     &                                '   nz = ', nz_spx(i,ib),
     &                                '   nr = ', nr_spx(i,ib),
     &                                '   ml = ', ml_spx(i,ib);
                  enddo
              enddo

          endif

          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  do jb = 1 , N_blocks
                      do j = 1 , id_spx(jb)
                          if( ib.eq.jb .and. i.eq.j ) CYCLE;

                          if( fg_spx(i,ib).eq.fg_spx(j,jb) .and.
     &                        nz_spx(i,ib).eq.nz_spx(j,jb) .and.
     &                        nr_spx(i,ib).eq.nr_spx(j,jb) .and.
     &                        ml_spx(i,ib).eq.ml_spx(j,jb)       ) then
                              stop 'Error: Redundancy in basis!';
                          endif

                      enddo
                  enddo
              enddo
          enddo

          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  fg  = fg_spx(i,ib);
                  nzz = nz_spx(i,ib);
                  nrr = nr_spx(i,ib);
                  mll = ml_spx(i,ib);
                  if( fg .eq. 'f' ) then
                      if( nzz + 2*nrr + abs(mll) .gt. n0f   ) then
                          stop 'Error: Basis vector out of shell!';
                      endif
                  endif
                  if( fg .eq. 'g' ) then
                      if( nzz + 2*nrr + abs(mll) .gt. n0f+1 ) then
                          stop 'Error: Basis vector out of shell!';
                      endif
                  endif
              enddo
          enddo

          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  fg  = fg_spx(i,ib);
                  nzz = nz_spx(i,ib);
                  nrr = nr_spx(i,ib);
                  mll = ml_spx(i,ib);

                  if( i .ne. index_of_vector(fg,nzz,nrr,mll,ib) ) then
                      stop 'Error: Wrong basis construction!';
                  endif

              enddo
          enddo

      endif






  200 format(i4,a,a,a,i3,a,i3,a,i3);


      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END base_simplex() **************************';
      endif

      return;
      end;
