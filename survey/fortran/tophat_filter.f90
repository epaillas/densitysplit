program tophat_filter
  use procedures
  use OMP_LIB
  implicit none
  
  real*8 :: rgrid, rfilter, gridmax, gridmin
  real*8 :: disx, disy, disz, dis2
  real*8 :: dim1_min, dim1_max, dim1_min2, dim1_max2
  
  integer*8 :: ndata1, ndata2, nrandom1, nrandom2
  integer*8 :: i, ii, ix, iy, iz
  integer*8 :: nrows, ncols
  integer*8 :: ipx, ipy, ipz, ndif
  integer*8 :: ngrid
  integer*4 :: nthreads, use_weights
  integer*8, dimension(:, :, :), allocatable :: lirst_data2, lirst_random2

  integer*8, dimension(:), allocatable :: ll_data2, ll_random2
  
  real*8, allocatable, dimension(:,:)  :: data1, data2, random1, random2
  real*8, dimension(:), allocatable :: D1D2, D1R2, R1D2, R1R2, xi_r
  real*8, dimension(:), allocatable :: weight_data1, weight_data2
  real*8, dimension(:), allocatable :: weight_random1, weight_random2

  character(20), external :: str
  character(len=500) :: data_filename2, data_filename1, output_filename
  character(len=500) :: random_filename1, random_filename2
  character(len=20) :: dmax_char, dmin_char, gridmin_char, gridmax_char
  character(len=20) :: ngrid_char, use_weights_char, rfilter_char, nthreads_char
  character(len=20) :: estimator, tracers_fileformat

  if (iargc() .ne. 15 ) then
    write(*,*) 'Some arguments are missing.'
    write(*,*) '1) data_filename1'
    write(*,*) '2) data_filename2'
    write(*,*) '3) random_filename1'
    write(*,*) '4) random_filename2'
    write(*,*) '5) output_filename'
    write(*,*) '6) dim1_min'
    write(*,*) '7) dim1_max'
    write(*,*) '8) rfilter'
    write(*,*) '9) ngrid'
    write(*,*) '10) gridmin'
    write(*,*) '11) gridmax'
    write(*,*) '12) estimator'
    write(*,*) '13) nthreads'
    write(*,*) '14) use_weights'
    write(*,*) '15) tracers_fileformat'
    write(*,*) ''
    stop
  end if

  ! read arguments from command line
  call getarg(1, data_filename1)
  call getarg(2, data_filename2)
  call getarg(3, random_filename1)
  call getarg(4, random_filename2)
  call getarg(5, output_filename)
  call getarg(6, dmin_char)
  call getarg(7, dmax_char)
  call getarg(8, rfilter_char)
  call getarg(9, ngrid_char)
  call getarg(10, gridmin_char)
  call getarg(11, gridmax_char)
  call getarg(12, estimator)
  call getarg(13, nthreads_char)
  call getarg(14, use_weights_char)
  call getarg(15, tracers_fileformat)
  
  ! convert string arguments to corresponding data types
  read(dmin_char, *) dim1_min
  read(dmax_char, *) dim1_max
  read(rfilter_char, *) rfilter
  read(ngrid_char, *) ngrid
  read(gridmin_char, *) gridmin
  read(gridmax_char, *) gridmax
  read(nthreads_char, *) nthreads
  read(use_weights_char, *) use_weights

    ! read files 
    if (trim(tracers_fileformat) == 'ascii') then
        if (use_weights == 1) then
            call read_catalogue_type2(data_filename1, data1, weight_data1, ndata1)
            call read_catalogue_type2(data_filename2, data2, weight_data2, ndata2)
            if (estimator .eq. 'LS') call read_catalogue_type2(random_filename1, random1, weight_random1, nrandom1)
            call read_catalogue_type2(random_filename2, random2, weight_random2, nrandom2)
        else
            call read_catalogue_type1(data_filename1, data1, weight_data1, ndata1)
            call read_catalogue_type1(data_filename2, data2, weight_data2, ndata2)
            if (estimator .eq. 'LS') call read_catalogue_type1(random_filename1, random1, weight_random1, nrandom1)
            call read_catalogue_type1(random_filename2, random2, weight_random2, nrandom2)

        end if
    else
        if (use_weights == 1) then
            call read_catalogue_type6(data_filename1, data1, weight_data1, ndata1)
            call read_catalogue_type6(data_filename2, data2, weight_data2, ndata2)
            if (estimator .eq. 'LS') call read_catalogue_type6(random_filename1, random1, weight_random1, nrandom1)
            call read_catalogue_type6(random_filename2, random2, weight_random2, nrandom2)
        else
            call read_catalogue_type5(data_filename1, data1, weight_data1, ndata1)
            call read_catalogue_type5(data_filename2, data2, weight_data2, ndata2)
            if (estimator .eq. 'LS') call read_catalogue_type5(random_filename1, random1, weight_random1, nrandom1)
            call read_catalogue_type5(random_filename2, random2, weight_random2, nrandom2)
        end if
    end if

  call linked_list(data2, ngrid, gridmin, gridmax, ll_data2, lirst_data2, rgrid)
  call linked_list(random2, ngrid, gridmin, gridmax, ll_random2, lirst_random2, rgrid)

  ! calculate number counts around each centre
  allocate(D1D2(ndata1))
  allocate(D1R2(ndata1))
  allocate(xi_r(ndata1))
  D1D2 = 0
  D1R2 = 0
  if (estimator .eq. 'LS') then
    allocate(R1R2(nrandom1))
    allocate(R1D2(nrandom1))
    R1R2 = 0
    R1D2 = 0
  end if
  
  ndif = int(dim1_max / rgrid + 1.)
  dim1_min2 = dim1_min ** 2
  dim1_max2 = dim1_max ** 2

  call OMP_SET_NUM_THREADS(nthreads)


  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ipx, ipy, ipz, &
  !$OMP& ix, iy, iz, ii, disx, disy, disz, dis2)
  do i = 1, ndata1
    ipx = int((data1(1, i) - gridmin) / rgrid + 1.)
    ipy = int((data1(2, i) - gridmin) / rgrid + 1.)
    ipz = int((data1(3, i) - gridmin) / rgrid + 1.)

    ! loop over cells around each centre
    do ix = ipx - ndif, ipx + ndif, 1
      do iy = ipy - ndif, ipy + ndif, 1
        do iz = ipz - ndif, ipz + ndif, 1
          if ((ix-ipx)**2 + (iy-ipy)**2 + (iz-ipz)**2 .gt. (ndif + 1)**2) cycle

          ! loop over data2 in each cell
          ii = lirst_data2(ix, iy, iz)
          if (ii .ne. 0) then
            do
              ii = ll_data2(ii)
              disx = data2(1, ii) - data1(1, i)
              disy = data2(2, ii) - data1(2, i)
              disz = data2(3, ii) - data1(3, i)

              dis2 = disx * disx + disy * disy + disz * disz

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                D1D2(i) = D1D2(i) + weight_data1(i) * weight_data2(ii)
              end if
  
              if(ii.eq.lirst_data2(ix, iy, iz)) exit
  
            end do
          end if

          ! loop over random2 in each cell
          ii = lirst_random2(ix, iy, iz)
          if (ii .ne. 0) then
            do
              ii = ll_random2(ii)
              disx = random2(1, ii) - data1(1, i)
              disy = random2(2, ii) - data1(2, i)
              disz = random2(3, ii) - data1(3, i)
 
              dis2 = disx * disx + disy * disy + disz * disz

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                D1R2(i) = D1R2(i) + weight_data1(i) * weight_random2(ii)
              end if
  
              if(ii .eq. lirst_random2(ix, iy, iz)) exit
  
            end do
          end if
        end do
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  if (estimator .eq. 'LS') then 
    ! Loop over randoms # 1
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, ipz, &
    !$OMP ix, iy, iz, disx, disy, disz, dis2)
    do i = 1, nrandom1
      ipx = int((random1(1, i) - gridmin) / rgrid + 1.)
      ipy = int((random1(2, i) - gridmin) / rgrid + 1.)
      ipz = int((random1(3, i) - gridmin) / rgrid + 1.)

      do ix = ipx - ndif, ipx + ndif, 1
        do iy = ipy - ndif, ipy + ndif, 1
          do iz = ipz - ndif, ipz + ndif, 1 
            if ((ix - ipx)**2 + (iy - ipy)**2 + (iz - ipz)**2 .gt. (ndif + 1)**2) cycle

            ii = lirst_random2(ix, iy, iz)
            if (ii .ne. 0) then
              do
                ii = ll_random2(ii)
                disx = random2(1, ii) - random1(1, i)
                disy = random2(2, ii) - random1(2, i)
                disz = random2(3, ii) - random1(3, i)

                dis2 = disx * disx + disy * disy + disz * disz

                if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                  R1R2(i) = R1R2(i) + weight_random1(i) * weight_random2(ii)
                end if

                  if(ii .eq. lirst_random2(ix, iy, iz)) exit

              end do
            end if

            ii = lirst_data2(ix, iy, iz)
            if (ii .ne. 0) then
              do
                ii = ll_data2(ii)
                disx = data2(1, ii) - random1(1, i)
                disy = data2(2, ii) - random1(2, i)
                disz = data2(3, ii) - random1(3, i)

                dis2 = disx * disx + disy * disy + disz * disz

                if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                  R1D2(i) = R1D2(i) + weight_random1(i) * weight_data2(ii)
                end if

                  if(ii .eq. lirst_data2(ix, iy, iz)) exit

              end do
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end if


  ! Normalize pair counts
  D1D2 = D1D2 * 1. / (SUM(weight_data1) * SUM(weight_data2))
  D1R2 = D1R2 * 1. / (SUM(weight_data1) * SUM(weight_random2))
  if (estimator .eq. 'LS') then
    R1R2 = R1R2 * 1. / (SUM(weight_random1) * SUM(weight_random2))
    R1D2 = R1D2 * 1. / (SUM(weight_random1) * SUM(weight_data2))
  end if

  ! Calculate density contrast
  if (estimator .eq. 'DP') then
    xi_r = (D1D2 / D1R2) - 1
  else if (estimator .eq. 'LS') then
    xi_r = (D1D2 - D1R2 - R1D2 + R1R2) / R1R2
  else
    write(*,*) 'Estimator for the correlation function was not recognized.'
    stop
  end if
 
  ! write output  
  open(12, file=output_filename, status='replace', form='unformatted')
  write(12) ndata1
  write(12) xi_r

  end program tophat_filter
  
  
