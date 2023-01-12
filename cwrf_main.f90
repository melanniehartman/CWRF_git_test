
PROGRAM cwrf_test

!HOLAAAAA Main program to test the CWRF/DayCent interface
! * calls GetWthNcFileDim() to get the dimensions of the netcdf weather files
! * calls allocGridCells() to allocate memory for all grid cells
! * calls ReadWthNcFile() each year for each weather variable to read meteorological drivers
! * calls Read100Files() to read in DayCent *.100 parameter files during initialization.
! * calls daycent_init() for each gridcell to initialize it
! * calls CreateDayCentNcFile_daily() to create output NetCDF file (optional)
! * calls daycent_driver each day for each gridcell
! * calls WriteDayCentNcFile_annually() to write to daily output NetCDF file once a year (optional)

!  NetCDF File Name                       Variable(dimensions)                        Description
!  
!  tmin_cwrf/air.2m.yyyy.min.out.cwrf.nc  double air(time, south_north, west_east)    daily minimum temperature at 2m (°C)
!  tmax_cwrf/air.2m.yyyy.max.out.cwrf.nc  double air(time, south_north, west_east)    daily maximum temperature at 2m (°C)
!  apcp_cwrf/apcp.yyyy.out.cwrf.nc        double apcp(time, south_north, west_east)   total daily precipitation (cm)
!  dswrf_cwrf/dswrf.yyyy.out.cwrf.nc      double dswrf(time, south_north, west_east)  total solar radiation input (KJ/m2/day)
!  rhum_cwrf/rhum.2m.yyyy.out.cwrf.nc     double rhum(time, south_north, west_east)   mean daily relative humidity (%)
!  twnd_cwrf/twnd.10m.yyyy.out.cwrf.nc    double twnd(time, south_north, west_east)   mean daily wind speed at 10m (m/sec)
!  vpd_cwrf/vpd.2m.yyyy.out.cwrf.nc       double vpd(time, south_north, west_east)    vapor pressure deficit (kPa)
!  weasd_cwrf/weasd.yyyy.out.cwrf.nc      double weasd(time, south_north, west_east)  Snow cover (cm SWE)

    ! USES
    USE modWeather                          !create_cwrf_weather.f90
    USE modConstants                        !fconstants.f90
    USE modWrtNcFiles                       !wrtncfiles.f90
    USE module_daycent_globalVars           !globalVariables.f90
    USE module_daycent_cwrf_interface       !dcinterface.f90 defining dcinput and dcoutput
    USE DayCentParamLib
    IMPLICIT NONE
 
    ! Local Variables
    !type (gridCell), allocatable :: gridpt(:)
    integer(8) :: ipt, nGridCells
    integer :: dayofmo, month
    integer :: ilat, ilon, iyr, idy               ! loop indices
    integer :: itime                              ! timestep count
    integer :: iFirstYr, iLastYr                  ! First and last years of the weather record
    integer :: id, indx                           ! temporary variables
    integer :: ntime                              ! number of times in each netCDF weather file (expecting 365 or 366 days/file)
    integer :: nlat                               ! number of latitudes (south to north)
    integer :: nlon                               ! number of longitudes (west to east)
    character(len=150) :: ncPathNameWth           ! Path to NetCDF input weather files
    character(len=150) :: ncFileNameWth           ! NetCDF input weather files
    character(len=150) :: ncFileNameOut           ! Daily DayCent NetCDF outpt file 
    character(len=150) :: fNameLandCells          ! cwrfLandCells.csv (input file with list of all land cells)
    character(len=6) :: scellID                   ! cellID converted to a string
    character(len=4) :: syear                     ! iyr converted to a string
    character(len=5) :: varname                   ! name of variable to retrieve from netcdf file
    real(4) :: lat, lon                           ! temporary variables
    !!real(8) :: solrad, sradWm2                  ! srad (KJ/m2/day) is converted to langleys/day and W/m2 in weather file
    real(8) :: windmph                            ! windspeed in m.p.h
    real(8) :: fracGridCell                       ! Grid Cell fraction of each gridpt (0.0-1.0)
    logical :: wrtNetCDFfile
    type (dcInputs)  :: dcinput
    type (dcOutputs) :: dcoutput

    external getmonth
    integer :: getmonth
 
    ! Soil properties stored in NetCDF files (full time series)
    integer, allocatable :: cellID(:,:)           ! cellID(nlon,nlat) = 1000*ilat + ilon
    integer, allocatable :: cellMissing(:,:)      ! cellMissing(nlon,nlat) 1=nodata, 0=data
    ! xlat and xlon in netCDF files are float (real(4)) 
    integer(8), allocatable :: luindex(:,:)       ! CWRF Land Use Index 
    real(4), allocatable :: xlat(:,:)             ! xlat(nlon,nlat) latitude (decimal degrees)
    real(4), allocatable :: xlon(:,:)             ! xlon(nlon,nlat) longitude (decimal degrees)
    ! mylat and mylon are read from fNameLandCells
    real(4), allocatable :: mylat(:,:)            ! mylat(nlon,nlat) latitude (decimal degrees)
    real(4), allocatable :: mylon(:,:)            ! mylon(nlon,nlat) longitude (decimal degrees)
    !meteorological variables in netCDF files are double (real(8))
    real(8), allocatable :: precip(:,:,:)         ! precip(nlon,nlat,ntime) total daily precipitation (cm/day)
    real(8), allocatable :: tmax(:,:,:)           ! tmax(nlon,nlat,ntime)   daily maximum temperature at 2m (°C)
    real(8), allocatable :: tmin(:,:,:)           ! tmin(nlon,nlat,ntime)   daily minimum temperature at 2m (°C)
    real(8), allocatable :: srad(:,:,:)           ! srad(nlon,nlat,ntime)   total solar radiation input (KJ/m2/day)
    real(8), allocatable :: windsp(:,:,:)         ! windsp(nlon,nlat,ntime) mean daily wind speed at 10m (m/sec)
    real(8), allocatable :: rhum(:,:,:)           ! rhum(nlon,nlat,ntime)   mean daily relative humidity (%)
    real(8), allocatable :: vpd(:,:,:)            ! vpd(nlon,nlat,ntime)    vapor pressure deficit (kPa)
    real(8), allocatable :: weasd(:,:,:)          ! weasd(nlon,nlat,ntime)  SWE (cm)

    ! Daily output variables for NetCDF file, -mdh 10/8/2020.
    real(8), allocatable:: aglivc(:,:,:)
    real(8), allocatable:: bglivc(:,:,:)
    real(8), allocatable:: rleavc(:,:,:)
    real(8), allocatable:: frootc(:,:,:)
    real(8), allocatable:: fbrchc(:,:,:)
    real(8), allocatable:: rlwodc(:,:,:)
    real(8), allocatable:: crootc(:,:,:)
    real(8), allocatable:: cropANPP(:,:,:)
    real(8), allocatable:: cropBNPP(:,:,:)
    real(8), allocatable:: cgrain(:,:,:)
    real(8), allocatable:: crmvst(:,:,:)
    real(8), allocatable:: treeANPP(:,:,:)
    real(8), allocatable:: treeBNPP(:,:,:)
    real(8), allocatable:: somsc(:,:,:)
    real(8), allocatable:: totsysc(:,:,:)
    real(8), allocatable:: Nuptake(:,:,:)
 
!   character(len=40) :: s1, s2
    character(len=40) :: FMT
    integer :: luindx, cltveg, natveg

    !Get nlon and nlat dimensions.  Any netcdf file will do.

    ncPathNameWth = "/data/uvb/common/DayCent/melannie/UV-B/NARR_CWRF_grid/"
    ncFileNameWth = trim(ncPathNameWth) // "apcp_cwrf/apcp.1979.out.cwrf.nc"
    call GetWthNcFileDim(ncFileNameWth, ntime, nlat, nlon)
 
    ! The size of these arrays does not change each year
    allocate(cellID(nlon,nlat))
    allocate(cellMissing(nlon,nlat))
    allocate(luindex(nlon,nlat))
    allocate(xlat(nlon,nlat))
    allocate(xlon(nlon,nlat))
    allocate(mylat(nlon,nlat))
    allocate(mylon(nlon,nlat))

    cellMissing(:,:) = 1  !Assume cell is missing unless otherwise indicated in fNameLandCells.
    luindex(:,:) = 0
    cellID(:,:) = 0
   
    ! Generate weather files for land cells only (cellMissing = 0)
    ! File fNameLandCells contains cellIDs and missing status for every cell in the grid
    
    fNameLandCells = 'testcwrfLandCells.csv'
    open(unit=20, FILE=fNameLandCells)
    read(20,*)   ! Read past header line
 
    ! Read contents of fNameLandCells.  Initialize cellID and cellMissing arrays
    nGridCells = 0
    222 continue
!       FMT = "(4(i8),2(f15.6),3(i3),2(a40))"
!       read(20,FMT,end=99) id,ilat,ilon,indx,lat,lon,luindx,natveg,cltveg
!       write(*,110) id, ilat, ilon, indx, lat, lon, luindx, natveg, cltveg
!110 format(i6,1x, i3,1x, i3,1x, i1,1x, f7.2,1x, f7.2,1x, i2,1x, i2)
        FMT = "(4(i8),2(f15.6),2(i3),2(a40))"
        read(20,FMT,end=99) id,ilat,ilon,indx,lat,lon,luindx,natveg
        write(*,110) id, ilat, ilon, indx, lat, lon, luindx, natveg
110 format(i6,1x, i3,1x, i3,1x, i1,1x, f7.2,1x, f7.2,1x, i2,1x, i2)
        if (indx == 0) nGridCells = nGridCells + 1
        cellID(ilon,ilat) = id
        cellMissing(ilon,ilat) = indx
        luindex(ilon,ilat) = luindx
        mylat(ilon,ilat) = lat
        mylon(ilon,ilat) = lon
        write(*,111) 'fNameLandCells cellID(', ilon, ',', ilat,')=', cellID(ilon,ilat)
        write(*,112) 'fNameLandCells luindex(', ilon, ',', ilat,')=', luindex(ilon,ilat)
        write(*,*)
111 format(a22,i3,a1,i3,a3,i6)
112 format(a23,i3,a1,i3,a3,i2)
    goto 222
    
    99 continue
    close(unit=20)
 
    call allocGridCells(nGridCells)

    iFirstYr = 1979  ! First year of the simulation
    iLastYr = 1980   ! Final year of the simulation
    !iLastYr = 2015
    wrtNetCDFfile = .TRUE. ! Write daily gridded output to a NetCDF file? One file for each year.
!    wrtNetCDFfile = .FALSE. ! Write daily gridded output to a NetCDF file? One file for each year.

    do iyr = iFirstYr, iLastYr
     
        write(syear,'(i4)') iyr
  
        ! Determine dimensions of netCDF file each year since ntime can change
  
        ncFileNameWth =  trim(ncPathNameWth) // "apcp_cwrf/apcp." // syear // ".out.cwrf.nc"
        call GetWthNcFileDim(ncFileNameWth, ntime, nlat, nlon)
     
        ! Allocate variables to store meteorological values
        ! Dimension order should be the opposite of that seen in the netCDF files
        ! since arrays in FORTRAN are stored in column major order.

        ! The size of these arrays changes each year
        allocate(precip(nlon,nlat,ntime))
        allocate(tmax(nlon,nlat,ntime))
        allocate(tmin(nlon,nlat,ntime))
        allocate(srad(nlon,nlat,ntime))
        allocate(windsp(nlon,nlat,ntime))
        allocate(rhum(nlon,nlat,ntime))
        allocate(vpd(nlon,nlat,ntime))
        allocate(weasd(nlon,nlat,ntime))
        allocate(aglivc(nlon,nlat,ntime))
        allocate(bglivc(nlon,nlat,ntime))
        allocate(rleavc(nlon,nlat,ntime))
        allocate(frootc(nlon,nlat,ntime))
        allocate(fbrchc(nlon,nlat,ntime))
        allocate(rlwodc(nlon,nlat,ntime))
        allocate(crootc(nlon,nlat,ntime))
        allocate(cropANPP(nlon,nlat,ntime))
        allocate(cropBNPP(nlon,nlat,ntime))
        allocate(cgrain(nlon,nlat,ntime))
        allocate(crmvst(nlon,nlat,ntime))
        allocate(treeANPP(nlon,nlat,ntime))
        allocate(treeBNPP(nlon,nlat,ntime))
        allocate(somsc(nlon,nlat,ntime))
        allocate(totsysc(nlon,nlat,ntime))
        allocate(Nuptake(nlon,nlat,ntime))
     
        precip(:,:,:) = dMISSING_VAL
        tmax(:,:,:) = dMISSING_VAL
        tmin(:,:,:) = dMISSING_VAL
        srad(:,:,:) = dMISSING_VAL
        windsp(:,:,:) = dMISSING_VAL
        rhum(:,:,:) = dMISSING_VAL
        vpd(:,:,:) = dMISSING_VAL
        weasd(:,:,:) = dMISSING_VAL
        aglivc(:,:,:) = dMISSING_VAL
        bglivc(:,:,:) = dMISSING_VAL
        rleavc(:,:,:) = dMISSING_VAL
        frootc(:,:,:) = dMISSING_VAL
        fbrchc(:,:,:) = dMISSING_VAL
        rlwodc(:,:,:) = dMISSING_VAL
        crootc(:,:,:) = dMISSING_VAL
        cropANPP(:,:,:) = dMISSING_VAL
        cropBNPP(:,:,:) = dMISSING_VAL
        cgrain(:,:,:) = dMISSING_VAL
        crmvst(:,:,:) = dMISSING_VAL
        treeANPP(:,:,:) = dMISSING_VAL
        treeBNPP(:,:,:) = dMISSING_VAL
        somsc(:,:,:) = dMISSING_VAL
        totsysc(:,:,:) = dMISSING_VAL
        Nuptake(:,:,:) = dMISSING_VAL
     
        ncFileNameWth =  trim(ncPathNameWth) // "apcp_cwrf/apcp." // syear // ".out.cwrf.nc"
        varname = "apcp"
        call ReadWthNcFile(ncFileNameWth, varname, ntime, nlat, nlon, xlat, xlon, precip)
     
        ncFileNameWth =  trim(ncPathNameWth) // "tmax_cwrf/air.2m." // syear // ".max.out.cwrf.nc"
        varname = "air"
        call ReadWthNcFile(ncFileNameWth, varname, ntime, nlat, nlon, xlat, xlon, tmax)
     
        ncFileNameWth =  trim(ncPathNameWth) // "tmin_cwrf/air.2m." // syear // ".min.out.cwrf.nc"
        varname = "air"
        call ReadWthNcFile(ncFileNameWth, varname, ntime, nlat, nlon, xlat, xlon, tmin)
     
        ncFileNameWth =  trim(ncPathNameWth) // "dswrf_cwrf/dswrf." // syear // ".out.cwrf.nc"
        varname = "dswrf"
        call ReadWthNcFile(ncFileNameWth, varname, ntime, nlat, nlon, xlat, xlon, srad)
     
        ncFileNameWth =  trim(ncPathNameWth) // "twnd_cwrf/twnd.10m." // syear // ".out.cwrf.nc"
        varname = "twnd"
        call ReadWthNcFile(ncFileNameWth, varname, ntime, nlat, nlon, xlat, xlon, windsp)
     
        ncFileNameWth =  trim(ncPathNameWth) // "rhum_cwrf/rhum.2m." // syear // ".out.cwrf.nc"
        varname = "rhum"
        call ReadWthNcFile(ncFileNameWth, varname, ntime, nlat, nlon, xlat, xlon, rhum)
     
        ncFileNameWth =  trim(ncPathNameWth) // "vpd_cwrf/vpd.2m." // syear // ".out.cwrf.nc"
        varname = "vpd"
        call ReadWthNcFile(ncFileNameWth, varname, ntime, nlat, nlon, xlat, xlon, vpd)

        ncFileNameWth =  trim(ncPathNameWth) // "weasd_cwrf/weasd." // syear // ".out.cwrf.nc"
        varname = "weasd"
        call ReadWthNcFile(ncFileNameWth, varname, ntime, nlat, nlon, xlat, xlon, weasd)

        if (iyr == iFirstYr) then

            ! Store parameters sets from *.100 files (except fix.100 and site.100) 
            ! in parameter library. These parameters are common to all gridcells. -mdh 8/22/2018
            call Read100Files()

            ipt = 0
            do ilon = 1, nlon
                do ilat = 1, nlat
                    if (cellMissing(ilon,ilat) .eq. 0) then
                        ipt = ipt + 1
                        write(scellID,'(i6)') cellID(ilon,ilat)
                        dcinput%curday = 1
                        dcinput%year = iyr 
                        fracGridCell = 1.0
                        !Create a 6-digit cellID
                        if (cellID(ilon,ilat) .lt. 100000)  scellID(1:1) = '0'
                        if (cellID(ilon,ilat) .lt. 10000)   scellID(2:2) = '0'
                        write(*,113) ' cwrf_main: cellID(', ilon,ilat,')=', cellID(ilon,ilat), 'scellID=', scellID
                        write(*,114) ' cwrf_main: mylat(', ilon, ilat,')=', mylat(ilon,ilat), &
                                                 'xlat(', ilon, ilat,')=', xlat(ilon,ilat)
                        write(*,114) ' cwrf_main: mylon(', ilon, ilat,')=', mylon(ilon,ilat), &
                                                'xlon(', ilon, ilat,')=', xlon(ilon,ilat)
113 format(a19,i3,',',i3,a2,i6,2x,a8,a6)
114 format(a18, i3,',',i3, a2,f9.4,2x, a5, i3,',',i3, a2,f9.4)
                        if (abs(xlat(ilon,ilat) - mylat(ilon,ilat)) > 0.001) then
                          write(*,*) 'ERROR. Latitudes between netcdf files and input file are not consistent.'
                          STOP
                        endif
                        if (abs(xlon(ilon,ilat) - mylon(ilon,ilat)) > 0.001) then
                          write(*,*) 'ERROR. Longitudes between netcdf files and input file are not consistent.'
                          STOP
                        endif
    
                        call daycent_init(luindex(ilon,ilat),fracGridCell,dcinput,gridpt(ipt),scellID)  
    
                    endif
                enddo
            enddo
        endif

        if (wrtNetCDFfile) then
            ! Create a new netcdf file each year
            ncFileNameOut = "./daycent_cwrf_" // syear // ".nc"
            call CreateDayCentNcFile_daily(ncFileNameOut,ntime,nGridCells,nlat,nlon,cellID,cellMissing,luindex,xlat,xlon)
        endif

        do idy = 1, ntime
            ipt = 0
            itime = idy  ! Timestep count
            do ilon = 1, nlon
                do ilat = 1, nlat
         
                    if (cellMissing(ilon,ilat) .eq. 0) then
    
                        ipt = ipt + 1

                        call GetMonthDay(iyr, idy, month, dayofmo)
                        !daylen = DayLength(idy, xlat(ilon,ilat))
                        !!solrad = srad(ilon,ilat,idy)*1000/(CM2_PER_M2*JOULES_PER_CALORIE)  ! Convert KJ/m2/day to langleys/day (cal/cm2/day)
                        !!sradWm2 = srad(ilon,ilat,idy)*1000/(daylen*SEC_PER_HOUR)           ! Convert KJ/m2/day to mean daytime input (W/m2)
                        windmph = windsp(ilon,ilat,idy)*SEC_PER_HOUR/METERS_PER_MILE         ! Convert m/s to mph

                        ! Inputs to daycent_driver
                        dcinput%curday = idy                 ! currect day of year (1-366)
                        dcinput%month = month                ! current month (1-12)
                        dcinput%year = iyr                   ! current calendar year
                        dcinput%tmin = tmin(ilon,ilat,idy)   ! daily minimum temperature at 2m (°C)
                        dcinput%tmax = tmax(ilon,ilat,idy)   ! daily maximum temperature at 2m (°C)
                        dcinput%ppt =  precip(ilon,ilat,idy) ! total daily precipitation (cm)
                        dcinput%sradKJ = srad(ilon,ilat,idy) ! total solar radiation input (KJ/m2/day)
                        dcinput%parKJ = 0.5*dcinput%sradKJ   ! photosynthetically active solar radiation input (KJ/m2/day)
                        dcinput%rhumid = rhum(ilon,ilat,idy) ! mean daily relative humidity (%)
                        dcinput%windsp = windmph             ! mean daily wind speed at 10m (mph)
                        dcinput%vpd = vpd(ilon,ilat,idy)     ! vapor pressure deficit (kPa) (Should be daytime average!)
                        dcinput%snow = weasd(ilon,ilat,idy)  ! Snow cover (cm SWE), not currently used within daycent_driver
                        dcinput%ambientCO2ppm = 350          ! ambient CO2 concentration (ppm), replace with value from CWRF
                        dcinput%Ndep = -99.0                 ! Atmospheric N deposition (gN/m2/day). When set to negative value, 
                                                             !   DayCent computes Ndep.

                        call daycent_driver(dcinput,gridpt(ipt),dcoutput)

                        if (wrtNetCDFfile) then
                            !call WriteDayCentNcFile_daily(ncFileNameOut, dcoutput, ilat, ilon, iyr, idy, itime)

                            aglivc(ilon,ilat,itime) = dcoutput%aglivc
                            bglivc(ilon,ilat,itime) = dcoutput%bglivc
                            rleavc(ilon,ilat,itime) = dcoutput%rleavc
                            frootc(ilon,ilat,itime) = dcoutput%frootc
                            fbrchc(ilon,ilat,itime) = dcoutput%fbrchc
                            rlwodc(ilon,ilat,itime) = dcoutput%rlwodc
                            crootc(ilon,ilat,itime) = dcoutput%crootc
                            cropANPp(ilon,ilat,itime) = dcoutput%cropANPp
                            cropBNPP(ilon,ilat,itime) = dcoutput%cropBNPP
                            cgrain(ilon,ilat,itime) = dcoutput%cgrain
                            crmvst(ilon,ilat,itime) = dcoutput%crmvst
                            treeANPP(ilon,ilat,itime) = dcoutput%treeANPP
                            treeBNPP(ilon,ilat,itime) = dcoutput%treeBNPP
                            somsc(ilon,ilat,itime) = dcoutput%somsc
                            totsysc(ilon,ilat,itime) = dcoutput%totsysc
                            Nuptake(ilon,ilat,itime) = dcoutput%Nuptake
                        endif

                    endif  !cellMissing .eq. 0
     
                enddo ! ilat
            enddo ! ilon
        enddo !idy

        ! Write a year's worth of daily values to the netCDF output file. -mdh 10/8/2020
        if (wrtNetCDFfile) then 
            write(*,*)
            write(*,*) 'Writing results to ', ncFileNameOut
            call WriteDayCentNcFile_annually(ncFileNameOut, nlat, nlon, ntime, iyr, &
                aglivc, bglivc, rleavc, frootc, fbrchc, rlwodc, crootc, cropANPP, cropBNPP, &
                cgrain, crmvst, treeANPP, treeBNPP, somsc, totsysc, Nuptake)
        endif
 
        !Deallocate weather & output arrays each year since ntime (# days per year) depends on calendar year
        deallocate(precip)
        deallocate(tmax)
        deallocate(tmin)
        deallocate(srad)
        deallocate(windsp)
        deallocate(rhum)
        deallocate(vpd)
        deallocate(weasd)
        deallocate(aglivc)
        deallocate(bglivc)
        deallocate(rleavc)
        deallocate(frootc)
        deallocate(fbrchc)
        deallocate(rlwodc)
        deallocate(crootc)
        deallocate(cropANPP)
        deallocate(cropBNPP)
        deallocate(cgrain)
        deallocate(crmvst)
        deallocate(treeANPP)
        deallocate(treeBNPP)
        deallocate(somsc)
        deallocate(totsysc)
        deallocate(Nuptake)
 
    enddo ! iyr

    write(*,*)
    write(*,*) 'Simulation complete.'

END PROGRAM cwrf_test
