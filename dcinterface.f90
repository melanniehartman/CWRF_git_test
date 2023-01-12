MODULE module_daycent_cwrf_interface
    
! This module
! - defines input type (dcInputs)
! - defines output type (dcOutputs)
! - contains subroutine daycent_driver 
! - contains subroutine daycent_init
!   + subroutine daycent_init requires the following files for
!     initialization (in the future some of these files may be with 
!     information coming directly from CWRF):
!     - spin.bin: binary file that is used to initialize all state variables.
!       This binary file is generated by a spinup simulation for each grid cell.
!     - test.sch: scheduling information
!     - fix.100:  "fixed" parameters that are usually consistent for each site.
!     - crop.100: crop parameters
!     - tree.100: tree parameters
!     - other *.100 files: required to implement management options in the schedule file.
!     - site.100: weather statitics and site-specific information
!     - sitepar.in: more site-specific information
!     - soils.in: soil properties and soil layer structure

    integer :: CWRFLYRS
    parameter (CWRFLYRS = 11)

    type dcInputs

        integer ::           &
            curday,          &! day of year (1..366)
            month,           &! month (1..12)
            year              ! year

        double precision ::    &
            Ndep,              &! N deposition (wet + dry) (gN/m2/day)
            tmin,              &! minimum temperature for the day (C)
            tmax,              &! maximum temperature for the day (C)
            tmean,             &! mean temperature for the day, derived from hourly temperatures (C)
            ppt,               &! total precipitation for the day (cm)
            sradKJ,            &! incoming solar radiation (KJ/m2/day)
            parKJ,             &! photosynthetically active radiation (KJ/m2/day) (Need Einsteins/m2/day)
            rhumid,            &! relative humidity (%)
            windsp,            &! wind speed (mph)
            uvb,               &! UV-b radiation (units?)
            vpd,               &! vapor pressure deficit (kPa)
            snow,              &! snow cover (cm SWE)
            ambientCO2ppm       ! ambient atmospheric CO2 concentration (ppm)

    end type dcInputs

    type dcOutputs

    ! ATTENTION: Need to add P variables eventually.
    ! Added somtc, srfcomtc, Ra, and GPP (-mdh 1/10/2022).

        double precision ::    &
            aglivc,            &! live crop above-ground biomass (gC/m2)
            bglivc,            &! live crop juvenile + mature fine root biomass (gC/m2)
            stdedc,            &! dead crop standing above-ground biomass (gC/m2)
            rleavc,            &! live tree leaf biomass (gC/m2)
            frootc,            &! live tree fine root biomass (gC/m2)
            fbrchc,            &! live tree fine branch biomass (gC/m2)
            rlwodc,            &! live tree large wood biomass (gC/m2)
            crootc,            &! live tree coarse root biomass (gC/m2)
            somsc,             &! total soil (below ground) organic matter C (excluding litter C) (gC/m2)
            somtc,             &! total soil (below ground) organic matter C (including litter C) (gC/m2)
            srfcomtc,          &! total surface (above ground) organic matter C (including litter C) (gC/m2)
            totsysc,           &! total system C (gC/m2)
            nitrate(CWRFLYRS), &! soil nitrate by layer (indicator of fertility)(gNO3-N/m2)
            ammonium,          &! soil ammonium pool (indicator of fertility) (gNH4-N/m2)
            somseN,            &! total soil organic matter N (indicator of fertility) (gN/m2)
            Nuptake,           &! net N uptake for crops and trees (gN/m2/day)
            cropANPP,          &! above-ground net carbon uptake for crops (gC/m2/day)
            cropBNPP,          &! below-ground net carbon uptake for crops (gC/m2/day)
            treeANPP,          &! above-ground net carbon uptake for trees (gC/m2/day)
            treeBNPP,          &! below-ground net carbon uptake for trees (gC/m2/day)
            cgrain,            &! amount of crop grain C harvest (gC/m2/day)
            crmvst,            &! amount of crop stover/straw C harvest (gC/m2/day)
            egrainN,           &! amount of crop grain N harvest (gN/m2/day)
            ermvstN,           &! amount of crop stover/straw N harvest (gN/m2/day)
            ctubes,            &! amount of below ground root harvest (gC/m2/day)
            N2Oflux,           &! N2O emissions from nitrification and denitrification (gN2O-N/m2/day)
            NOflux,            &! NOx emissions from nitrification and denitrification (gNO-N/m2/day)
            CH4oxid,           &! CH4 oxidation (g CH4-C/m2/day)
            CH4prod,           &! CH4 production (g CH4-C/m2/day)
            GPP,               &! gross primary production a.k.a. gross canopy photosynthesis (gC/m2/day)
            Ra,                &! autotrophic respiration (gC/m2/day)
            Rh,                &! heterotrophic respiration (gC/m2/day)
            cropLAI,           &! LAI of live shoot biomass (calculated from aglivc) (m2/m2)
            treeLAI,           &! LAI of live tree leaf biomass (calculated from rleavc) (m2/m2)
            canopyHeight,      &! canopy height (m)
            dswc(CWRFLYRS),    &! change in soil water content by layer from previous day (+ is gain) (cm)
            dsnow,             &! change in SWE from previous day (+ is net accumulation) (cm)
            dtotsysC,          &! change in total system C from previous day (+ is gain)(gC/m2)
            dsomC,             &! change in soil organic matter C from previous day (+ is gain)(gC/m2)
            outflow,           &! runoff + baseflow (cm/day)
            intercpt,          &! amount of water intercepted by litter and vegetation, assumed to evaporate (cm/day)
            evapdly,           &! evaporation from H2O from soil (cm/day)
            trandly,           &! transpiration (cm/day)
            sublim,            &! sublimation (cm/day)
            volpl,             &! volatilization of N from plants (for indirect N2O) (gN/m2/day)
            leachNO3N,         &! amount of NO3-N leached into stream (for indirect N2O) (gN/m2/day)
            rootdens(CWRFLYRS)  ! frac. of roots in each soil layer (0.0-1.0)
       
    end type dcOutputs

CONTAINS

!--------------------------------------------------------------------------------
    SUBROUTINE daycent_init(luindex,fracGridCell,dcinput,dcgridpt,scellID)

        USE module_daycent_globalVars
        implicit none

        integer, INTENT(IN) :: luindex                 ! CWRF land use index
        character(len=6), INTENT(IN)   :: scellID      ! Cell ID
        real(8), INTENT(IN)            :: fracGridCell ! fraction of grid cell represented by this point
        type (dcInputs), INTENT(IN)    :: dcinput
        type (gridCell), INTENT(INOUT) :: dcgridpt

        !Local Variables
        integer          :: doinit, endOfDay
        doinit = 1
        endOfDay = 1  !Trigger common blocks to be saved to global data structure

        write(*,*) '----------------------------------------------------------------------------'
        write(*,110) 'INITIALIZE:  LUindex=', luindex, '   cellID=', scellID,  &
                      '   curday=', dcinput%curday, '   year=', dcinput%year
110 format(a21,i2, a10,a6, a10,i3, a8,i4) 

        dcgridpt%scellID = scellID
        dcgridpt%luindex = luindex
        dcgridpt%fracGridCell = fracGridCell

        call getsetglobvars(dcinput,dcgridpt,doinit,endOfDay) 

        !Reinitialize these variables after spinup
        dcgridpt%plt1%cgrain = 0.0
        dcgridpt%plt1%crmvst = 0.0
        dcgridpt%plt1%egrain(:) = 0.0
        dcgridpt%plt1%ermvst(:) = 0.0


    END SUBROUTINE daycent_init

!--------------------------------------------------------------------------------
    SUBROUTINE daycent_driver(dcinput,dcgridpt,dcoutput)

        USE module_daycent_globalVars
        implicit none

    ! Subroutine arguments
    !  doinit = do model initialization from spinup run

        type (dcInputs),  INTENT(IN)    :: dcinput 
        type (gridCell),  INTENT(INOUT) :: dcgridpt 
        type (dcOutputs), INTENT(OUT)   :: dcoutput 
    
    !Local Variables
    !  endOfDay: 0=copy global vars into common block vars at the beginning of the day
    !            1=copy common block vars into global vars at the end of the day 

    integer endOfDay, doinit

    write(*,*) '---------------------------------------------------'
    write(*,111) 'cellID=', dcgridpt%scellID, ' curday=', dcinput%curday, ' month=', dcinput%month, ' year=', dcinput%year
111 format(a7,a6, a8,i3, a7,i2, a6,i4) 

    !Copy values from global data structure into common block variables
    doinit = 0
    endOfDay = 0

    call getsetglobvars(dcinput,dcgridpt,doinit,endOfDay)

    call daycent_eachday(dcinput%curday,dcinput%Ndep)

    !Save common block values to global data structure
    doinit = 0
    endOfDay = 1

    call getsetglobvars(dcinput,dcgridpt,doinit,endOfDay) 

    !Assign values to dcoutput variables to send back to CWRF
    call assignoutputvars(dcgridpt,dcoutput)

    return

    END SUBROUTINE daycent_driver

END MODULE module_daycent_cwrf_interface

