      subroutine getsetglobvars(dcinput,dcgridpt,doinit,endOfDay)

c   ANOTHER MODIF
c   Contains subroutine getsetglobvars that calls other subroutines for 
c   getting/setting global variables at the beginning/end of each daily timestep.
c   Subroutine getsetglobvars in called by daycent_init and daycent_driver

      USE module_daycent_globalVars
      USE module_daycent_cwrf_interface

      implicit none
c     include 'const.inc'
      include 'chrvar.inc'
      include 'comput.inc'
      include 'dovars.inc'
      include 'dynam.inc'
      include 'fertil.inc'
c     include 'forrem.inc'
      include 'isovar.inc'
      include 'ligvar.inc'
      include 'jday.inc'
      include 'monprd.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'photosyn.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
c     include 'potent.inc'
      include 'seq.inc'
      include 'schvar.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 't0par.inc'
      include 'waterbal.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

c ... Formal parameters (INTENT(IN))
c       doinit      = do model initialization from spinup run
c       endOfDay:   0=copy global vars into common block vars at the beginning of the day
c                   1=copy common block vars into global vars at the end of the day 
c       dcinput
c         curday   = day of year (1..366)
c         month    = current month (1..12)
c         year     = simulation year
c         tmin     = minimum daily temperature (C)
c         tmax     = maximum daily temperature (C)
c         ppt      = precipitation (cm)
c         solrad   = solar radiation (langleys/day)
c         srad     = solar radiation (W/m2, average over daylight period)
c         sradKJ   = solar radiation input (KJ/m2/day)
c         parKJ    = photosynthetically active radiation input (KJ/m2/day)
c         rhumida  = relative humidity (%)
c         windsp   = wind speed (m.p.h.)
c         vpd      = vapor pressure deficit (kPa)
c         snow     = snow cover (cm SWE)
c         ambientCO2ppm = ambient atmospheric CO2 concentration (ppm)

      type(dcInputs) :: dcinput
      type(gridCell) :: dcgridpt
      integer doinit
      integer endOfDay

c ....Interfaces to C functions

      INTERFACE

        SUBROUTINE save_c_vars(
     &    endOfDay, numlyrs, nelyrs, ubnd, lbnd, width, depth, dpthmn,
     &    dpthmx, bulkd2, fieldc, wiltpt, ecoeff, tcoeff, sandfrac,
     &    clayfrac, orgfrac, swclimit, satcond, pH2, swcfc, swcwp,
     &    swc, swcmin, minpot, wfps, sumecoeff, thetas, thetas_bd,
     &    psis, b, binverse, stmtemp, soiltavg, soiltmin, soiltmax,
     &    usexdrvrs, texture, hours_rain, drainlag, jdayStart, jdayEnd, 
     &    SnowFlag,watertable, sublimscale, reflec, albedo, fswcinit, 
     &    dmpflux, hpotdeep, ksatdeep, rlatitude, cldcov, tbotmn,
     &    tbotmx, dmp, timlag, Ncoeff, N2Oadjust_fc, N2Oadjust_wp,
     &    MaxNitAmt, netmn_to_no3, wfpsdnitadj, n2n2oadj, noxn2oadj, 
     &    elevation,
     &    site_slope, aspect, ehoriz, whoriz, dDO_fc, dDO_wp, 
     &    daylength,indewpt,soiltavgdy,par,isoilwkptr)

          INTEGER endOfDay
          INTEGER numlyrs
          INTEGER nelyrs
          INTEGER ubnd(*)
          INTEGER lbnd(*)
          DOUBLE PRECISION width(*)
          DOUBLE PRECISION depth(*)
          DOUBLE PRECISION dpthmn(*)
          DOUBLE PRECISION dpthmx(*)
          DOUBLE PRECISION bulkd2(*)
          DOUBLE PRECISION fieldc(*)
          DOUBLE PRECISION wiltpt(*)
          DOUBLE PRECISION ecoeff(*)
          DOUBLE PRECISION tcoeff(*)
          DOUBLE PRECISION sandfrac(*)
          DOUBLE PRECISION clayfrac(*)
          DOUBLE PRECISION orgfrac(*)
          DOUBLE PRECISION swclimit(*)
          DOUBLE PRECISION satcond(*)
          DOUBLE PRECISION pH2(*)
          DOUBLE PRECISION swcfc(*)
          DOUBLE PRECISION swcwp(*)
          DOUBLE PRECISION swc(*)
          DOUBLE PRECISION swcmin(*)
          DOUBLE PRECISION minpot(*)
          DOUBLE PRECISION wfps(*)
          DOUBLE PRECISION sumecoeff
          DOUBLE PRECISION thetas(*)
          DOUBLE PRECISION thetas_bd(*)
          DOUBLE PRECISION psis(*)
          DOUBLE PRECISION b(*)
          DOUBLE PRECISION binverse(*)
          DOUBLE PRECISION stmtemp(*)
          DOUBLE PRECISION soiltavg(*)
          DOUBLE PRECISION soiltmin(*)
          DOUBLE PRECISION soiltmax(*)
          INTEGER usexdrvrs
          INTEGER texture
          INTEGER hours_rain
          INTEGER drainlag
          INTEGER jdayStart
          INTEGER jdayEnd
          INTEGER SnowFlag
          INTEGER watertable(*)
          DOUBLE PRECISION sublimscale
          DOUBLE PRECISION reflec
          DOUBLE PRECISION albedo
          DOUBLE PRECISION fswcinit
          DOUBLE PRECISION dmpflux
          DOUBLE PRECISION hpotdeep
          DOUBLE PRECISION ksatdeep
          DOUBLE PRECISION rlatitude
          DOUBLE PRECISION cldcov(*)
          DOUBLE PRECISION tbotmn
          DOUBLE PRECISION tbotmx
          DOUBLE PRECISION dmp
          DOUBLE PRECISION timlag
          DOUBLE PRECISION Ncoeff
          DOUBLE PRECISION N2Oadjust_fc
          DOUBLE PRECISION N2Oadjust_wp
          DOUBLE PRECISION MaxNitAmt
          DOUBLE PRECISION netmn_to_no3
          DOUBLE PRECISION wfpsdnitadj
          DOUBLE PRECISION n2n2oadj
          DOUBLE PRECISION noxn2oadj
          DOUBLE PRECISION elevation
          DOUBLE PRECISION site_slope
          DOUBLE PRECISION aspect
          DOUBLE PRECISION ehoriz
          DOUBLE PRECISION whoriz
          DOUBLE PRECISION dDO_fc
          DOUBLE PRECISION dDO_wp
          DOUBLE PRECISION daylength(*)
          INTEGER indewpt
          DOUBLE PRECISION soiltavgdy(*)
          DOUBLE PRECISION par
          INTEGER isoilwkptr

        END SUBROUTINE save_c_vars

      END INTERFACE


c ... LOCAL VARIABLES
      logical ext
      integer imo, idy, iwrthdr, ilyr, iwrtcsv
      character*100 newbin, oldbin, schnam, soilnam
      character*100 sitnam
      character*100 filpath
      character*6 sID
      character*2 sLU

      data iwrthdr /1/
      save iwrthdr

c ... Local Variables that are part of the global data structure
c ... but not in a common block.  These are all(?) variables
c ... stored in C language structs.
      integer numlyrs,nelyrs,jdayStart,jdayEnd,hours_rain,
     &  drainlag,watertable(MONTHS+1),SnowFlag,
     &  lbnd(CMXLYR),ubnd(CMXLYR),indewpt,isoilwkptr
      double precision 
     &  sublimscale,reflec,albedo,timlag,
     &  elevation,site_slope,aspect,fswcinit,
     &  ehoriz,whoriz,hpotdeep,ksatdeep,
     &  dmp,dmpflux,cldcov(MONTHS+1),tbotmn,tbotmx,
     &  Ncoeff,N2Oadjust_fc,N2Oadjust_wp,MaxNitAmt, 
     &  netmn_to_no3,n2n2oadj,noxn2oadj,wfpsdnitadj,
     &  width(SWMXLYR),depth(SWMXLYR),bulkd2(SWMXLYR),
     &  dpthmn(SWMXLYR),dpthmx(SWMXLYR),
     &  fieldc(SWMXLYR),wiltpt(SWMXLYR),
     &  ecoeff(SWMXLYR),tcoeff(SWMXLYR),
     &  sandfrac(SWMXLYR),clayfrac(SWMXLYR),orgfrac(SWMXLYR),
     &  swclimit(SWMXLYR),satcond(SWMXLYR),pH2(SWMXLYR),
     &  swcfc(SWMXLYR),swcwp(SWMXLYR),swcmin(SWMXLYR),
     &  swc(SWMXLYR),wfps(SWMXLYR),minpot(SWMXLYR),sumecoeff,
     &  thetas(SWMXLYR),thetas_bd(SWMXLYR),
     &  psis(SWMXLYR),b(SWMXLYR),binverse(SWMXLYR),
     &  soiltmin(SWMXLYR),soiltmax(SWMXLYR),soiltavg(SWMXLYR),
     &  stmtemp(MAXSTLYR),soiltavgdy(7),dDO_fc,dDO_wp,
     &  rlatitude,par

c ... Set to 1 to write to dcoutput.csv file
c     if (dcgridpt%scellID .eq. '075095' .OR. 
c    &    dcgridpt%scellID .eq. '083073') then
c         iwrtcsv = 1
c     else
c         iwrtcsv = 0
c     endif
c     iwrtcsv = 0
      iwrtcsv = 1

      if (doinit .eq. 1) then

c ..... Model Initialization at first time step
c ..... ATTENTION: these are assigned in detiv_cwrf also.
        ext = .true.
        filpath = ' '
        write(sLU,'(i2)') dcgridpt%luindex
        !Create a 2-digit LU index
        if (dcgridpt%luindex .lt. 10) sLU(1:1) = '0'
        sID = dcgridpt%scellID
        write(schnam,'(a17,a6,a4)') './TRANS_SCH/test_', sID, '.sch'
        write(newbin,'(a18)') 'XXXXXXXXXXXXXX.bin'
        write(oldbin,'(a18,a6,a4)') './SPINUP_BIN/spin_',sID,'.bin'
        write(soilnam,'(a16,a6,a3)') './SOILSIN/soils_', sID, '.in'
        write(sitnam,'(a16,a6,a4)') './SITE_100/site_', sID, '.100'
        write(*,*) 'getsetglbvars: schnam =', schnam
        write(*,*) 'getsetglbvars: oldbin =', oldbin
        write(*,*) 'getsetglbvars: soilnam =', soilnam
        write(*,*) 'getsetglbvars: sitnam =', sitnam

c ..... Read spinup binary file, the schedule file (header and 
c ..... first block only), fix.100, site.100, crop.100, tree.100, 
c ..... and soils.in files
        call detiv_cwrf(ext,schnam,newbin,oldbin,soilnam,
     &                  sitnam,filpath)
 
c ..... Adjust parameters from crop.100 and fix.100 for daily production
c ..... if necessary. -mdh 1/95
c ..... I don't think this call is needed for CWRF coupling. -mdh 5/25/2016
C       call adjustpar

c ATTENTION: uncomment these lines when compiling without csa_main.f
        data (idysimo(imo), imo=1,12) 
     &       /31,28,31,30,31,30,31,31,30,31,30,31/
        data (ilstdy(imo),  imo=1,12) 
     &       /31,59,90,120,151,181,212,243,273,304,334,365/
        data (ifrstdy(imo), imo=1,12) 
     &       /1,32,60,91,121,152,182,213,244,274,305,335/

c ..... Initialize air temperature arrays for computing
c ..... running average temperatures at the start of simulation.
        prcann = 0.0
        do 10 imo = 1, 12
          do 20 idy = ifrstdy(imo), ilstdy(imo)
            tempmin(idy) = tmn2m(imo)
            tempmax(idy) = tmx2m(imo)
            avgtemp(idy) = 0.5 * (tmn2m(imo) + tmx2m(imo))
c           soiltavgdy(idy) = tmn2m(imo) 
20        continue
          prcann = prcann + precip(imo)
10      continue

c ..... tavewkqueue(7) and tavemthqueue(30) are the average air temperatures over
c ....  the past 7 and 30 days, respectively.  Initialize these circular
c ..... queues with mean annual air temperature for December (previous 7 and 30 days). 
c ..... -mdh 11/16/2022.
c ..... soiltavgdy(7) is the average soil temperature in the second soil layer over 
c ..... the past 7 days. -mdh 11/21/2022
        iwkptr = 0
        imthptr = 0
        isoilwkptr = 0
        do 15 idy = 1, 7
          tavewkqueue(idy) = 0.5 * (tmn2m(12) + tmx2m(12))
          soiltavgdy(idy) = tmn2m(12) 
15      continue
        do 25 idy = 1, 30
          tavemthqueue(idy) = 0.5 * (tmn2m(12) + tmx2m(12))
25      continue

        if (iwrthdr .eq. 1 .and. iwrtcsv .eq. 1) then

c ...     Write output file header just once. This file contains output for 
c ...     all grid cells.

          open(unit=123, file='dcoutput.csv', status='replace')
c                   10        20        30        40        50
c            123456789012345678901234567890123456789012345678901234567890
          write(123, 125) 
     &      'time,curday,month,ID,ptagc,ptbgc,aglivc,bglivcj,rleavc,', 
     &      'frootcj,rlwodc,crootc,metabc(1),metabc(2),strucc(1),',
     &      'strucc(2),som1c(1),som1c(2),som2c(1),som2c(2),som3c,',
     &      'agdefac,bgdefac,ammonium,nitrat(1),nitrat(2),nitrat(3),',
     &      'vswc(1),vswc(2),vswc(3),aglive(1),bglivej(1),soiltavg(1),',
     &      'soiltavg(2),soiltavg(3),tmin,tmax,ppt,sradKJ,solrad,',
     &      'rhumid,windsp,vpd,snow,evapdly,trandly,intrcpt,',
     &      'sublim,outflow,petdly,stemp,carbostg,tminrlN,minerlN1,',
     &      'minerlN2,minerlN3,minerlN4,minerlN5,minerlN6,minerlN7,',
     &      'cgrain,crmvst,egrainN,ermvstN'
125       format(a55,a52,a52,a58,a57,a52,a47,a54,a54,a29)
          close(123)
          iwrthdr = 0
        endif

        do ilyr = 1, SWMXLYR
          dcgridpt%prevdy%swc(ilyr) = 0.0
        end do
        dcgridpt%prevdy%snow = 0.0
        dcgridpt%prevdy%somsc = 0.0
        dcgridpt%prevdy%totsysc = 0.0

c ... End doinit==1
      endif

      month = dcinput%month
      par = dcinput%parKJ

c     Get/Set variables in all C structures. 
c     Function initsw must be called first (from detiv_cwrf)
c ... Dynamic variables (swc, wfps, soiltave, soiltmin, soiltmax, stmtemp) 

      if (endOfDay .eq. 1) then

c ..... Retrieve variables from C structures, return these to local
c ..... variables so they can be saved to global data structure
        call save_c_vars(
     &    endOfDay,numlyrs,nelyrs,ubnd,lbnd,width,depth,dpthmn,
     &    dpthmx,bulkd2,fieldc,wiltpt,ecoeff,tcoeff,sandfrac,
     &    clayfrac,orgfrac,swclimit,satcond,pH2,swcfc,swcwp,
     &    swc,swcmin,minpot,wfps,sumecoeff,thetas,thetas_bd,
     &    psis,b,binverse,stmtemp,soiltavg,soiltmin,soiltmax,
     &    usexdrvrs,texture,hours_rain,drainlag,jdayStart,jdayEnd,
     &    SnowFlag,watertable,sublimscale,reflec,albedo,fswcinit,
     &    dmpflux,hpotdeep,ksatdeep,rlatitude,cldcov,tbotmn,
     &    tbotmx,dmp,timlag,Ncoeff,N2Oadjust_fc,N2Oadjust_wp,
     &    MaxNitAmt,netmn_to_no3,wfpsdnitadj,n2n2oadj,noxn2oadj,
     &    elevation,
     &    site_slope,aspect,ehoriz,whoriz,dDO_fc,dDO_wp,
     &    daylength,indewpt,soiltavgdy,par,isoilwkptr)

      endif

c ... Static parameter values
c ... These are set by detiv, but need to be copied from common blocks
c ... to global data structure at the end of the day, and need to be
c ... copied from global structure to common blocks at the beginning 
c ....of the day except during initialization 

      call SaveGlobalFixParams(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  idef=idef,nsnfix=nsnfix,ntspm=ntspm,seed=seed,adep=adep,
     &  agppa=agppa,agppb=agppb,aneref=aneref,animpt=animpt,
     &  bgppa=bgppa,bgppb=bgppb,co2ppm=co2ppm,co2rmp=co2rmp,
     &  damr=damr,damrmn=damrmn,dec1=dec1,dec2=dec2,dec3=dec3,
     &  dec4=dec4,dec5=dec5,deck5=deck5,dligdf=dligdf,dresp=dresp,
     &  edepth=edepth,elitst=elitst,enrich=enrich,favail=favail,
     &  fleach=fleach,fwloss=fwloss,fxmca=fxmca,fxmcb=fxmcb,
     &  fxmxs=fxmxs,fxnpb=fxnpb,gremb=gremb,lhzf=lhzf,minlch=minlch,
     &  omlech=omlech,p1co2a=p1co2a,p1co2b=p1co2b,p2co2=p2co2,
     &  p3co2=p3co2,pabres=pabres,peftxa=peftxa,peftxb=peftxb,
     &  phesp=phesp,pligst=pligst,pmco2=pmco2,pmnsec=pmnsec,
     &  pmntmp=pmntmp,pmxbio=pmxbio,pmxtmp=pmxtmp,pparmn=pparmn,
     &  pprpts=pprpts,ps1co2=ps1co2,ps1s3=ps1s3,ps2s3=ps2s3,
     &  psecmn=psecmn,psecoc1=psecoc1,psecoc2=psecoc2,rcestr=rcestr,
     &  rictrl=rictrl,riint=riint,rsplig=rsplig,spl=spl,
     &  strmax=strmax,texepp=texepp,texesp=texesp,teff=teff,
     &  tmelt=tmelt,varat11=varat11,varat12=varat12,varat21=varat21,
     &  varat22=varat22,varat3=varat3,vlossg_m=vlossg_m,
     &  sfavail=sfavail)

      call SaveGlobalSiteParams(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  ivauto=ivauto,labtyp=labtyp,labyr=labyr,micosm=micosm,
     &  nelem=nelem,Ninput=Ninput,nlayer=nlayer,Nstart=Nstart,
     &  OMADinput=OMADinput,OMADstart=OMADstart,phsys=phsys,phtm=phtm,
     &  swflag=swflag,claypg_const=claypg_const,tlaypg=tlaypg,
     &  usexdrvrs=usexdrvrs,jdayStart=jdayStart,jdayEnd=jdayEnd,
     &  SnowFlag=SnowFlag,hours_rain=hours_rain,drainlag=drainlag,
     &  watertable=watertable,timlag=timlag,basef=basef,co2ipr=co2ipr,
     &  co2ice=co2ice,co2irs=co2irs,co2itr=co2itr,co2sys=co2sys,
     &  co2tm=co2tm,drain=drain,epnfa=epnfa,epnfs=epnfs,
     &  ckmrspmx=ckmrspmx,cmrspnpp=cmrspnpp,fkmrspmx=fkmrspmx,
     &  fmrsplai=fmrsplai,flodeff=flodeff,ppdf=ppdf,prdx=prdx,
     &  ps2mrsp=ps2mrsp,
     &  precip=precip,pslsrb=pslsrb,rcelit=rcelit,rces1=rces1,
     &  rces2=rces2,rces3=rces3,remwsd=remwsd,rock=rock,
     &  satmos=satmos,snfxmx=snfxmx,sorpmx=sorpmx,stamt=stamt,
     &  stormf=stormf,cmix=cmix,tmix=tmix,ststart=ststart,
     &  stsys=stsys,wscoeff=wscoeff,sublimscale=sublimscale,
     &  reflec=reflec,albedo=albedo,fswcinit=fswcinit,dmpflux=dmpflux,
     &  hpotdeep=hpotdeep,ksatdeep=ksatdeep,cldcov=cldcov, 
     &  tbotmn=tbotmn,tbotmx=tbotmx,dmp=dmp,Ncoeff=Ncoeff,
     &  N2Oadjust_fc=N2Oadjust_fc,N2Oadjust_wp=N2Oadjust_wp,
     &  MaxNitAmt=MaxNitAmt,netmn_to_no3=netmn_to_no3,
     &  wfpsdnitadj=wfpsdnitadj,n2n2oadj=n2n2oadj,
     &  noxn2oadj=noxn2oadj,
     &  elevation=elevation,site_slope=site_slope,aspect=aspect,
     &  ehoriz=ehoriz,whoriz=whoriz,dDO_fc=dDO_fc,dDO_wp=dDO_wp,
     &  rlatitude=rlatitude,sradadj=sradadj,
     &  indewpt=indewpt,
     &  tminintercept=tminintercept,tminslope=tminslope,
     &  maxphoto=maxphoto,bioabsorp=bioabsorp,
     &  mti_mx_incr_r=mti_mx_incr_r,mti_mn_sr=mti_mn_sr,
     &  mti_mx_sr=mti_mx_sr,mdr_mn_redc_r=mdr_mn_redc_r,
     &  mdr_mn_sr=mdr_mn_sr,mdr_mx_sr=mdr_mx_sr,
     &  photo_co2_fraction=photo_co2_fraction,
     &  maxphoto_lig_slp=maxphoto_lig_slp,sitlat=sitlat,sitlng=sitlng,
     &  sitpot=sitpot,sand=sand,silt=silt,clay=clay,sitpot_m=sitpot_m)

      call SaveGlobalSoilProperties(endOfDay=endOfDay, 
     &  dcgridpt=dcgridpt,
     &  numlyrs=numlyrs,nelyrs=nelyrs,ubnd=ubnd,lbnd=lbnd,
     &  width=width,depth=depth,dpthmn=dpthmn,dpthmx=dpthmx,
     &  bulkd=bulkd2,fieldc=fieldc,wiltpt=wiltpt,
     &  ecoeff=ecoeff,tcoeff=tcoeff,sandfrac=sandfrac,
     &  clayfrac=clayfrac,orgfrac=orgfrac,swclimit=swclimit,
     &  satcond=satcond,pH=pH2,swcfc=swcfc,swcwp=swcwp,
     &  swcmin=swcmin,minpot=minpot,sumecoeff=sumecoeff, 
     &  thetas=thetas,thetas_bd=thetas_bd,psis=psis,b=b,
     &  binverse=binverse)

c ... Dynamic variables

      call SaveGlobalChrVars(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  cmdary=cmdary,curcrp=curcrp,curtre=curtre, 
     &  initcp=initcp,initre=initre,wthr=wthr,
     &  curcult=curcult,curfert=curfert,curfire=curfire, 
     &  curgraz=curgraz,curharv=curharv,curirri=curirri,
     &  curomad=curomad,curtrm=curtrm,typary=typary)

c     write(*,*) 'getsetglobalvars: bglcisj(1) =', bglcisj(1)
c     write(*,*) 'getsetglobalvars: bglcism(1) =', bglcism(1)

      call SaveGlobalPlotVars(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  claypg=claypg,texture=texture,aglcis=aglcis,aglive=aglive,
     &  aminrl=aminrl,annet=annet,asmos=asmos,
     &  bglcisj=bglcisj,bglcism=bglcism,
     &  bglivej=bglivej,bglivem=bglivem,cgrain=cgrain,crmvst=crmvst,
     &  egrain=egrain,ermvst=ermvst,co2cce=co2cce,co2crs=co2crs,
     &  co2cpr=co2cpr,co2ctr=co2ctr,crpstg=crpstg,crpval=crpval,
     &  elimit=elimit,metabe=metabe,metcis=metcis,minerl=minerl,
     &  occlud=occlud,parent=parent,secndy=secndy,snlq=snlq,snow=snow,
     &  som1ci=som1ci,som1e=som1e,som2ci=som2ci,som2e=som2e,
     &  som3ci=som3ci,som3e=som3e,stdcis=stdcis,strcis=strcis,
     &  strlig=strlig,struce=struce,stdede=stdede,carbostg=carbostg,
     &  cltfac=cltfac,csrsnk=csrsnk,esrsnk=esrsnk,croote=croote,
     &  crtcis=crtcis,fbrche=fbrche,fbrcis=fbrcis,forstg=forstg,
     &  frootej=frootej,frootem=frootem,frtcisj=frtcisj,
     &  frtcism=frtcism,rleave=rleave,rlvcis=rlvcis,rlwode=rlwode,
     &  rlwcis=rlwcis,w1lig=w1lig,w2lig=w2lig,w3lig=w3lig,
     &  wd1cis=wd1cis,wd2cis=wd2cis,wd3cis=wd3cis,
     &  wood1e=wood1e,wood2e=wood2e,wood3e=wood3e,eftext=eftext,
     &  fps1s3=fps1s3,fps2s3=fps2s3,orglch=orglch,p1co2=p1co2,
     &  ratnew1=ratnew1,ratnew2=ratnew2,cisofr=cisofr,cisotf=cisotf,
     &  nitrate=nitrate,ammonium=ammonium,
     &  frac_nh4_fert=frac_nh4_fert,frac_no3_fert=frac_no3_fert,
     &  vlossg=vlossg,hpttr=hpttr,htran=htran,psloss=psloss,
     &  satmt=satmt,sirri=sirri,trbasl=trbasl,swc=swc,wfps=wfps,
     &  stmtemp=stmtemp,soiltavg=soiltavg,soiltmin=soiltmin,
     &  soiltmax=soiltmax,ptagc=ptagc,ptbgc=ptbgc)

      call sumcar
      call savarp

      call SaveGlobalCropParams(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  auirri=auirri,bioflg=bioflg,crpgrw=crpgrw,flghrv=flghrv,
     &  grzeff=grzeff,himon=himon,omadtyp=omadtyp,seedl=seedl,
     &  frtcindx=frtcindx,curgdys=curgdys,isagri=isagri,
     &  aglivb=aglivb,astrec=astrec,aglrem=aglrem,astgc=astgc,
     &  astlbl=astlbl,astlig=astlig,awhc=awhc,basfc2=basfc2,
     &  bglrem=bglrem,biok5=biok5,biomax=biomax,cfrtcn=cfrtcn,
     &  cfrtcw=cfrtcw,clteff=clteff,cmxturn=cmxturn,crprtf=crprtf,
     &  cultra=cultra,cwlit=cwlit,cwstcr=cwstcr,cwstress=cwstress,
     &  efrgrn=efrgrn,fallrt=fallrt,fawhc=fawhc,fdfrem=fdfrem,
     &  fdgrem=fdgrem,fecf=fecf,feclig=feclig,flfrem=flfrem,
     &  flgrem=flgrem,fligni=fligni,fligni11=fligni11,fnue=fnue,
     &  fret=fret,frtc=frtc,frtsh=frtsh,fsdeth=fsdeth,fulcan=fulcan,
     &  gfcret=gfcret,gret=gret,grwprc=grwprc,gwstress=gwstress,
     &  hibg=hibg,himax=himax,hiwsf=hiwsf,irramt=irramt,
     &  irraut=irraut,mrtfrac=mrtfrac,pltmrf=pltmrf,pltlig=pltlig,
     &  pramn=pramn,pramx=pramx,prbmn=prbmn,prbmx=prbmx,
     &  rdrj=rdrj,rdrm=rdrm,rdsrfc=rdsrfc,rmvstr=rmvstr,rtdtmp=rtdtmp,
     &  sdethc=sdethc,sfclit=sfclit,stdead=stdead,stcrlai=stcrlai,
     &  twhc=twhc,vlossp=vlossp,tmpgerm=tmpgerm,ddbase=ddbase,
     &  tmpkill=tmpkill,mnddhrv=mnddhrv,mxddhrv=mxddhrv,
     &  clsgres=clsgres)

      call SaveGlobalForestParams(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  decid=decid,forgrw=forgrw,basfct=basfct,btolai=btolai,
     &  ccefor=ccefor,cerfor=cerfor,decw1=decw1,decw2=decw2,
     &  decw3=decw3,fcfrac=fcfrac,forrtf=forrtf,klai=klai,
     &  laitop=laitop,ldrmlt=ldrmlt,leafdr=leafdr,maxlai=maxlai,
     &  maxldr=maxldr,maxnp=maxnp,sapk=sapk,swold=swold,
     &  tfrtcn=tfrtcn,tfrtcw=tfrtcw,tree_cfrac=tree_cfrac,
     &  tmxturn=tmxturn,twstress=twstress,wdlig=wdlig,
     &  wmrtfrac=wmrtfrac,woodb=woodb,wooddr=wooddr,
     &  wrdsrfc=wrdsrfc)

c     write(*,*) 'getsetglobalvars: dograz =', dograz
c     write(*,*) 'getsetglobalvars: grazcnt =', grazcnt
c     write(*,*) 'getsetglobalvars: grazday =', grazday

      call SaveGlobalSchedlVars(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  cultday=cultday,erodday=erodday,fertday=fertday,
     &  fireday=fireday,flstday=flstday,foneday=foneday,
     &  frstday=frstday,grazday=grazday,hrvtday=hrvtday,
     &  irriday=irriday,lastday=lastday,omadday=omadday,
     &  plntday=plntday,seneday=seneday,tremday=tremday,
     &  cultcnt=cultcnt,fertcnt=fertcnt,erodcnt=erodcnt,
     &  grazcnt=grazcnt,irricnt=irricnt,plntcnt=plntcnt,
     &  senecnt=senecnt,savefrstday=savefrstday,
     &  saveplntday=saveplntday,evtptr=evtptr,rptyrs=rptyrs,
     &  timary=timary,ttlind=ttlind,cursys=cursys,decsys=decsys,
     &  prevsys=prevsys,docult=docult,doerod=doerod,dofert=dofert,
     &  dofire=dofire,doflst=doflst,dofone=dofone,dofrst=dofrst,
     &  dograz=dograz,dohrvt=dohrvt,doirri=doirri,dolast=dolast,
     &  doomad=doomad,doplnt=doplnt,dosene=dosene,dotrem=dotrem,
     &  frstschd=frstschd,harvschd=harvschd,plntschd=plntschd,
     &  senmschd=senmschd,fltary=fltary)

c ... Remove Nscalar and OMADscalar from global data structure. -mdh 11/11/2022
c     call SaveGlobalFertilVars(endOfDay=endOfDay,dcgridpt=dcgridpt,
c    &  ninhtm=ninhtm,aufert=aufert,feramt=feramt,
c    &  Nscalar=Nscalar,ninhib=ninhib,nreduce=nreduce,
c    &  OMADscalar=OMADscalar,savedfert=savedfert)

      call SaveGlobalFertilVars(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  ninhtm=ninhtm,aufert=aufert,feramt=feramt,
     &  ninhib=ninhib,nreduce=nreduce,
     &  savedfert=savedfert)

      call SaveGlobalPhenoVars(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  cgrwdys=cgrwdys,cpsndys=cpsndys,curgdys=curgdys,
     &  drpdys=drpdys,fgrwdys=fgrwdys,fpsndys=fpsndys,
     &  frstmth=frstmth,furgdys=furgdys,grnfldys=grnfldys,
     &  greenUpCnt=greenUpCnt,numGreenupPeriods=numGreenupPeriods,
     &  plntd=plntd,tfstmth=tfstmth,basetemp=basetemp,
     &  clsgres=clsgres,dayhrs=dayhrs,ddbase=ddbase,
     &  flsgres=flsgres,mnddhrv=mnddhrv,mxddhrv=mxddhrv,
     &  soiltavewk=soiltavewk,thermunits=thermunits,
     &  tmpgerm=tmpgerm,tmpkill=tmpkill,tmplff=tmplff,tmplfs=tmplfs,
     &  accumdd=accumdd,decidgrow=decidgrow,drpdlv=drpdlv,
     &  grnfill=grnfill,grnhrvt=grnhrvt,hrsinc=hrsinc,
     &  plntkill=plntkill,startd=startd,
     &  soiltavgdy=soiltavgdy,isoilwkptr=isoilwkptr)

      call SaveGlobalPhotosynVars(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  aMax=aMax,aMaxFrac=aMaxFrac,
     &  aMaxScalar1=aMaxScalar1,aMaxScalar2=aMaxScalar2,
     &  aMaxScalar3=aMaxScalar3,aMaxScalar4=aMaxScalar4,
     &  attenuation=attenuation,baseFolRespFrac=baseFolRespFrac,
     &  cFracLeaf=cFracLeaf,dVpdExp=dVpdExp,dVpdSlope=dVpdSlope,
     &  growthDays1=growthDays1,growthDays2=growthDays2,
     &  growthDays3=growthDays3,growthDays4=growthDays4,
     &  halfSatPar=halfSatPar,leafCSpWt=leafCSpWt,psnTMin=psnTMin,
     &  psnTOpt=psnTOpt)

c ... Remove daily arrays (dimension NDAY+1) that don't need to be 
c ... saved to the global variable data structure every day.
c ... This includes ppt, srad, sradKJ, parKJ, solrad,
c ... vpd, rhumid, windsp, tempmin, tempmax, avgtemp, 
c ... Internally in DayCent, they will keep their dimension to
c ... avoid massive amounts of code changes. This is safe because
c ... only the array element of "curday" is accessed in the code.
c ... Also remove precscalar(12), tminscalar(12), and tmaxscalar(12).
c ... Add in temperature queues and ptrs for computing tavewk and tavemth.
c ... -mdh 11/16/2022
c     call SaveWthVars(endOfDay=endOfDay,dcgridpt=dcgridpt,
c    &  wthinput=wthinput,wthstart=wthstart,tmn2m=tmn2m,tmx2m=tmx2m,
c    &  maxt=maxt,precscalar=precscalar,tminscalar=tminscalar,
c    &  tmaxscalar=tmaxscalar,srad=srad,sradKJ=sradKJ,parKJ=parKJ,
c    &  vpd=vpd,avgtemp=avgtemp,tempmin=tempmin,tempmax=tempmax,
c    &  ppt=ppt,solrad=solrad,rhumid=rhumid,windsp=windsp)

      call SaveWthVars(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  wthinput=wthinput,wthstart=wthstart,tmn2m=tmn2m,tmx2m=tmx2m,
     &  maxt=maxt,tavewkqueue=tavewkqueue,
     &  tavemthqueue=tavemthqueue,iwkptr=iwkptr,imthptr=imthptr)

      call SaveTimeVars(endOfDay=endOfDay,dcgridpt=dcgridpt,
     &  leapyr=leapyr,month=month,strtyr=strtyr,blktnd=blktnd,
     &  decodt=decodt,dt=dt,dtpl=dtpl,strplt=strplt,tend=tend,
     &  time=time,tplt=tplt)
 
      month = dcinput%month

      if (endOfDay .eq. 0) then

c ..... Set variables in C structures from values stored in global data structure.
        call save_c_vars(
     &    endOfDay,numlyrs,nelyrs,ubnd,lbnd,width,depth,dpthmn,
     &    dpthmx,bulkd2,fieldc,wiltpt,ecoeff,tcoeff,sandfrac,
     &    clayfrac,orgfrac,swclimit,satcond,pH2,swcfc,swcwp,
     &    swc,swcmin,minpot,wfps,sumecoeff,thetas,thetas_bd,
     &    psis,b,binverse,stmtemp,soiltavg,soiltmin,soiltmax,
     &    usexdrvrs,texture,hours_rain,drainlag,jdayStart,jdayEnd,
     &    SnowFlag,watertable,sublimscale,reflec,albedo,fswcinit,
     &    dmpflux,hpotdeep,ksatdeep,rlatitude,cldcov,tbotmn,
     &    tbotmx,dmp,timlag,Ncoeff,N2Oadjust_fc,N2Oadjust_wp,
     &    MaxNitAmt,netmn_to_no3,wfpsdnitadj,n2n2oadj,noxn2oadj,
     &    elevation,
     &    site_slope,aspect,ehoriz,whoriz,dDO_fc,dDO_wp,
     &    daylength,indewpt,soiltavgdy,par,isoilwkptr)

        !At the beginning of each day, insert weather driver values from CWRF into 
        !wthdaily.inc daily arrays at beginning of day

        tempmin(dcinput%curday) = dcinput%tmin
        tempmax(dcinput%curday) = dcinput%tmax
        ppt(dcinput%curday)     = dcinput%ppt
        sradKJ(dcinput%curday)  = dcinput%sradKJ
        parKJ(dcinput%curday)   = dcinput%parKJ
c ..... srad = W/m2, average over daylight period (used by snowmelt and photosynthesis functions)
c ..... srad should be average W/m2 over daylight hours, not over 24 hours. -mdh 12/30/2022
c ....  srad(dcinput%curday)    = dcinput%sradKJ / W2KJ 
        srad(dcinput%curday)    = dcinput%sradKJ * 1000 
     &                            / (daylength(dcinput%curday) * 3600)
c ..... solrad = langleys/day: (langley=cal/cm2; 1000 J/KJ; m2=10000 cm2; cal=4.184J)
        solrad(dcinput%curday)  = dcinput%sradKJ*1000/(10000*4.184)
        rhumid(dcinput%curday)  = dcinput%rhumid
        windsp(dcinput%curday)  = dcinput%windsp
        vpd(dcinput%curday)     = dcinput%vpd
        snow                    = dcinput%snow

c ..... Set variables to allow CRWF CO2 concentration to override fix.100 values
c ....  -mdh 11/9/2022
        co2sys = 2
        co2ppm(1) = dcinput%ambientCO2ppm
        co2ppm(2) = dcinput%ambientCO2ppm

c ....  Don't save weather scalars in global variable data structure. Instead, 
c ....  set them to their default values here (no scaling) since they are
c ....  unlikely to be used, but are still in the DayCent code. -mdh 11/11/2022
        do 30 imo = 1, 12
          precscalar(imo) = 1.0
          tmaxscalar(imo) = 0.0
          tminscalar(imo) = 0.0
          Nscalar(imo) = 1.0
          OMADscalar(imo) = 1.0
          pHscalar(imo) = 1.0
30      continue

      endif 

      if (endOfDay .eq. 1 .and. doinit .eq. 0 
     &    .and. iwrtcsv .eq. 1) then
        open(unit=123, file='dcoutput.csv',status='old',access='append')
        write(123, 124) time,dcinput%curday,month,dcgridpt%scellID,
     &    ptagc,ptbgc,aglivc,bglivcj,rleavc,frootcj,rlwodc,crootc,
     &    metabc(1),metabc(2),strucc(1),strucc(2),
     &    som1c(1),som1c(2),som2c(1),som2c(2),som3c,
     &    agdefac,bgdefac,ammonium,nitrate(1),nitrate(2),nitrate(3),
     &    swc(1)/width(1),swc(2)/width(2),swc(3)/width(3),
     &    aglive(1),bglivej(1),soiltavg(1),
     &    soiltavg(2),soiltavg(3),dcinput%tmin,dcinput%tmax,
     &    dcinput%ppt,dcinput%sradKJ,solrad(dcinput%curday),
     &    dcinput%rhumid,dcinput%windsp,dcinput%vpd,dcinput%snow,
     &    evapdly,trandly,intrcpt,sublim,outflow,petdly,stemp,
     &    carbostg(1,1),tminrl(1),minerl(1,1),minerl(2,1),
     &    minerl(3,1),minerl(4,1),minerl(5,1),minerl(6,1),
     &    minerl(7,1),cgrain,crmvst,egrain(1),ermvst(1)
124     format(f8.2,2(',',i3),',',a6,60(',',f10.4))
        close(123)
      endif

      if (endOfDay .eq. 0) then
        do ilyr = 1, SWMXLYR
          dcgridpt%prevdy%swc(ilyr) = swc(ilyr)
        end do
        dcgridpt%prevdy%snow = snow
        dcgridpt%prevdy%somsc = somsc
        dcgridpt%prevdy%totsysc = totsysc
      endif

      return
      end


c ======================================================================
      integer function getmonth(curday)
        integer curday

        getmonth=-1
        if (curday .ge. 1 .and. curday .le. 31) then
          getmonth = 1
        else if (curday .ge. 32 .and. curday .le. 59) then
          getmonth = 2
        else if (curday .ge. 60 .and. curday .le. 90) then
          getmonth = 3
        else if (curday .ge. 91 .and. curday .le. 120) then
          getmonth = 4
        else if (curday .ge. 121 .and. curday .le. 151) then
          getmonth = 5
        else if (curday .ge. 152 .and. curday .le. 181) then
          getmonth = 6
        else if (curday .ge. 182 .and. curday .le. 212) then
          getmonth = 7
        else if (curday .ge. 213 .and. curday .le. 243) then
          getmonth = 8
        else if (curday .ge. 244 .and. curday .le. 273) then
          getmonth = 9
        else if (curday .ge. 274 .and. curday .le. 304) then
          getmonth = 10
        else if (curday .ge. 305 .and. curday .le. 334) then
          getmonth = 11
        else if (curday .ge. 335 .and. curday .le. 366) then
          getmonth = 12
        endif

        return
        end
