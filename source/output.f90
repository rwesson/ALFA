module mod_output
use mod_functions
use mod_globals
implicit none
contains

subroutine write_plaintext(realspec,fittedspectrum,continuum,skyspectrum,fittedlines,redshiftguess_overall,resolutionguess,normalisation,hbetaflux,outputbasename,totallines)

  type(linelist), dimension(:), allocatable :: fittedlines
  type(spectrum), dimension(:), allocatable :: realspec,fittedspectrum, skyspectrum, continuum
  character(len=512) :: outputbasename
  real :: redshiftguess_overall,resolutionguess
  integer :: totallines
  real :: normalisation, hbetaflux
  integer :: omp_get_thread_num,threadnumber
  integer :: dpfmt
  character(len=1) :: dpfmtch

  threadnumber=omp_get_thread_num()

! set output format sensibly
  dpfmt=max(2,5-floor(log10(fittedspectrum(1)%wavelength)))
  write(dpfmtch,"(I1)") dpfmt

  open(100+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_fit")

  write (100+threadnumber,*) "#alfa version ",VERSION
  write (100+threadnumber,*) "#fit generated using: ",trim(commandline)
  write (100+threadnumber,*) "#""wavelength""  ""input spectrum ""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""sky lines""  ""residuals""  ""uncertainty"""
  do i=1,spectrumlength
    write(100+threadnumber,"(F9."//dpfmtch//", 7(ES12.3))") fittedspectrum(i)%wavelength,realspec(i)%flux + continuum(i)%flux, fittedspectrum(i)%flux + continuum(i)%flux + skyspectrum(i)%flux, realspec(i)%flux, continuum(i)%flux, skyspectrum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux, realspec(i)%uncertainty
  enddo

  close(100+threadnumber)

! now write out the line list.

!  if (maxval(fittedlines%peak) .lt. 0.1 .or. maxval(fittedlines%peak) .gt. 1.e7) then
    fluxformat="ES12.3"
!  else
!    fluxformat="F12.3"
!  endif

  if (messages) print *,gettime(),"writing output files ",trim(outputdirectory),trim(outputbasename),"_lines and ",trim(outputdirectory),trim(outputbasename),"_fit"

  open(200+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_lines")
  write (200+threadnumber,*) "#wlen (obs)  (rest)      flux     uncert.        peak        FWHM   Ion     Multiplet   UpperTerm LowerTerm g1  g2"

!check whether to write out continuum fluxes

  writeb1=0
  writeb2=0
  writep1=0
  writep2=0

  if (minval(continuum%wavelength)-3630.0 .lt.0..and. maxval(continuum%wavelength)-3630. .gt.0.) writeb1=minloc(abs(fittedlines%wavelength-3630.0),1)
  if (minval(continuum%wavelength)-3700.0 .lt.0..and. maxval(continuum%wavelength)-3700. .gt.0.) writeb2=minloc(abs(fittedlines%wavelength-3700.0),1)
  if (minval(continuum%wavelength)-8100.0 .lt.0..and. maxval(continuum%wavelength)-8100. .gt.0.) writep1=minloc(abs(fittedlines%wavelength-8100.0),1)
  if (minval(continuum%wavelength)-8400.0 .lt.0..and. maxval(continuum%wavelength)-8400. .gt.0.) writep2=minloc(abs(fittedlines%wavelength-8400.0),1)

  do i=1,totallines

    if (i.eq.writeb1) then
      write (200+threadnumber,"(2(F9."//dpfmtch//"),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 3630.0, 3630.0, normalisation*continuum(minloc(abs(continuum%wavelength-3630.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
    endif

    if (i.eq.writeb2) then
      write (200+threadnumber,"(2(F9."//dpfmtch//"),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 3700.0, 3700.0, normalisation*continuum(minloc(abs(continuum%wavelength-3700.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty
    endif

    if (i.eq.writep1) then
      write (200+threadnumber,"(2(F9."//dpfmtch//"),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 8100.0, 8100.0, normalisation*continuum(minloc(abs(continuum%wavelength-8100.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty
    endif

    if (i.eq.writep2) then
      write (200+threadnumber,"(2(F9."//dpfmtch//"),"//fluxformat//","//fluxformat//",'         0.0         0.0')") 8400.0, 8400.0, normalisation*continuum(minloc(abs(continuum%wavelength-8400.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty
    endif

    if (fittedlines(i)%blended .eq. 0 .and. fittedlines(i)%uncertainty .gt. detectionlimit) then
      write (200+threadnumber,"(2(F9."//dpfmtch//"),4("//fluxformat//"),3X,A12,A12,A12,A12,I5,I5)") fittedlines(i)%wavelength*fittedlines(i)%redshift, fittedlines(i)%wavelength, gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty,fittedlines(i)%peak/normalisation, fittedlines(i)%wavelength/fittedlines(i)%resolution * 2.35482,fittedlines(i)%ion,fittedlines(i)%multiplet,fittedlines(i)%lowerterm,fittedlines(i)%upperterm,fittedlines(i)%g1,fittedlines(i)%g2
    elseif (fittedlines(i)%blended .ne. 0) then
      if (fittedlines(fittedlines(i)%blended)%uncertainty .gt. detectionlimit) then
        write (200+threadnumber,"(F9."//dpfmtch//",F9."//dpfmtch//",'           *           *           *           *   ',A12,A12,A12,A12,I5,I5)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,fittedlines(i)%ion,fittedlines(i)%multiplet,fittedlines(i)%lowerterm,fittedlines(i)%upperterm,fittedlines(i)%g1,fittedlines(i)%g2
      endif
! write out 3 sigma upper limit for non-detections if upperlimits flag is set
    elseif (fittedlines(i)%uncertainty .le. detectionlimit .and. upperlimits .eqv. .true.) then
      write (200+threadnumber,"(2(F9."//dpfmtch//"),"//fluxformat//",' upper limit',12X,12X,3X,A12,A12,A12,A12,I5,I5)") fittedlines(i)%wavelength*fittedlines(i)%redshift, fittedlines(i)%wavelength, 3.*gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty,fittedlines(i)%ion,fittedlines(i)%multiplet,fittedlines(i)%lowerterm,fittedlines(i)%upperterm,fittedlines(i)%g1,fittedlines(i)%g2
    endif
  enddo

!write out measured Hbeta flux if normalisation was applied, and resolution and velocity

  if (hbetaflux .gt. 0.d0 .and. normalisation .ne. 1.d0) then
    write (200+threadnumber,"(A,ES8.2)") "# measured flux of Hbeta: ",hbetaflux
  endif
  write (200+threadnumber,"(A,ES8.2)") "# spectral resolution: ",resolutionguess
  write (200+threadnumber,"(A,F6.2)") "# radial velocity: ",c*(redshiftguess_overall-1)

!done, close files

  close(200+threadnumber)

end subroutine write_plaintext




subroutine write_latex(realspec,fittedspectrum,continuum,skyspectrum,fittedlines,redshiftguess_overall,resolutionguess,normalisation,hbetaflux,outputbasename,totallines)

  type(linelist), dimension(:), allocatable :: fittedlines
  type(spectrum), dimension(:), allocatable :: realspec,fittedspectrum, skyspectrum, continuum
  character(len=512) :: outputbasename
  real :: redshiftguess_overall,resolutionguess
  integer :: totallines
  real :: normalisation, hbetaflux
  integer :: omp_get_thread_num,threadnumber
  integer :: dpfmt
  character(len=1) :: dpfmtch

  threadnumber=omp_get_thread_num()

! set output format sensibly
  dpfmt=max(2,5-floor(log10(fittedspectrum(i)%wavelength)))
  write(dpfmtch,"(I1)") dpfmt

  open(100+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_fit.tex")

  write (100+threadnumber,*) "%alfa version ",VERSION
  write (100+threadnumber,*) "%fit generated using: ",trim(commandline)
  write (100+threadnumber,*) "#""wavelength""  ""input spectrum ""  ""fitted spectrum""  ""cont-subbed orig"" ""continuum""  ""sky lines""  ""residuals""  ""uncertainty"""
  do i=1,spectrumlength
    write(100+threadnumber,"(F9."//dpfmtch//", 7(ES12.3))") fittedspectrum(i)%wavelength,realspec(i)%flux + continuum(i)%flux, fittedspectrum(i)%flux + continuum(i)%flux + skyspectrum(i)%flux, realspec(i)%flux, continuum(i)%flux, skyspectrum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux, realspec(i)%uncertainty
  enddo

  close(100+threadnumber)

! now write out the line list.

  if (maxval(fittedlines%peak) .lt. 0.1 .or. maxval(fittedlines%peak) .gt. 1.e7) then
    fluxformat="ES12.3"
  else
    fluxformat="F12.3"
  endif

  if (messages) print *,gettime(),"writing output files ",trim(outputdirectory),trim(outputbasename),"_lines.tex and ",trim(outputdirectory),trim(outputbasename),"_fit.tex"

  open(200+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_lines.tex")

  write(200+threadnumber,*) "\\ \hline"
  write(200+threadnumber,*) "Observed wavelength & Rest wavelength & Flux & Uncertainty & Peak & FWHM & Ion & Multiplet & Lower term & Upper term & g$_1$ & g$_2$ \\"
  write(200+threadnumber,*) "\hline"

!check whether to write out continuum fluxes

  writeb1=0
  writeb2=0
  writep1=0
  writep2=0

  if (minval(continuum%wavelength)-3630.0 .lt.0..and. maxval(continuum%wavelength)-3630. .gt.0.) writeb1=minloc(abs(fittedlines%wavelength-3630.0),1)
  if (minval(continuum%wavelength)-3700.0 .lt.0..and. maxval(continuum%wavelength)-3700. .gt.0.) writeb2=minloc(abs(fittedlines%wavelength-3700.0),1)
  if (minval(continuum%wavelength)-8100.0 .lt.0..and. maxval(continuum%wavelength)-8100. .gt.0.) writep1=minloc(abs(fittedlines%wavelength-8100.0),1)
  if (minval(continuum%wavelength)-8400.0 .lt.0..and. maxval(continuum%wavelength)-8400. .gt.0.) writep2=minloc(abs(fittedlines%wavelength-8400.0),1)

  do i=1,totallines

    if (i.eq.writeb1) then
      write (200+threadnumber,"(F9."//dpfmtch//",' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Balmer cont.\\')") 3630.0, normalisation*continuum(minloc(abs(continuum%wavelength-3630.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
    endif

    if (i.eq.writeb2) then
      write (200+threadnumber,"(F9."//dpfmtch//",' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Paschen cont.\\')") 3700.0, normalisation*continuum(minloc(abs(continuum%wavelength-3700.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty
    endif

    if (i.eq.writep1) then
      write (200+threadnumber,"(F9."//dpfmtch//",' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Paschen cont.\\')") 8100.0, normalisation*continuum(minloc(abs(continuum%wavelength-8100.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty
    endif

    if (i.eq.writep2) then
      write (200+threadnumber,"(F9."//dpfmtch//",' &           & ',"//fluxformat//",' & ',"//fluxformat//",' & Brackett cont.\\')") 8400.0, normalisation*continuum(minloc(abs(continuum%wavelength-8400.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty
    endif

    if (fittedlines(i)%blended .eq. 0 .and. fittedlines(i)%uncertainty .gt. detectionlimit) then
      write (200+threadnumber,"(F9."//dpfmtch//",' & ',F9."//dpfmtch//",' & ',"//fluxformat//",' & ',"//fluxformat//",' & ',"//fluxformat//",' & ',"//fluxformat//",A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty, fittedlines(i)%peak/normalisation, fittedlines(i)%wavelength/fittedlines(i)%resolution * 2.35482, fittedlines(i)%ion
    elseif (fittedlines(i)%blended .ne. 0) then
      if (fittedlines(fittedlines(i)%blended)%uncertainty .gt. detectionlimit) then
        write (200+threadnumber,"(F9."//dpfmtch//",' & ',F9."//dpfmtch//",' &            * &            * &            * &            *',A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,fittedlines(i)%ion
      endif
! write out 3 sigma upper limit for non-detections if upperlimits flag is set
    elseif (fittedlines(i)%uncertainty .le. detectionlimit .and. upperlimits .eqv. .true.) then
      write (200+threadnumber,"(F9."//dpfmtch//",' & ',F9."//dpfmtch//",' & ',"//fluxformat//",' & upper limit ',A85)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength, 3.*gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty, fittedlines(i)%ion
    endif
  enddo

!write out measured Hbeta flux to latex table, if normalisation was applied and if output is required

  if (hbetaflux .gt. 0.d0 .and. normalisation .ne. 1.d0) then
    write (200+threadnumber,*) "\hline"
    write (200+threadnumber,"(A,ES8.2,A)") "\multicolumn{10}{l}{Measured flux of H$\beta$: ",hbetaflux,"} \\"
  endif

  write (200+threadnumber,"(A,ES8.2,A)") "\multicolumn{10}{l}{Spectral resolution: ",resolutionguess,"} \\"
  write (200+threadnumber,"(A,F6.2,A)") "\multicolumn{10}{l}{Radial velocity: ",c*(redshiftguess_overall-1),"} \\"

  write (200+threadnumber,*) "\hline"

!done, close files

  close(200+threadnumber)

end subroutine write_latex




subroutine write_csv(realspec,fittedspectrum,continuum,skyspectrum,fittedlines,redshiftguess_overall,resolutionguess,normalisation,hbetaflux,outputbasename,totallines)

  type(linelist), dimension(:), allocatable :: fittedlines
  type(spectrum), dimension(:), allocatable :: realspec,fittedspectrum, skyspectrum, continuum
  character(len=512) :: outputbasename
  real :: redshiftguess_overall,resolutionguess
  integer :: totallines
  real :: normalisation, hbetaflux
  integer :: omp_get_thread_num,threadnumber
  integer :: dpfmt
  character(len=1) :: dpfmtch

  threadnumber=omp_get_thread_num()

! set output format sensibly
  dpfmt=max(2,5-floor(log10(fittedspectrum(i)%wavelength)))
  write(dpfmtch,"(I1)") dpfmt

  open(100+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_fit.csv")

  write (100+threadnumber,*) "#alfa version ",VERSION
  write (100+threadnumber,*) "#fit generated using: ",trim(commandline)
  write (100+threadnumber,*) "#wavelength,input spectrum,fitted spectrum,cont-subbed orig,continuum,sky lines,residuals,uncertainty"
  do i=1,spectrumlength
    write(100+threadnumber,"(F9."//dpfmtch//", ',', 6(ES12.3,','), ES12.3)") fittedspectrum(i)%wavelength,realspec(i)%flux + continuum(i)%flux, fittedspectrum(i)%flux + continuum(i)%flux + skyspectrum(i)%flux, realspec(i)%flux, continuum(i)%flux, skyspectrum(i)%flux, realspec(i)%flux - fittedspectrum(i)%flux, realspec(i)%uncertainty
  enddo

  close(100+threadnumber)

! now write out the line list.

  if (maxval(fittedlines%peak) .lt. 0.1 .or. maxval(fittedlines%peak) .gt. 1.e7) then
    fluxformat="ES12.3"
  else
    fluxformat="F12.3"
  endif

  if (messages) print *,gettime(),"writing output files ",trim(outputdirectory),trim(outputbasename),"_lines.csv and ",trim(outputdirectory),trim(outputbasename),"_fit.csv"

  open(200+threadnumber,file=trim(outputdirectory)//trim(outputbasename)//"_lines.csv")
  write (200+threadnumber,*) "#wlen (obs),(rest),flux,uncertainty,peak,FWHM"

!check whether to write out continuum fluxes

  writeb1=0
  writeb2=0
  writep1=0
  writep2=0

  if (minval(continuum%wavelength)-3630.0 .lt.0..and. maxval(continuum%wavelength)-3630. .gt.0.) writeb1=minloc(abs(fittedlines%wavelength-3630.0),1)
  if (minval(continuum%wavelength)-3700.0 .lt.0..and. maxval(continuum%wavelength)-3700. .gt.0.) writeb2=minloc(abs(fittedlines%wavelength-3700.0),1)
  if (minval(continuum%wavelength)-8100.0 .lt.0..and. maxval(continuum%wavelength)-8100. .gt.0.) writep1=minloc(abs(fittedlines%wavelength-8100.0),1)
  if (minval(continuum%wavelength)-8400.0 .lt.0..and. maxval(continuum%wavelength)-8400. .gt.0.) writep2=minloc(abs(fittedlines%wavelength-8400.0),1)

  do i=1,totallines

    if (i.eq.writeb1) then
      write (200+threadnumber,"(2(F9."//dpfmtch//",','),"//fluxformat//",',',"//fluxformat//",',         0.0,         0.0')") 3630.0, 3630.0, normalisation*continuum(minloc(abs(continuum%wavelength-3630.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty
    endif

    if (i.eq.writeb2) then
      write (200+threadnumber,"(2(F9."//dpfmtch//",','),"//fluxformat//",',',"//fluxformat//",',         0.0,         0.0')") 3700.0, 3700.0, normalisation*continuum(minloc(abs(continuum%wavelength-3700.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty
    endif

    if (i.eq.writep1) then
      write (200+threadnumber,"(2(F9."//dpfmtch//",','),"//fluxformat//",',',"//fluxformat//",',         0.0,         0.0')") 8100.0, 8100.0, normalisation*continuum(minloc(abs(continuum%wavelength-8100.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty
    endif

    if (i.eq.writep2) then
      write (200+threadnumber,"(2(F9."//dpfmtch//",','),"//fluxformat//",',',"//fluxformat//",',         0.0,         0.0')") 8400.0, 8400.0, normalisation*continuum(minloc(abs(continuum%wavelength-8400.)))%flux, normalisation*realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty
    endif

    if (fittedlines(i)%blended .eq. 0 .and. fittedlines(i)%uncertainty .gt. detectionlimit) then
      write (200+threadnumber,"(2(F9."//dpfmtch//",','),3("//fluxformat//",','),"//fluxformat//")") fittedlines(i)%wavelength*fittedlines(i)%redshift, fittedlines(i)%wavelength, gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution)), gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty,fittedlines(i)%peak/normalisation, fittedlines(i)%wavelength/fittedlines(i)%resolution * 2.35482
    elseif (fittedlines(i)%blended .ne. 0) then
      if (fittedlines(fittedlines(i)%blended)%uncertainty .gt. detectionlimit) then
        write (200+threadnumber,"(F9."//dpfmtch//",','F9."//dpfmtch//",',',3(A12,','),A12)") fittedlines(i)%wavelength*fittedlines(i)%redshift,fittedlines(i)%wavelength,"*","*","*","*"
      endif
! write out 3 sigma upper limit for non-detections if upperlimits flag is set
    elseif (fittedlines(i)%uncertainty .le. detectionlimit .and. upperlimits .eqv. .true.) then
      write (200+threadnumber,"(2(F9."//dpfmtch//",','),"//fluxformat//",' upper limit')") fittedlines(i)%wavelength*fittedlines(i)%redshift, fittedlines(i)%wavelength, 3.*gaussianflux(fittedlines(i)%peak,(fittedlines(i)%wavelength/fittedlines(i)%resolution))/fittedlines(i)%uncertainty
    endif
  enddo

!write out measured Hbeta flux if normalisation was applied

  if (hbetaflux .gt. 0.d0 .and. normalisation .ne. 1.d0) then
    write (200+threadnumber,"(A,ES8.2)") "# measured flux of Hbeta: ",hbetaflux
  endif
  write (200+threadnumber,"(A,ES8.2)") "# spectral resolution: ",resolutionguess
  write (200+threadnumber,"(A,F6.2)") "# radial velocity: ",c*(redshiftguess_overall-1)

!done, close files

  close(200+threadnumber)


end subroutine write_csv

subroutine write_fits(realspec,fittedspectrum,continuum,skyspectrum,fittedlines,redshiftguess_overall,resolutionguess,normalisation,hbetaflux,outputbasename,totallines)
!create a single fits file with two extensions, one with the fit and one with the linelist

  implicit none
  integer :: status,unit,readwrite,blocksize,tfields,varidat
  character(len=16) :: extname
  character(len=16),dimension(8) :: ttype_fit,tform_fit,tunit_fit
  character(len=16),dimension(12) :: ttype_lines,tform_lines,tunit_lines
  character(len=16),dimension(4) :: ttype_qc,tform_qc,tunit_qc
  real,dimension(:),allocatable :: linefluxes,linesigmas,linelambdas
  type(linelist), dimension(:), allocatable :: fittedlines
  type(spectrum), dimension(:), allocatable :: realspec,fittedspectrum, skyspectrum, continuum
  character(len=512) :: outputbasename
  real :: redshiftguess_overall,resolutionguess
  integer :: totallines,detectedlines,i,rownumber
  real :: normalisation, hbetaflux
  character(len=8) :: writevalue
  character(len=30) :: cfitsioerror
  integer :: omp_get_thread_num,threadnumber

  threadnumber=omp_get_thread_num()

  status=0
  readwrite=1

! initialise

!  call ftgiou(unit,status)
  unit=51+threadnumber
  call ftinit(unit,"!"//trim(outputdirectory)//trim(outputbasename)//"_fit.fits",blocksize,status)
  call ftphps(unit,16,0,0,status)

! add history keyword to primary HDU

  call ftpdat(unit,status)
  call ftphis(unit,trim(commandline),status)

! first extension for the fit

  ttype_fit=(/"Wavelength      ","InputSpec       ","FittedSpec      ","ContSubbedInput ","Continuum       ","SkyLines        ","Residuals       ","Uncertainty     "/)
  tform_fit=(/"1E","1E","1E","1E","1E","1E","1E","1E"/)
  tunit_fit=(/"Angstrom        ","Flux            ","Flux            ","Flux            ","Flux            ","Flux            ","Flux            ","Flux            "/)

  extname='FIT'
  varidat=0
  tfields=8

  call ftibin(unit,spectrumlength,tfields,ttype_fit,tform_fit,tunit_fit,extname,varidat,status)

! record version, date and setting used in comments

  call ftpcom(unit,"Produced by alfa version "//VERSION,status)
  call ftpcom(unit,"Command line: '"//trim(commandline)//"'",status)
  call ftpcom(unit,"input file: "//trim(spectrumfile),status)
!  if (len(trim(imagesection)).gt.0) print *,"                fitting section:               ",imagesection
  if (.not.normalise) then
  call ftpcom(unit,"normalisation:                    using measured value of Hb",status)
  else
    if (normalisation.eq.0.d0) then
  call ftpcom(unit,"normalisation:                    no normalisation",status)
    else
  write(writevalue,"(F8.3)") normalisation
  call ftpcom(unit,"normalisation:                    to Hb="//writevalue,status)
    endif
  endif
  if (subtractcontinuum) then
  call ftpcom(unit,"continuum fitting:                enabled",status)
  write(writevalue,"(I8)") continuumwindow
  call ftpcom(unit,"continuum window:                 "//writevalue,status)
  else
  call ftpcom(unit,"continuum fitting:                disabled",status)
  endif
  write(writevalue,"(F8.2)") baddata
  call ftpcom(unit,"spectrum fitted if max value >    "//writevalue,status)
  write(writevalue,"(F8.2)") wavelengthscaling
  call ftpcom(unit,"Angstroms per wavelength unit:    "//writevalue,status)
  if (tablewavelengthcolumn.ne.1) then
  write(writevalue,"(I2)") tablewavelengthcolumn
  call ftpcom(unit,"table wavelength column:          "//writevalue,status)
  endif
  if (tablefluxcolumn.ne.2) then
  write(writevalue,"(I2)") tablefluxcolumn
  call ftpcom(unit,"table flux column:                "//writevalue,status)
  endif
  if (collapse) then
  call ftpcom(unit,"multiple spectra:                  collapsed to 1D",status)
  else
  call ftpcom(unit,"multiple spectra:                  fitted individually",status)
  endif
!  call ftpcom(unit,"velocity guess:                   ",status)
!  call ftpcom(unit,redshiftguess_initial,status)
!  call ftpcom(unit,"resolution guess:                 ",status)
!  call ftpcom(unit,resolutionguess_initial,status)
  write(writevalue,"(F8.2)") vtol1*c
  call ftpcom(unit,"first pass velocity tolerance:    "//writevalue,status)
  write(writevalue,"(F8.2)") vtol2*c
  call ftpcom(unit,"second pass velocity tolerance:   "//writevalue,status)
  write(writevalue,"(F8.2)") rtol1
  call ftpcom(unit,"first pass resolution tolerance:  "//writevalue,status)
  write(writevalue,"(F8.2)") rtol2
  call ftpcom(unit,"second pass resolution tolerance: "//writevalue,status)
  if (subtractsky) then
  call ftpcom(unit,"sky line fitting:                 enabled",status)
  call ftpcom(unit,"sky line catalogue:               "//trim(skylinelistfile),status)
  else
  call ftpcom(unit,"sky line fitting:                 disabled",status)
  endif
  call ftpcom(unit,"strong line catalogue:            "//trim(stronglinelistfile),status)
  call ftpcom(unit,"deep line catalogue:              "//trim(deeplinelistfile),status)
!  if (exclusioncount .gt. 0) then
!  call ftpcom(unit,"lines excluded from fitting:      "//exclusions,status)
!  endif
  if (rebinfactor .gt. 1) then
  write(writevalue,"(I8)") rebinfactor
  call ftpcom(unit,"spectra rebinned by factor of:    "//writevalue,status)
  endif
  write(writevalue,"(I8)") generations
  call ftpcom(unit,"number of generations:            "//writevalue,status)
  write(writevalue,"(I8)") popsize
  call ftpcom(unit,"population size:                  "//writevalue,status)
  write(writevalue,"(F8.2)") pressure
  call ftpcom(unit,"pressure factor:                  "//writevalue,status)
  call ftpcom(unit,"output directory:                 "//trim(outputdirectory),status)
  call ftpcom(unit,"output format:                    "//outputformat,status)

  call ftpcle(unit,1,1,1,spectrumlength,fittedspectrum%wavelength,status)
  call ftpcle(unit,2,1,1,spectrumlength,realspec%flux + continuum%flux,status)
  call ftpcle(unit,3,1,1,spectrumlength,fittedspectrum%flux + continuum%flux + skyspectrum%flux,status)
  call ftpcle(unit,4,1,1,spectrumlength,realspec%flux,status)
  call ftpcle(unit,5,1,1,spectrumlength,continuum%flux,status)
  call ftpcle(unit,6,1,1,spectrumlength,skyspectrum%flux,status)
  call ftpcle(unit,7,1,1,spectrumlength,realspec%flux - fittedspectrum%flux,status)
  call ftpcle(unit,8,1,1,spectrumlength,realspec%uncertainty,status)

  if (status .gt. 0) then
    call ftgerr(status,cfitsioerror)
    print *,gettime(),"[107] CFITSIO error: ",status,cfitsioerror
    print *,gettime(),"thread ",threadnumber,", unit ",unit
    call exit(107)
  else
    print *,gettime(),"Wrote fit to extension FIT of output file ",trim(outputdirectory)//trim(outputbasename)//"_fit.fits"
  endif

! second extension for the lines

  status=0
  extname="LINES"
  tfields=12

  ttype_lines=(/"WlenObserved    ","WlenRest        ","Flux            ","Uncertainty     ","Peak            ","FWHM            ","Ion             ","Multiplet       ","LowerTerm       ","Upperterm       ","g1              ","g2              "/)
  tform_lines=(/"1E ","1E ","1E ","1E ","1E ","1E ","12A","12A","12A","12A","1I ","1I "/)
  tunit_lines=(/"Angstrom        ","Angstrom        ","Flux            ","Flux            ","Flux            ","Angstrom        ","                ","                ","                ","                ","                ","                "/)

  call ftibin(unit,totallines,tfields,ttype_lines,tform_lines,tunit_lines,extname,varidat,status)

  allocate(linefluxes(totallines))
  allocate(linesigmas(totallines))
  allocate(linelambdas(totallines))

! replicates gaussianflux function. why can't that be vectorised?
! should replace this with calculation at time of fitting?
  linefluxes=fittedlines%peak*(fittedlines%wavelength/fittedlines%resolution)*(2*3.14159265359)**0.5
  linesigmas=linefluxes/fittedlines%uncertainty
  linelambdas=fittedlines%wavelength*fittedlines%redshift

! deal with blends and upper limits

  detectedlines=totallines
  do i=1,totallines
    if (fittedlines(i)%blended .ne. 0) then
        linelambdas(i)=0
      if (fittedlines(fittedlines(i)%blended)%uncertainty .gt. detectionlimit) then
        linefluxes(i)=0
        linesigmas(i)=0
      else ! blends of non-detections. todo: save with null value
        linelambdas(i)=0
        linefluxes(i)=-999.
        linesigmas(i)=-999.
        detectedlines=detectedlines-1
      endif
! write out 3 sigma upper limit as a negative flux for non-detections
    elseif (fittedlines(i)%uncertainty .le. detectionlimit) then
      linefluxes(i)=-detectionlimit*linesigmas(i)
      linesigmas(i)=-linesigmas(i)
      detectedlines=detectedlines-1
    endif
  enddo

  call ftpcle(unit,1,1,1,totallines,linelambdas,status)
  call ftpcle(unit,2,1,1,totallines,fittedlines%wavelength,status)
  call ftpcle(unit,3,1,1,totallines,linefluxes,status)
  call ftpcle(unit,4,1,1,totallines,linesigmas,status)
  call ftpcle(unit,5,1,1,totallines,fittedlines%peak/normalisation,status)
  call ftpcle(unit,6,1,1,totallines,fittedlines%wavelength/fittedlines%resolution * 2.35482,status)
  call ftpcls(unit,7,1,1,totallines,fittedlines%ion,status)
  call ftpcls(unit,8,1,1,totallines,fittedlines%multiplet,status)
  call ftpcls(unit,9,1,1,totallines,fittedlines%lowerterm,status)
  call ftpcls(unit,10,1,1,totallines,fittedlines%upperterm,status)
  call ftpclj(unit,11,1,1,totallines,fittedlines%g1,status)
  call ftpclj(unit,12,1,1,totallines,fittedlines%g2,status)

!write out continuum fluxes if measured
  rownumber=totallines

  if (minval(continuum%wavelength)-3630.0 .lt.0..and. maxval(continuum%wavelength)-3630. .gt.0.) then
    call ftirow(unit,rownumber,1,status)
    rownumber=rownumber+1
    call ftpcle(unit,1,rownumber,1,1,(/3630.0/),status)
    call ftpcle(unit,2,rownumber,1,1,(/3630.0/),status)
    call ftpcle(unit,3,rownumber,1,1,normalisation*continuum(minloc(abs(continuum%wavelength-3630.)))%flux,status)
    call ftpcle(unit,4,rownumber,1,1,normalisation*realspec(minloc(abs(continuum%wavelength-3630.)))%uncertainty,status)
  endif

  if (minval(continuum%wavelength)-3700.0 .lt.0..and. maxval(continuum%wavelength)-3700. .gt.0.) then
    call ftirow(unit,rownumber,1,status)
    rownumber=rownumber+1
    call ftpcle(unit,1,rownumber,1,1,(/3700.0/),status)
    call ftpcle(unit,2,rownumber,1,1,(/3700.0/),status)
    call ftpcle(unit,3,rownumber,1,1,normalisation*continuum(minloc(abs(continuum%wavelength-3700.)))%flux,status)
    call ftpcle(unit,4,rownumber,1,1,normalisation*realspec(minloc(abs(continuum%wavelength-3700.)))%uncertainty,status)
  endif

  if (minval(continuum%wavelength)-8100.0 .lt.0..and. maxval(continuum%wavelength)-8100. .gt.0.) then
    call ftirow(unit,rownumber,1,status)
    rownumber=rownumber+1
    call ftpcle(unit,1,rownumber,1,1,(/8100.0/),status)
    call ftpcle(unit,2,rownumber,1,1,(/8100.0/),status)
    call ftpcle(unit,3,rownumber,1,1,normalisation*continuum(minloc(abs(continuum%wavelength-8100.)))%flux,status)
    call ftpcle(unit,4,rownumber,1,1,normalisation*realspec(minloc(abs(continuum%wavelength-8100.)))%uncertainty,status)
  endif

  if (minval(continuum%wavelength)-8400.0 .lt.0..and. maxval(continuum%wavelength)-8400. .gt.0.) then
    call ftirow(unit,rownumber,1,status)
    rownumber=rownumber+1
    call ftpcle(unit,1,rownumber,1,1,(/8400.0/),status)
    call ftpcle(unit,2,rownumber,1,1,(/8400.0/),status)
    call ftpcle(unit,3,rownumber,1,1,normalisation*continuum(minloc(abs(continuum%wavelength-8400.)))%flux,status)
    call ftpcle(unit,4,rownumber,1,1,normalisation*realspec(minloc(abs(continuum%wavelength-8400.)))%uncertainty,status)
  endif

  if (status .gt. 0) then
    call ftgerr(status,cfitsioerror)
    print *,gettime(),"[107] CFITSIO error: ",status,cfitsioerror
    call exit(107)
  else
    print *,gettime(),"Wrote line list to extension LINES"
  endif

! third extension for nlines, f(hb), rv, resolution

  status=0
  extname="QC"
  tfields=4

  ttype_qc=(/"NumberOfLines   ","HbetaFlux       ","RadialVelocity  ","Resolution      "/)
  tform_qc=(/"1J","1E","1E","1E"/)
  tunit_qc=(/"                ","SameAsInputFlux ","kms-1           ","                "/)

  call ftibin(unit,1,tfields,ttype_qc,tform_qc,tunit_qc,extname,varidat,status)

  call ftpclj(unit,1,1,1,1,(/detectedlines/),status)
  call ftpcle(unit,2,1,1,1,(/hbetaflux/),status)
  call ftpcle(unit,3,1,1,1,(/c*(redshiftguess_overall-1)/),status)
  call ftpcle(unit,4,1,1,1,(/resolutionguess/),status)

  if (status .gt. 0) then
    call ftgerr(status,cfitsioerror)
    print *,gettime(),"[107] CFITSIO error: ",status,cfitsioerror
    call exit(107)
  else
    print *,gettime(),"Wrote QC info to extension QC"
  endif

! close

  call ftclos(unit, status)
  call ftfiou(unit, status)

end subroutine write_fits

end module mod_output
