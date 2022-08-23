program norminv2
!
!*****************************************************************************80
! Program NorMinV2 accepts partial extraction chemical analyses and
! calculates a corresponding normative copper-iron sulfide content. This
! program was written in FORTRAN 90 for Minera Escondida Ltda. using the
! GNU Fortran compiler (MinGW.org GCC-8.2.0-3.
! User documentation is provided in the accompanying report 
!   < NorMinV2 Programming and User Guide >
!
! Copyright (C) 2020  Richard K. Preece
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details at <https://www.gnu.org/licenses/>.
!
! File:   NorMinV2.f90
! Author: Richard Preece
!
! VERSION 1 Created August 16, 2019, NORMINESC to calculate 3-component 
! sulfide assemblages. Two different mineral models were available for
! calculation of CP-CC-CV and CP-BN-(CC-CV) mineral models. The program
! reproduced algorithm, capabilities and options previously coded as 
! Microsoft EXCEL spreadsheets. Additional capabilities include
!     - user selection of sulfide mineral model (with or without bornite)
!     - capability to read variably-formatted csv input files
!     - capability to write csv output file
!     - calculation of volume percent mineral content
! VERSION 2 Created February 2020, NORMINV2 to add capabilities for calculating
! 4-component oxide-sulfide mineral models. Program automatically selects mineral
! model based on the Minera Escondida Mineral Zone coding structure and certain 
! analytical criteria. 
! VERSION 2.1 Modified February 2021, Corrected bug in get_real subroutine that
! failed to read input data when in the last value field of the line   
! VERSION 2.2 Modified June 2021, Revised MinZone codes for Leach Cap and Oxide 
! for consistency with standard MEL geological codes
!           Leach Cap: Old Code= 1  New Code= 0
!           Oxide:     Old Code= 2  New Code= 1    
! VERSION 2.3 Modified August 2022, Added TCu, Fe and S2 values to output file    
!
! Define variable types and values
!
    integer err_flag
    real (kind = 8), DIMENSION(3,3) :: m_cpcccv, m_cpccbn, m_cpcvbn
    real (kind = 8), dimension(4,4) :: m_goccbrxc, m_cpccbrxc, m_cpcccvbr
    real (kind = 8), dimension(4,4) :: m_min1, m_min2, m_min3, m_min4
    real (kind = 8) :: d_cpcccv, d_cpccbn, d_cpcvbn
    real (kind = 8) :: d_goccbrxc, d_cpccbrxc, d_cpcccvbr
    real (kind = 8) :: d_min1, d_min2, d_min3, d_min4
    real (kind = 8) :: m33det, m44det

    character (len = 25) :: paramfile, inputcsv, outputcsv
    character (len = 25) :: dhname, sampid
    character (len = 1) :: yesno
    character (len = 80) :: dumline
    character (len = 256) :: record, rec_out

! Current program version number 23 Aug 2022
    vers = 2.3
    idebug = 0        ! Debug switch. Set = 1 to print raw & final CSR values

!==============================================================================
!    START DATA INPUT SECTION
!    To be replaced by direct IO to data base
!==============================================================================
! Welcome message to screen
    call timestamp()
    write ( *, '(a)') ' '
    write ( *, '(a, f5.1,/)') &
         'ESCONDIDA NORMATIVE MINERAL CALCULATIONS version:', vers

! Ask for parameter file name, stop program if not found
    write ( *, '(a)') '  Enter name of parameter file'
    read (*, '(a)') paramfile
    open (52, file = paramfile, status = "old", iostat = ios)
    if (ios > 0) then
        write (*, '(2a)') paramfile, ' is not found'
        stop
    end if
    
! Read parameter file for mineral definition parameters
    write (*, '(a)') '    Reading mineral data from parameter file'
    read (52, '(a)') dumline
    read (52, 10) cpxcu,cpxfe,cpxs,cpxas,cprfs,cprcn,cpras
    read (52, 10) ccxcu,ccxfe,ccxs,ccxas,ccrfs,ccrcn,ccras
    read (52, 10) cvxcu,cvxfe,cvxs,cvxas,cvrfs,cvrcn,cvras
    read (52, 10) bnxcu,bnxfe,bnxs,bnxas,bnrfs,bnrcn,bnras
    read (52, 10) enxcu,enxfe,enxs,enxas,enrfs,enrcn,enras
    read (52, 10) brxcu,brxfe,brxs,brxas,brrfs,brrcn,brras
    read (52, 10) xcxcu,xcxfe,xcxs,xcxas,xcrfs,xcrcn,xcras
    read (52, 10) goxcu,goxfe,goxs,goxas,gorfs,gorcn,goras
10  format (15x,7f10.3)

! Reading data file structure defined in parameter file
    write (*, '(a)') '      Now reading data file structure'
! Skip header lines
    do i = 1, 3
        read (52, '(a)') dumline
    end do
    read (52, 11) nhdr
    read (52, 11) idh
    read (52, 11) ismp
    read (52, 11) ifrom
    read (52, 11) ito
    read (52, 11) itcu
    read (52, 11) iscu
    read (52, 11) icncu
    read (52, 11) ifscu
    read (52, 11) ife
    read (52, 11) is2
    read (52, 11) imz
    read (52, 12) ias, ind_as
    read (52, 12) ibn, idef_bn
    read (52, 13) idens, dens_def
11 format (15x,i10)
12 format (15x, 2i10)
13 format (15x, i10,f10.0)
    close (52)

! Ask for data input file name, stop program if not found
    write ( *, '(a)') '  Enter name of input data file'
    read (*, '(a)') inputcsv
    open (51, file = inputcsv, status = "old", iostat = ios)
    if (ios > 0) then
        write (*, '(2a)') inputcsv, ' is not found'
        stop
    end if

! Ask for output file name, check if file currently exists
100 write ( *, '(a)') '  Enter name of output data file'
    read (*, '(a)') outputcsv
    open (61, file = outputcsv, status = "new", iostat = ios, &
              err = 101)

! Error recovery if data file currently exists
101 continue
    if (ios > 0) then
        write (*, '(2a)') outputcsv, ' already exists.'
        write (*, fmt = '(a)', advance = 'NO') &
        'Enter Y to overwrite current file, N to enter new file name >'
        read (*,'(a1)') yesno
        if (yesno == 'Y' .or. yesno == 'y') then
            open (61, file = outputcsv, status = "unknown", &
              iostat = ios, err = 100)
        else
            go to 100
        end if
    end if
    if (idebug .eq. 1) open (62, file = "debuglist.dat", status = "unknown")

! Write header record for output file
!       Version 2.3, Added TCU, Fe, S2 fields to header record    
    write (61,'(6a)') &
     'dhname,f_int,t_int,sampid,err_flag,minzon,ind_bn,tcu,fe,s2,', &
     'csrcp,csrcc,csrcv,', &
     'csrbn,csrgo,csrbr,csrxc,cspcp,cspcc,cspcv,cspbn,cspgo,cspbr,cspxc,', &
     'cspen,cp_wtpct,cc_wtpct,cv_wtpct,bn_wtpct,go_wtpct,br_wtpct,xc_wtpct,', &
     'en_wtpct,py_best,xfe,cp_volpct,cc_volpct,cv_volpct,bn_volpct,', &
     'go_volpct,br_volpct,xc_volpct,en_volpct,py_volpct'

!==============================================================================
!     END OF SECTION FOR SETTING UP DATA INPUT FILES
!      Begin section for initializing variables and reading data records
!==============================================================================
     ! Initialize record counter and move to first record of input file
    ndata = 0
    do i = 1, nhdr
        read (51,'(a80)') dumline
    end do
!
! Set matrix for end-member extraction rates of cp-cc-cv, calculate determinant
    do j = 1, 3
        m_cpcccv (1,j) = 1.0
    end do
    m_cpcccv (2,1) = cprcn
    m_cpcccv (2,2) = ccrcn
    m_cpcccv (2,3) = cvrcn
    m_cpcccv (3,1) = cprfs
    m_cpcccv (3,2) = ccrfs
    m_cpcccv (3,3) = cvrfs
    d_cpcccv = M33DET (m_cpcccv)

! Set matrix for end-member extraction rates of cp-cc-bn, calculate determinant
    do j = 1, 3
        m_cpccbn (1,j) = 1.0
    end do
    m_cpccbn (2,1) = cprcn
    m_cpccbn (2,2) = ccrcn
    m_cpccbn (2,3) = bnrcn
    m_cpccbn (3,1) = cprfs
    m_cpccbn (3,2) = ccrfs
    m_cpccbn (3,3) = bnrfs
    d_cpccbn = M33DET (m_cpccbn)

! Set matrix for end-member extraction rates of cp-cv-bn, calculate determinant
    do j = 1, 3
        m_cpcvbn (1,j) = 1.0
    end do
    m_cpcvbn (2,1) = cprcn
    m_cpcvbn (2,2) = cvrcn
    m_cpcvbn (2,3) = bnrcn
    m_cpcvbn (3,1) = cprfs
    m_cpcvbn (3,2) = cvrfs
    m_cpcvbn (3,3) = bnrfs
    d_cpcvbn = M33DET (m_cpcvbn)

! Set matrix for end-member extraction rates of CuGoe-cc-br-chrys, calculate determinant
    do j = 1, 4
        m_goccbrxc (1,j) = 1.0
    end do
    m_goccbrxc (2,1) = gorcn
    m_goccbrxc (2,2) = ccrcn
    m_goccbrxc (2,3) = brrcn
    m_goccbrxc (2,4) = xcrcn
    m_goccbrxc (3,1) = gorfs
    m_goccbrxc (3,2) = ccrfs
    m_goccbrxc (3,3) = brrfs
    m_goccbrxc (3,4) = xcrfs
    m_goccbrxc (4,1) = goras
    m_goccbrxc (4,2) = ccras
    m_goccbrxc (4,3) = brras
    m_goccbrxc (4,4) = xcras
    d_goccbrxc = M44DET (m_goccbrxc)

! Set matrix for end-member extraction rates of cp-cc-br-chrys, calculate determinant
    do j = 1, 4
        m_cpccbrxc (1,j) = 1.0
    end do
    m_cpccbrxc (2,1) = cprcn
    m_cpccbrxc (2,2) = ccrcn
    m_cpccbrxc (2,3) = brrcn
    m_cpccbrxc (2,4) = xcrcn
    m_cpccbrxc (3,1) = cprfs
    m_cpccbrxc (3,2) = ccrfs
    m_cpccbrxc (3,3) = brrfs
    m_cpccbrxc (3,4) = xcrfs
    m_cpccbrxc (4,1) = cpras
    m_cpccbrxc (4,2) = ccras
    m_cpccbrxc (4,3) = brras
    m_cpccbrxc (4,4) = xcras
    d_cpccbrxc = M44DET (m_cpccbrxc)

! Set matrix for end-member extraction rates of cp-cc-cv-br, calculate determinant
    do j = 1, 4
        m_cpcccvbr (1,j) = 1.0
    end do
    m_cpcccvbr (2,1) = cprcn
    m_cpcccvbr (2,2) = ccrcn
    m_cpcccvbr (2,3) = cvrcn
    m_cpcccvbr (2,4) = brrcn
    m_cpcccvbr (3,1) = cprfs
    m_cpcccvbr (3,2) = ccrfs
    m_cpcccvbr (3,3) = cvrfs
    m_cpcccvbr (3,4) = brrfs
    m_cpcccvbr (4,1) = cpras
    m_cpcccvbr (4,2) = ccras
    m_cpcccvbr (4,3) = cvras
    m_cpcccvbr (4,4) = brras
    d_cpcccvbr = M44DET (m_cpcccvbr)
!
!------------------------------------------------------------------
!    Program reads and processes one data record at a time
!    SECTION FOR READING DATA RECORD FROM CVS FILE
!      Replace with code for reading data base record
!------------------------------------------------------------------
! Read new record
200 ndata = ndata + 1
    read (51, 201, end = 900) record
201 format(a)

! Parse record to retrieve data
! Find number of values (n_count) and record length (l_count)
    call csv_value_count(record, istatus, n_count, l_count)
! Call subroutine to extract data, identified by counting number of values
    ! Different subroutines used based on data type (real, integer, character)
    call get_char (dhname,idh,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in DH Name on line ', ndata,'Error code = ',ios
    call get_char (sampid,ismp,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in Sample ID on line ', ndata,'Error code = ',ios
    call get_real (f_int,ifrom,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in From on line ', ndata,'Error code = ',ios
    call get_real (t_int,ito,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in From on line ', ndata,'Error code = ',ios
    call get_real (tcu,itcu,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in TCu on line ', ndata,'Error code = ',ios
    call get_real (scu,iscu,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in SCu on line ', ndata,'Error code = ',ios
    call get_real (cncu,icncu,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in CNCu on line ', ndata,'Error code = ',ios
    call get_real (fscu,ifscu,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in FSCu on line ', ndata,'Error code = ',ios
    call get_real (fe,ife,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in Fe on line ', ndata,'Error code = ',ios
    call get_real (s2,is2,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in S2 on line ', ndata,'Error code = ',ios
    call get_real (as,ias,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in As on line ', ndata,'Error code = ',ios
    call get_integer (minzon,imz,record,n_count,l_count,ios)
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in MinZone on line ', ndata,'Error code = ',ios
    call get_integer (ind_bn,ibn,record,n_count,l_count,ios)
    if (ind_bn .lt. 0) ind_bn = idef_bn                     ! Set default flag value
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in Bornite Flag on line ', ndata,'Error code = ',ios
    call get_real (density,idens,record,n_count,l_count,ios)
    if (density .lt. 0) density = dens_def         ! Set default density value
    if (ios .gt. 0) write (*,'(2(a,i6))') &
        'Error in Density on line ', ndata,'Error code = ',ios
!-----------------------------------------------------------------------------
!   END SECTION FOR READING DATA RECORD
!        Begin section for processing PtXt data and calculation of 
!         normative mineralogy
!------------------------------------------------------------------------------
! Check for complete assay set
    if (min(tcu,scu,cncu,fscu).lt. 0.0) then
        write (*,'(a,i6,/,a)') 'Assay data missing in line number ',ndata, &
                               ' Record not processed '
        go to 200
    end if
! Check for valid mineral zone
    if (minzon .lt. 0 .or. minzon .gt. 10) then
        write (*,'(a,i6,/,a)') 'Invalid mineral zone in line number ',ndata, &
                               ' Record not processed '
        go to 200
    end if

    if(idebug .eq. 1) &
    write(62,'(i6,1x,a,2f10.1,i5)') ndata,dhname,f_int,t_int,minzon
! *****************************************************************************
! Start normative mineral calculations
! Six potential mineral associations are defined, 3 for sulfide and 3 for oxide/mixed
! Selection of the mineral association depends on: 
!    - logged mineral zone 
!         MinZones 0, 1, 5 routed to oxide/mixed 4-component associations
!         Minzones 4, 6, 7, 8, 10 routed to sulfide 3-component associations
!    - analytical criteria
!         Sulfide sulfur <0.15% S or SCu >= CNCu routed to oxide only association
!    - normative mineral errors
!         Normative mineral CSR <-0.1 or > 1.1 routed to alternative mineral association
! *****************************************************************************
    if (minzon .le. 2 .or. minzon .eq. 5) go to 400
    
! -----------------------------------------------------------------------------
! Start Sulfide association
! -----------------------------------------------------------------------------
! Adjust assays for As-bearing mineral, depending on ind_as flag setting
! Skip adjustment in oxide mineral zones
    if (as .lt. 0 .or. ind_as .eq. 0) then   ! ind_as = 0 don't process As values
        as = 0.0
        en_wtpct = 0.0
        go to 210
    else if (ind_as == 1) then        ! ind_as = 1 value input as ppm As, convert to wt. %
        as = as / 10000
    end if                            ! ind_as = 2 value as wt. %

! Calculate wt % arsenic mineral from As assay
! Adjust As mineral content for those samples with low TCu and elevated As
    en_wtpct = min(as/enxas,(tcu-0.005)/enxcu) ! Check for samples with
210 cspen = en_wtpct*enxcu                     ! Calculate CSPen from mineral and Cu content
!
! Adjust PtXt for contribution by As mineral
    tcumod = tcu - cspen
    fscumod = max(0.001,fscu - (cspen*enrfs))
    cncumod = max(0.001,cncu - (cspen*enrcn))
    fscu_rat = fscumod/tcumod
    cncu_rat = cncumod/tcumod

! Calculation of chalcopyrite-chalcocite-covellite association, executed for
! samples with ind_bn = 0
300 continue
    if (ind_bn == 1) go to 320
    csrcp_raw =((ccrcn*cvrfs)-(ccrfs*cvrcn)+(cvrcn-ccrcn)*fscu_rat   &
                +(ccrfs-cvrfs)*cncu_rat)/d_cpcccv
    csrcc_raw =((cvrcn*cprfs)-(cvrfs*cprcn)+(cvrfs-cprfs)*cncu_rat   &
                +(cprcn-cvrcn)*fscu_rat)/d_cpcccv
    csrcv_raw = 1-csrcp_raw-csrcc_raw
    csrbn_raw = 0.0

! Set flags for calculated csr value that exceed tolerances of 10%
    if (max(csrcp_raw,csrcc_raw,csrcv_raw).le. 1.100 .and. &
      min(csrcp_raw,csrcc_raw,csrcv_raw).ge. -0.100) then
        err_flag = 0
        go to 310
    end if
    if (tcu < 0.15 .and. min(csrcp_raw,csrcc_raw,csrcv_raw)< -0.100) then
        err_flag = 1         ! Low Grade flag, likely round-off errors or high As value
        if (csrcv_raw < -0.100) then
            if(idebug .eq. 1) &
        write(62,'(a,3f7.3)') ' cpcccv raw ->',csrcp_raw,csrcc_raw,csrcv_raw
            go to 400
        end if
    else if (csrcv_raw < -0.100) then
        err_flag = 2         ! Elevated FSCu relative to CNCu, potentially oxidized
        if(idebug .eq. 1) &
        write(62,'(a,3f7.3)') ' cpcccv raw ->',csrcp_raw,csrcc_raw,csrcv_raw
        go to 400
    else if (csrcp_raw < -0.100) then
        err_flag = 3         ! Elevated CNCu relative to FSCu
    else
        err_flag = 4         ! Undefined error
    end if
!
! Adjust mineral csr for csrcp < 0 and csrcv < 0
310 continue
    if (csrcp_raw .lt. 0.000) then
        csrcp1 = 0.0000
    else if (csrcv_raw .lt. 0.000) then
        csrcp1 = csrcp_raw/(csrcp_raw + csrcc_raw)
    else
        csrcp1 = csrcp_raw
    end if

    if (csrcv_raw .lt. 0.000) then
        csrcc1 = csrcc_raw/(csrcp_raw + csrcc_raw)
        csrcv1 = 0.0000
    else
        csrcc1 = csrcc_raw
        csrcv1 = csrcv_raw
    end if
!
! Second-stage adjustment, for csrcp1 > 1 and csrcc1 < 0
    if (csrcp1 .gt. 1.000) then
        csrcp = csrcp1/(csrcp1 + csrcv1)
    else
        csrcp = csrcp1
    end if

    if (csrcc1 .lt. 0.000) then
        csrcc = 0.0000
    else
        csrcc = csrcc1/(csrcp1 + csrcc1 + csrcv1)
    end if

    csrcv = 1 - csrcp - csrcc
    csrbn = 0.0000
    csrgo = 0.0000
    csrbr = 0.0000
    csrxc = 0.0000
    
    if(idebug .eq. 1) &
    write(62,'(a,3f7.3,a,3f7.3)') ' cpcccv raw ->',csrcp_raw,csrcc_raw,csrcv_raw, &
    ' cpcccv final -> ', csrcp, csrcc, csrcv
    go to 600
!
! Calculation of chalcopyrite-bornite-chalcocite association, executed for
! samples with ind_bn = 1
320 continue
    csrcp_raw =((ccrcn*bnrfs)-(ccrfs*bnrcn)+(bnrcn-ccrcn)*fscu_rat   &
                +(ccrfs-bnrfs)*cncu_rat)/d_cpccbn
    csrcc_raw =((bnrcn*cprfs)-(bnrfs*cprcn)+(bnrfs-cprfs)*cncu_rat   &
                +(cprcn-bnrcn)*fscu_rat)/d_cpccbn
    csrbn_raw = 1-csrcp_raw-csrcc_raw
    csrcv_raw = 0.0000

! Check for consistency with bn-cc assumption, calculate bn-cv if negative cc
    if (csrcc_raw .le. -0.005) go to 330

! Set error flag for calculated csr value that exceed tolerances of +/-10%
    if (max(csrcp_raw,csrcc_raw,csrbn_raw).le. 1.1 .and. &
      min(csrcp_raw,csrcc_raw,csrbn_raw).ge. -0.1) then
        err_flag = 0
        go to 321
    end if

    if (tcu < 0.15 .and. min(csrcp_raw,csrcc_raw,csrbn_raw)< -0.100) then
        err_flag = 1         ! Low Grade flag, likely round-off errors or high As value
        if (csrbn_raw < -0.100) then
            if (idebug .eq. 1) &
          write(62,'(a,3f7.3)') ' cpbncc raw ->',csrcp_raw,csrbn_raw,csrcc_raw
            go to 400
        end if
    else if (csrbn_raw < -0.100) then
        err_flag = 2         ! Elevated FSCu relative to CNCu, potentially oxidized
        if (idebug .eq. 1) &
          write(62,'(a,3f7.3)') ' cpbncc raw ->',csrcp_raw,csrbn_raw,csrcc_raw
        go to 400
    else if (csrcp_raw < -0.100) then
        err_flag = 3         ! Elevated CNCu relative to FSCu
    else
        err_flag = 4         ! Undefined error
    end if
!
! Adjust mineral csr for csrcp < 0 and csrbn < 0
321 continue
    if (csrcp_raw .lt. 0.000) then
        csrcp1 = 0.0000
    else if (csrbn_raw .lt. 0.000) then
        csrcp1 = csrcp_raw/(csrcp_raw + csrcc_raw)
    else
        csrcp1 = csrcp_raw
    end if

    if (csrbn_raw .lt. 0.000) then
        csrcc1 = csrcc_raw/(csrcp_raw + csrcc_raw)
        csrbn1 = 0.0000
    else
        csrcc1 = csrcc_raw
        csrbn1 = csrbn_raw
    end if
!
! Second-stage adjustment, for csrcp1 > 1 and csrcc1 < 0
    if (csrcp1 .gt. 1.000) then
        csrcp = csrcp1/(csrcp1 + csrbn1)
    else
        csrcp = csrcp1
    end if

    if (csrcc1 .lt. 0.000) then
        csrcc = 0.0000
    else
        csrcc = csrcc1/(csrcp1 + csrcc1 + csrbn1)
    end if

    csrbn = 1 - csrcp - csrcc
    csrcv = 0.0000
    csrgo = 0.0000
    csrbr = 0.0000
    csrxc = 0.0000
    
    if(idebug .eq. 1) &
    write(62,'(a,3f7.3,a,3f7.3)') ' cpbncc raw ->',csrcp_raw,csrbn_raw,csrcc_raw, &
    ' cpbncc final -> ', csrcp, csrbn, csrcc
    go to 600
!
! Calculation of chalcopyrite-bornite-covellite association, executed for
! samples with ind_bn = 1 and negative chalcocite from previous section
330 continue
    csrcp_raw =((cvrcn*bnrfs)-(cvrfs*bnrcn)+(bnrcn-cvrcn)*fscu_rat   &
                +(cvrfs-bnrfs)*cncu_rat)/d_cpcvbn
    csrcv_raw =((bnrcn*cprfs)-(bnrfs*cprcn)+(bnrfs-cprfs)*cncu_rat   &
                +(cprcn-bnrcn)*fscu_rat)/d_cpcvbn
    csrbn_raw = 1-csrcp_raw-csrcv_raw
    csrcc_raw = 0.0

! Set error flag for calculated csr value that exceed tolerances of 10%
    if (max(csrcp_raw,csrcv_raw,csrbn_raw).le. 1.1 .and. &
      min(csrcp_raw,csrcv_raw,csrbn_raw).ge. -0.1) then
        err_flag = 0
        go to 331
    end if

    if (tcu < 0.15 .and. min(csrcp_raw,csrcv_raw,csrbn_raw)< -0.1) then
        err_flag = 1         ! Low Grade flag, likely round-off errors or high As value
    else if (csrbn_raw > 1.1) then
        err_flag = 2         ! Elevated FSCu relative to CNCu, in cp-cc-bn association
    else if (csrcp_raw < -0.1) then
        err_flag = 3         ! Elevated CNCu relative to FSCu
    else
        err_flag = 4         ! Undefined error
    end if
!
! Adjust mineral csr for csrcp < 0 and csrbn < 0
331 continue
    if (csrcp_raw .lt. 0.0) then
        csrcp1 = 0.0
    else if (csrcv_raw .lt. 0.0) then
        csrcp1 = csrcp_raw/(csrcp_raw + csrbn_raw)
    else
        csrcp1 = csrcp_raw
    end if

    if (csrcv_raw .lt. 0.0) then
        csrbn1 = csrbn_raw/(csrcp_raw + csrbn_raw)
        csrcv1 = 0.0
    else
        csrcv1 = csrcv_raw
        csrbn1 = csrbn_raw
    end if
!
! Second-stage adjustment, for csrcp1 > 1 and csrbn1 < 0
    if (csrcp1 .gt. 1.0) then
        csrcp = csrcp1/(csrcp1 + csrcv1)
     else
        csrcp = csrcp1
    end if

    if (csrbn1 .lt. 0.0) then
        csrbn = 0.0
    else
        csrbn = csrbn1/(csrcp1 + csrcv1 + csrbn1)
    end if

    csrcv = 1.0000 - csrcp - csrbn
    csrcc = 0.0000
    csrgo = 0.0000
    csrbr = 0.0000
    csrxc = 0.0000
    
    if(idebug .eq. 1) &
    write(62,'(a,3f7.3,a,3f7.3)') ' cpbncv raw ->',csrcp_raw,csrbn_raw,csrcv_raw, &
    ' cpbncv final -> ', csrcp, csrbn, csrcv
    go to 600
! -----------------------------------------------------------------------------
! End of Sulfide mineral associations

! Start Oxide/Mixed mineral associations
! -----------------------------------------------------------------------------
400 continue
! Initialize modified assay and calculate PtXt ratios
    scu_rat = scu/tcu
    fscu_rat = fscu/tcu
    cncu_rat = cncu/tcu
    tcumod = tcu
! Initialize bornite and enargite components
    csrbn = 0.0000
    en_wtpct = 0.000
    cspen = 0.000

! Check for sulfide sulfur, Oxide to goe-cc-br-xc
    if (s2 .le. 0.15 .and. s2 .ge. 0) go to 410

! Check for cpy as insol mineral, go to cp-cc-br-xc
    if (cncu > scu) go to 420

! Calculation of normative minerals, goe-cc-br-xc mineral association
410 continue
! Initialize normative mineral matrices
! Mineral 1 = Cugoethite, Mineral 2 = Chalcocite, Mineral 3 = Brochantite
! Mineral 4 = Chrysocolla
    m_min1 = m_goccbrxc
    m_min2 = m_goccbrxc
    m_min3 = m_goccbrxc
    m_min4 = m_goccbrxc
! Insert PtXt ratios into the end-member mineral matrices and find determinant
    m_min1(2,1) = cncu_rat
    m_min1(3,1) = fscu_rat
    m_min1(4,1) = scu_rat
    d_min1 = m44det (m_min1)

    m_min2(2,2) = cncu_rat
    m_min2(3,2) = fscu_rat
    m_min2(4,2) = scu_rat
    d_min2 = m44det (m_min2)

    m_min3(2,3) = cncu_rat
    m_min3(3,3) = fscu_rat
    m_min3(4,3) = scu_rat
    d_min3 = m44det (m_min3)

    m_min4(2,4) = cncu_rat
    m_min4(3,4) = fscu_rat
    m_min4(4,4) = scu_rat
    d_min4 = m44det (m_min4)
    
! Calculation of raw CSR values
    csrgo_raw = d_min1/d_goccbrxc
    csrcc_raw = d_min2/d_goccbrxc
    csrbr_raw = d_min3/d_goccbrxc
    csrxc_raw = d_min4/d_goccbrxc
    
! Initialize error flag and check for valid csr calculations 
    err_flag = 0
    csr_max = max(csrgo_raw,csrcc_raw,csrbr_raw,csrxc_raw)
    csr_min = min(csrgo_raw,csrcc_raw,csrbr_raw,csrxc_raw)
    if (csr_max .lt. 1.00001 .and. csr_min .gt. -0.00001) then
        csrgo1 = csrgo_raw
        csrcc1 = csrcc_raw
        csrbr1 = csrbr_raw
        csrxc1 = csrxc_raw
        go to 411
    end if
    
! Check for error in each end-member correct CSR and set error flag if necessary
    if (csrbr_raw .lt. 0.000) then
        if (csr_max .gt. 1.100 .or. csr_min .lt. -0.100) err_flag = 5
        csrbr1 = 0.0000
        csrcc1 = max (0.0000, csrbr_raw + csrcc_raw)
        csrxc1 = max (0.0000, csrxc_raw)
        csrgo1 = max (0.0000, csrgo_raw)
        go to 411
    end if
    
    if (csrcc_raw .lt. 0.000 .or. csrxc_raw .lt. 0.000) then
        if (csr_max .gt. 1.100 .or. csr_min .lt. -0.100) err_flag = 6
        csrgo1 = max (0.0000, csrgo_raw)
        csrcc1 = max (0.0000, csrcc_raw)
        csrxc1 = max (0.0000, csrxc_raw)
        csrbr1 = max (0.0000, csrbr_raw)
        go to 411    
    end if

    if (csrgo_raw .lt. 0.000) then
        if (csr_max .gt. 1.100 .or. csr_min .lt. -0.100) err_flag = 7
        csrgo1 = 0.0000
        csrcc1 = max (0.000, csrcc_raw)
        csrxc1 = max (0.000, csrxc_raw)
        csrbr1 = max (0.000, csrbr_raw)
        go to 411    
    end if
    
    if (csrmax .gt. 1.100) err_flag = 8
    csrgo1 = csrgo_raw
    csrcc1 = csrcc_raw
    csrxc1 = csrxc_raw
    csrbr1 = csrbr_raw
   
411 continue
! Normalize adjusted csr value and set missing normative minerals to zero
    sumcsr1 = csrgo1 + csrcc1 + csrxc1 + csrbr1
    csrgo = csrgo1/sumcsr1
    csrcc = csrcc1/sumcsr1
    csrbr = csrbr1/sumcsr1
    csrxc = csrxc1/sumcsr1
    csrcp = 0.0000
    csrcv = 0.0000
  
    if(idebug .eq. 1) &
    write(62,'(a,4f7.3,a,4f7.3)') ' goccbrxc raw ->',csrgo_raw,csrcc_raw,csrbr_raw, &
    csrxc_raw,' goccbrxc final -> ', csrgo, csrcc, csrbr,csrxc
    go to 600
 
! Calculation of normative minerals, cp-cc-br-xc mineral association
420 continue
! Initialize normative mineral matrices
! Mineral 1 = Chalcopyrite, Mineral 2 = Chalcocite, Mineral 3 = Brochantite
! Mineral 4 = Chrysocolla
    m_min1 = m_cpccbrxc
    m_min2 = m_cpccbrxc
    m_min3 = m_cpccbrxc
    m_min4 = m_cpccbrxc

! Insert PtXt ratios into the end-member mineral matrices and find determinant
    m_min1(2,1) = cncu_rat
    m_min1(3,1) = fscu_rat
    m_min1(4,1) = scu_rat
    d_min1 = m44det (m_min1)

    m_min2(2,2) = cncu_rat
    m_min2(3,2) = fscu_rat
    m_min2(4,2) = scu_rat
    d_min2 = m44det (m_min2)

    m_min3(2,3) = cncu_rat
    m_min3(3,3) = fscu_rat
    m_min3(4,3) = scu_rat
    d_min3 = m44det (m_min3)

    m_min4(2,4) = cncu_rat
    m_min4(3,4) = fscu_rat
    m_min4(4,4) = scu_rat
    d_min4 = m44det (m_min4)
    
! Calculation of raw CSR values
    csrcp_raw = d_min1/d_cpccbrxc
    csrcc_raw = d_min2/d_cpccbrxc
    csrbr_raw = d_min3/d_cpccbrxc
    csrxc_raw = d_min4/d_cpccbrxc
    
! Initialize error flag and check for valid csr calculations 
    err_flag = 0
    csr_max = max(csrcp_raw,csrcc_raw,csrbr_raw,csrxc_raw)
    csr_min = min(csrcp_raw,csrcc_raw,csrbr_raw,csrxc_raw)
    if (csr_max .lt. 1.000 .and. csr_min .gt. 0.000) then
        csrcp1 = csrcp_raw
        csrcc1 = csrcc_raw
        csrbr1 = csrbr_raw
        csrxc1 = csrxc_raw
        go to 421
    end if
    
! Check for error in each end-member correct CSR and set error flag if necessary
! Samples with negative chrysocolla is diagnostic of normative CV, take cpcccvbr branch
    if (csrxc_raw .lt. -0.010) then
        if (idebug .eq. 1) &
        write(62,'(a,4f7.3)') ' cpccbrxc raw ->',csrcp_raw,csrcc_raw,csrbr_raw,csrxc_raw 
        go to 430
    end if
     
    if (csrcp_raw .lt. 0.000) then
        if (csrbr_raw .lt. 0.000) then
             if (csr_max .gt. 1.100 .or. csr_min .lt. -0.100) err_flag = 9
             csrxc1 = max (0.0000, csrxc_raw + csrcp_raw)
             csrcc1 = max (0.0000, csrcc_raw + csrbr_raw)
             csrcp1 = 0.0000
             csrbr1 = 0.0000
             go to 421
        end if
        csrcp1 = 0.0000
        csrcc1 = max (0.0000, csrcc_raw)
        csrbr1 = max (0.0000, csrbr_raw)
        csrxc1 = max (0.0000, csrxc_raw)
        go to 421
    end if
    
    if (csr_max .gt. 1.100 .or. csr_min .lt. -0.100) err_flag = 10
    csrcp1 = max (0.000, csrcp_raw)
    csrcc1 = max (0.000, csrcc_raw)
    csrxc1 = max (0.000, csrxc_raw)
    csrbr1 = max (0.000, csrbr_raw)
   
421 continue
! Normalize adjusted csr value and set missing normative minerals to zer
    sumcsr1 = csrcp1 + csrcc1 + csrxc1 + csrbr1
    csrcp = csrcp1/sumcsr1
    csrcc = csrcc1/sumcsr1
    csrbr = csrbr1/sumcsr1
    csrxc = csrxc1/sumcsr1
    csrgo = 0.0000
    csrcv = 0.0000
    
    if(idebug .eq. 1) &
    write(62,'(a,4f7.3,a,4f7.3)') ' cpccbrxc raw ->',csrcp_raw,csrcc_raw,csrbr_raw, &
    csrxc_raw,' cpccbrxc final -> ', csrcp, csrcc, csrbr,csrxc
    go to 600
 
! Calculation of normative minerals, cp-cc-cv-br mineral association
430 continue
! Initialize normative mineral matrices
! Mineral 1 = Chalcopyrite, Mineral 2 = Chalcocite, Mineral 3 = Covellite
! Mineral 4 = Brochantite
    m_min1 = m_cpcccvbr
    m_min2 = m_cpcccvbr
    m_min3 = m_cpcccvbr
    m_min4 = m_cpcccvbr
    
! Insert PtXt ratios into the end-member mineral matrices and find determinant
    m_min1(2,1) = cncu_rat
    m_min1(3,1) = fscu_rat
    m_min1(4,1) = scu_rat
    d_min1 = m44det (m_min1)

    m_min2(2,2) = cncu_rat
    m_min2(3,2) = fscu_rat
    m_min2(4,2) = scu_rat
    d_min2 = m44det (m_min2)

    m_min3(2,3) = cncu_rat
    m_min3(3,3) = fscu_rat
    m_min3(4,3) = scu_rat
    d_min3 = m44det (m_min3)

    m_min4(2,4) = cncu_rat
    m_min4(3,4) = fscu_rat
    m_min4(4,4) = scu_rat
    d_min4 = m44det (m_min4)
    
! Calculation of raw CSR values
    csrcp_raw = d_min1/d_cpcccvbr
    csrcc_raw = d_min2/d_cpcccvbr
    csrcv_raw = d_min3/d_cpcccvbr
    csrbr_raw = d_min4/d_cpcccvbr
     
! Initialize error flag and check for valid csr calculations 
    err_flag = 0
    csr_max = max(csrcp_raw,csrcc_raw,csrbr_raw,csrcv_raw)
    csr_min = min(csrcp_raw,csrcc_raw,csrbr_raw,csrcv_raw)
    if (csr_max .lt. 1.000 .and. csr_min .gt. 0.000) then
        csrcp1 = csrcp_raw
        csrcc1 = csrcc_raw
        csrcv1 = csrcv_raw
        csrbr1 = csrbr_raw
        go to 431
    end if
    
! Check for error in each end-member correct CSR and set error flag if necessary

    if (csr_max .gt. 1.100 .or. csr_min .lt. -0.100) err_flag = 11
    csrcp1 = max (0.000, csrcp_raw)
    csrcc1 = max (0.000, csrcc_raw)
    csrcv1 = max (0.000, csrcv_raw)
    csrbr1 = max (0.000, csrbr_raw)
   
431 continue
! Normalize adjusted csr value and set missing normative minerals to zer
    sumcsr1 = csrcp1 + csrcc1 + csrcv1 + csrbr1
    csrcp = csrcp1/sumcsr1
    csrcc = csrcc1/sumcsr1
    csrbr = csrbr1/sumcsr1
    csrcv = csrcv1/sumcsr1
432 continue
    csrgo = 0.0000
    csrxc = 0.0000
    
    if(idebug .eq. 1) &
    write(62,'(a,4f7.3,a,4f7.3)') ' cpcccvbr raw ->',csrcp_raw,csrcc_raw,csrcv_raw, &
    csrbr_raw,' cpcccvbr final -> ', csrcp, csrcc, csrcv,csrbr
    go to 600
!  End of Oxide normative CSR calculations
!------------------------------------------------------------------------------
!
! Calculation of CSP, copper mineral abundance, and pyrite
600 continue
! Set Low Grade error flag
    if (err_flag > 0 .and. tcu .lt. 0.15) err_flag = 1

! Calculation of Copper Source Percent (CSP))
    cspcp = tcumod * csrcp
    cspcc = tcumod * csrcc
    cspcv = tcumod * csrcv
    cspbn = tcumod * csrbn
    cspgo = tcumod * csrgo
    cspxc = tcumod * csrxc
    cspbr = tcumod * csrbr
! Calculation of mineral abundance (wt. %)
    cp_wtpct = cspcp / cpxcu
    cc_wtpct = cspcc / ccxcu
    cv_wtpct = cspcv / cvxcu
    bn_wtpct = cspbn / bnxcu
    go_wtpct = cspgo / goxcu
    xc_wtpct = cspxc / xcxcu
    br_wtpct = cspbr / brxcu

! Calculation of pyrite based on Fe content
    if (fe .lt. 0) then
        py_fe = -99.
        xfe = -99.
    else
        cusul_fe = cp_wtpct*cpxfe + cc_wtpct*ccxfe + cv_wtpct*cvxfe &
                 + bn_wtpct*bnxfe + en_wtpct*enxfe
        py_fe = max (0.000, (fe - cusul_fe)/0.466)
    end if

! Calculation of pyrite based on S content
   if (s2 .lt. 0) then
        py_best = -99.             ! Sample missing both S and Fe assay
        xfe = -99.
        go to 610
   else
        cusul_s = cp_wtpct*cpxs + cc_wtpct*ccxs + cv_wtpct*cvxs  &
                 + bn_wtpct*bnxs + en_wtpct*enxs
        py_s = max (0.0, (s2 - cusul_s)/0.534)
   end if

! Select minimum pyrite value as best estimate.  Pooled analytical uncertainty is 11%
    
   if (py_fe .lt. 0.0)then
        py_best = py_s                   ! Check for missing Fe assay only
        xfe = -99.
        go to 610
   end if
    
    ave_est = (py_fe+py_s)/2.0
    if (abs(ave_est) .lt. 0.001) then
        py_best = 0.000                    ! Check for divide by zero
        go to 610
    end if

    if ((abs(py_fe-py_s)/ave_est).le. 0.11) then
        py_best = ave_est
    else
        py_best = min(py_fe,py_s)
    end if

! Calculate excess iron (non-sulfide Fe)
    xfe = fe - cusul_fe - py_best*0.466
    if (xfe .lt. 0.0) xfe = 0.0
    
610 continue

! Calculate volume percent sulfides from user-supplied whole rock density and
! mineral densities published by mindat.org
!
! Assignment of mineral densities (g/cc)from www.mindat.com 
    sg_py = 5.01     ! Calculated value pyrite
    sg_cp = 4.18     ! Calculated value chalcopyrite
    sg_cc = 5.75     ! Calculated for djurleite, near top of measured cc range
    sg_cv = 4.60     ! Calculated value covellite
    sg_bn = 5.09     ! Calculated value bornite
    sg_en = 4.40     ! Calculated value enargite
    sg_xc = 2.10     ! Mid-range of measured values, chrysocolla
    sg_br = 3.97     ! Measured value brochantite
    sg_go = 4.28     ! Mid-range of measured values, goethite
    
!
! Rock density provided by data base, with default value assigned where measured
! value is missing
!
! If default density is set to zero, skip calculations
    if (density .lt. 0.01)then
        cp_volpct = -99.
        cc_volpct = -99.
        cv_volpct = -99.
        bn_volpct = -99.
        en_volpct = -99.
        py_volpct = -99.
        xc_volpct = -99.
        br_volpct = -99.
        go_volpct = -99.
        go to 620
    end if      
!
! Calculation of volume percent
    cp_volpct = cp_wtpct * density / sg_cp
    cc_volpct = cc_wtpct * density / sg_cc
    cv_volpct = cv_wtpct * density / sg_cv
    bn_volpct = bn_wtpct * density / sg_bn
    en_volpct = en_wtpct * density / sg_en
    xc_volpct = xc_wtpct * density / sg_xc
    br_volpct = br_wtpct * density / sg_br
    go_volpct = go_wtpct * density / sg_go
    if (py_best .lt. 0) then
        py_volpct = -99.
    else
        py_volpct = py_best * density / sg_py
    end if
!
!================================================================================
! Calculations complete 
!  Print normative mineralogy to file using cvs formatting functions
! THIS SECTION TO BE REPLACED BY DIRECT IO TO DATABASE
!================================================================================
! Modified Aug 2022, Version 2.3, Added tcu, fe, and s2 to output record    
620 continue
    rec_out = ' '

    call csv_record_append_s  (dhname, rec_out)
    call csv_record_append_r4 (f_int, rec_out)
    call csv_record_append_r4 (t_int , rec_out)
    call csv_record_append_s  (sampid, rec_out)
    call csv_record_append_i4 (err_flag, rec_out)
    call csv_record_append_i4 (minzon, rec_out)
    call csv_record_append_i4 (ind_bn, rec_out)
    call csv_record_append_r4 (tcu, rec_out)
    call csv_record_append_r4 (fe, rec_out)
    call csv_record_append_r4 (s2, rec_out)
    call csv_record_append_r4 (csrcp, rec_out)
    call csv_record_append_r4 (csrcc, rec_out)
    call csv_record_append_r4 (csrcv, rec_out)
    call csv_record_append_r4 (csrbn, rec_out)
    call csv_record_append_r4 (csrgo, rec_out)
    call csv_record_append_r4 (csrbr, rec_out)
    call csv_record_append_r4 (csrxc, rec_out)
    call csv_record_append_r4 (cspcp, rec_out)
    call csv_record_append_r4 (cspcc, rec_out)
    call csv_record_append_r4 (cspcv, rec_out)
    call csv_record_append_r4 (cspbn, rec_out)
    call csv_record_append_r4 (cspgo, rec_out)
    call csv_record_append_r4 (cspbr, rec_out)
    call csv_record_append_r4 (cspxc, rec_out)
    call csv_record_append_r4 (cspen, rec_out)
    call csv_record_append_r4 (cp_wtpct, rec_out)
    call csv_record_append_r4 (cc_wtpct, rec_out)
    call csv_record_append_r4 (cv_wtpct, rec_out)
    call csv_record_append_r4 (bn_wtpct, rec_out)
    call csv_record_append_r4 (go_wtpct, rec_out)
    call csv_record_append_r4 (br_wtpct, rec_out)
    call csv_record_append_r4 (xc_wtpct, rec_out)
    call csv_record_append_r4 (en_wtpct, rec_out)
    call csv_record_append_r4 (py_best, rec_out)
    call csv_record_append_r4 (xfe, rec_out)
    call csv_record_append_r4 (cp_volpct, rec_out)
    call csv_record_append_r4 (cc_volpct, rec_out)
    call csv_record_append_r4 (cv_volpct, rec_out)
    call csv_record_append_r4 (bn_volpct, rec_out)
    call csv_record_append_r4 (go_volpct, rec_out)
    call csv_record_append_r4 (br_volpct, rec_out)
    call csv_record_append_r4 (xc_volpct, rec_out)
    call csv_record_append_r4 (en_volpct, rec_out)
    call csv_record_append_r4 (py_volpct, rec_out)
    
    call csv_file_record_write (outputcsv, 61, rec_out )
!================================================================================
! Record output complete 
!================================================================================
    
! Write progress indicator to screen
    if (mod(ndata,50).eq. 0) write (*,'(i6,a)') ndata,' records processed'

! Read next line
      go to 200

! Data file read. End program
900 continue
    ndata = ndata - 1
    write (*,'(i6,a)') ndata, " records processed. Run completed."
    close (51)
    close (61)
    if (idebug .eq. 1) close (62)
    stop
END PROGRAM norminv2
!
!=================================================================================
!=================================================================================
! Begin subroutine and function programming. 
! These are subprograms that are called in the MAIN program
!---------------------------------------------------------------------------------
! Begin subroutines for input and output of CSV formatted files
!---------------------------------------------------------------------------------
! Subroutines to parse data record from input record
! Separate subroutines for character, integer, and R4 real data types
! Minor error checking - when error is encountered, variable is set to 'missing'
! indicator, record_status is set to value > 0, and control is returned to calling
! statement.  
!---------------------------------------------------------------------------------

subroutine get_char (xdata,icol,csv_record,csv_value,csv_len,record_status)
! Subroutine to extract caracter data type from csv record
! icol is location in reoord, from 1 to value count as per parameter file input
! Written by Richard K. Preece
! August 19, 2019
   
    
    character (len = 25) xdata, cdata
    character csv_char
    character csv_char_old
    integer ( kind = 4) csv_len
    integer ( kind = 4) csv_loc
    character ( len = *) csv_record
    integer ( kind = 4) record_status
    integer ( kind = 4) value_count
    integer ( kind = 4) csv_value
    integer ( kind = 4) word_length

! Initiate counters
    value_count = 0
    record_status = 0
    word_length = 0
    csv_char_old = ','
    
! Check to determine if variable is missing from the input file     
    if (icol .eq. 0) then
        xdata = '-99'
        return
    end if

    do csv_loc = 1, csv_len
        csv_char = csv_record(csv_loc:csv_loc)
        !
        !  Check for start of new value field
        if (csv_char_old == ',') then
            value_count = value_count + 1
            word_length = 1
        else
             word_length = word_length + 1
        end if
            if (csv_char == ',') then
                if (value_count == icol) go to 100
            end if
        csv_char_old = csv_char
    end do

    if (csv_loc == csv_len) then
        ! data column not found, set data to -99 and return with error code
        record_status = 1
        xdata = ''
        return
    else
        !Unknown reason for early termination, set value to -99 and return with error code
        record_status = 9
        xdata = ''
        go to 101
    end if

    ! Data column found. Calculate column location, read real value, and return
100 continue
    mloc = csv_loc
    nloc = mloc - (word_length-1)
    read (csv_record(nloc:mloc), *,err=101,iostat=record_status) cdata
    xdata = trim(cdata) 
    
101 continue
    
return
end

subroutine get_integer (xdata,icol,csv_record,csv_value,csv_len,record_status)
! Subroutine to extract integer data type from csv record
! icol is location in record, from 1 to value count as per parameter file input
! Written by Richard K. Preece
! August 19, 2019
   
    
    integer (kind = 4) xdata
    character csv_char
    character csv_char_old
    integer ( kind = 4) csv_len
    integer ( kind = 4) csv_loc
    character ( len = *) csv_record
    integer ( kind = 4) record_status
    integer ( kind = 4) value_count
    integer ( kind = 4) csv_value
    integer ( kind = 4) word_length

! Initiate counters
    value_count = 0
    record_status = 0
    word_length = 0
    csv_char_old = ','
    
! Check to determine if variable is missing from the input file     
    if (icol .eq. 0) then
        xdata = -99
        return
    end if

    do csv_loc = 1, csv_len
        csv_char = csv_record(csv_loc:csv_loc)
        !
        !  Check for start of new value field
        if (csv_char_old == ',') then
            value_count = value_count + 1
            word_length = 1
        else
             word_length = word_length + 1
        end if
            if (csv_char == ',') then
                if (value_count == icol) go to 100
            end if
        csv_char_old = csv_char
    end do

    if (csv_loc == csv_len) then
        ! data column not found, set data to -99 and return with error code
        record_status = 1
        xdata = -99
        return
    else
        !Unknown reason for early termination, set value to -99 and return with error code
        record_status = 9
        xdata = -99
        go to 101
    end if

    ! Data column found. Calculate column location, read real value, and return
100 continue
    mloc = csv_loc
    nloc = mloc - (word_length-1)
    read (csv_record(nloc:mloc), *,err=101,iostat=record_status) zdata
    xdata = int(zdata)

101 continue
    
return
end

subroutine get_real (xdata,icol,csv_record,csv_value,csv_len,record_status)
! Subroutine to extract real data type from csv record
! icol is location in reoord, from 1 to value count as per parameter file input
! Written by Richard K. Preece
! August 16, 2019
   
    
    real (kind = 4) xdata
    character csv_char
    character csv_char_old
    integer ( kind = 4) csv_len
    integer ( kind = 4) csv_loc
    character ( len = *) csv_record
    integer ( kind = 4) record_status
    integer ( kind = 4) value_count
    integer ( kind = 4) csv_value
    integer ( kind = 4) word_length

! Initiate counters
    value_count = 0
    record_status = 0
    word_length = 0
    csv_char_old = ','
    
! Check to determine if variable is missing from the input file     
    if (icol .eq. 0) then
        xdata = -99.0
        return
    end if

    do csv_loc = 1, csv_len
        csv_char = csv_record(csv_loc:csv_loc)
        !
        !  Check for start of new value field
        if (csv_char_old == ',') then
            value_count = value_count + 1
            word_length = 1
        else
             word_length = word_length + 1
        end if
            if (csv_char == ',') then
                if (value_count == icol) go to 100                 
            end if
        csv_char_old = csv_char
    end do
    
    ! Check if item is in last value on line
    if (value_count == icol) then
        word_length = word_length + 1
        go to 100
    end if

    if (csv_loc == csv_len) then
        ! data column not found, set data to -99 and return with error code
        record_status = 1
        xdata = -99.0
        return
    else
        !Unknown reason for early termination, set value to -99 and return with error code
        record_status = 9
        xdata = -99.0
        go to 101
    end if

    ! Data column found. Calculate column location, read real value, and return
100 continue
    mloc = csv_loc
    nloc = mloc - (word_length-1)
    read (csv_record(nloc:mloc), *,err=101,iostat=record_status) xdata

101 continue
    
return
end
!----------------------------------------------------------------------------------
! Public domain subroutines to manage CSV-formatted input and output files
!----------------------------------------------------------------------------------
subroutine csv_value_count(csv_record, csv_record_status, value_count, csv_len)

    !*****************************************************************************80
    !
    !! CSV_COUNT counts the number of values in a CSV record.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    implicit none

    character csv_char
    character csv_char_old
    integer ( kind = 4) csv_len
    integer ( kind = 4) csv_loc
    character ( len = *) csv_record
    integer ( kind = 4) csv_record_status
    character :: TAB = achar(9)
    integer ( kind = 4) value_count
    integer ( kind = 4) word_length
    !
    !  No values so far.
    !
    value_count = 0
    !
    !  We begin in "unquoted" status.
    !
    csv_record_status = 0
    !
    !  How many characters in the record?
    !
    csv_len = len_trim(csv_record)
    !
    !  Count number of characters in each word.
    !
    word_length = 0
    !
    !  Consider each character.
    !
    csv_char_old = ','

    do csv_loc = 1, csv_len

        csv_char = csv_record(csv_loc:csv_loc)
        !
        !  Each comma divides one value from another.
        !
        if (csv_char_old == ',') then

            value_count = value_count + 1
            word_length = 0
            !
            !  For quotes, try using CSV_RECORD_STATUS to count the number of
            !  quoted characters.
            !
        else if (csv_char == '"') then

            if (0 < csv_record_status) then
            csv_record_status = 0
        else
            csv_record_status = csv_record_status + 1
            end if
            !
            !  Ignore blanks
            !
        else if (csv_char == ' ' .or. csv_char == TAB) then
            !
            !  Add character to length of word.
            !
        else

            word_length = word_length + 1

            if (value_count == 0) then
            value_count = 1
            end if

        end if

        csv_char_old = csv_char

    end do

    return
end
subroutine csv_file_close_read(csv_file_name, csv_file_unit)

    !*****************************************************************************80
    !
    !! CSV_FILE_CLOSE_READ closes a CSV file for reading.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
    !
    !    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
    !
    implicit none

    character ( len = *) csv_file_name
    integer ( kind = 4) csv_file_unit

    close ( unit = csv_file_unit)

    return
end
subroutine csv_file_close_write(csv_file_name, csv_file_unit)

    !*****************************************************************************80
    !
    !! CSV_FILE_CLOSE_WRITE closes a CSV file for writing.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
    !
    !    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
    !
    implicit none

    character ( len = *) csv_file_name
    integer ( kind = 4) csv_file_unit

    close ( unit = csv_file_unit)

    return
end
subroutine csv_file_header_write(csv_file_name, csv_file_unit, header)

    !*****************************************************************************80
    !
    !! CSV_FILE_HEADER_WRITE writes a header to a CSV file.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
    !
    !    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
    !
    !    Input, character ( len = * ) HEADER, the header.
    !
    implicit none

    character ( len = *) csv_file_name
    integer ( kind = 4) csv_file_unit
    character ( len = *) header

    write ( csv_file_unit, '(a)') trim(header)

    return
end
subroutine csv_file_line_count(csv_file_name, line_num)

    !*****************************************************************************80
    !
    !! CSV_FILE_LINE_COUNT counts the number of lines in a CSV file.
    !
    !  Discussion:
    !
    !    This routine does not try to distinguish the possible header line,
    !    blank lines, or cases where a single CSV record extends over multiple
    !    lines.  It simply counts the number of lines.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
    !
    !    Output, integer ( kind = 4 ) LINE_NUM, the number of lines.
    !
    implicit none

    character ( len = *) csv_file_name
    integer ( kind = 4) ierror
    integer ( kind = 4) input_status
    integer ( kind = 4) input_unit
    character ( len = 1023) line
    integer ( kind = 4) line_num

    line_num = -1

    call get_unit(input_unit)

    open ( unit = input_unit, file = csv_file_name, status = 'old', &
    iostat = input_status)

    if (input_status /= 0) then
        write ( *, '(a)') ' '
        write ( *, '(a)') 'CSV_FILE_LINE_COUNT - Fatal error!'
        write ( *, '(a,i8)') '  Could not open "' // trim(csv_file_name) // '".'
        stop
    end if

    line_num = 0

    do

        read ( input_unit, '(a)', iostat = input_status) line

        if (input_status /= 0) then
            ierror = line_num
            exit
        end if

        line_num = line_num + 1

    end do

    close ( unit = input_unit)

    return
end
subroutine csv_file_record_write(csv_file_name, csv_file_unit, record)

    !*****************************************************************************80
    !
    !! CSV_FILE_RECORD_WRITE writes a record to a CSV file.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
    !
    !    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
    !
    !    Input, character ( len = * ) RECORD, the record.
    !
    implicit none

    character ( len = *) csv_file_name
    integer ( kind = 4) csv_file_unit
    character ( len = *) record

    write ( csv_file_unit, '(a)') trim(record)

    return
end
subroutine csv_file_open_read(csv_file_name, csv_file_unit)

    !*****************************************************************************80
    !
    !! CSV_FILE_OPEN_READ opens a CSV file for reading.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
    !
    !    Output, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
    !
    implicit none

    character ( len = *) csv_file_name
    integer ( kind = 4) csv_file_status
    integer ( kind = 4) csv_file_unit

    call get_unit(csv_file_unit)

    open ( unit = csv_file_unit, file = csv_file_name, status = 'old', &
    iostat = csv_file_status)

    if (csv_file_status /= 0) then
        write ( *, '(a)') ' '
        write ( *, '(a)') 'CSV_FILE_OPEN_READ - Fatal error!'
        write ( *, '(a,i8)') '  Could not open "' // trim(csv_file_name) // '".'
        csv_file_unit = -1
        stop
    end if

    return
end
subroutine csv_file_open_write(csv_file_name, csv_file_unit)

    !*****************************************************************************80
    !
    !! CSV_FILE_OPEN_WRITE opens a CSV file for writing.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
    !
    !    Output, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
    !
    implicit none

    character ( len = *) csv_file_name
    integer ( kind = 4) csv_file_status
    integer ( kind = 4) csv_file_unit

    call get_unit(csv_file_unit)

    open ( unit = csv_file_unit, file = csv_file_name, status = 'replace', &
    iostat = csv_file_status)

    if (csv_file_status /= 0) then
        write ( *, '(a)') ' '
        write ( *, '(a)') 'CSV_FILE_OPEN_WRITE - Fatal error!'
        write ( *, '(a,i8)') '  Could not open "' // trim(csv_file_name) // '".'
        stop
    end if

    return
end
subroutine csv_record_append_i4(i4, record)

    !*****************************************************************************80
    !
    !! CSV_RECORD_APPEND_I4 appends an I4 to a CSV record.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    24 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I4, the integer to be appended
    !
    !    Input/output, character ( len = * ) RECORD, the CSV record.
    !
    implicit none

    character ( len = 5) fmat
    integer ( kind = 4) i
    integer ( kind = 4) i4
    integer ( kind = 4) i4_len
    integer ( kind = 4) i4_width
    character ( len = *) record
    !
    !  Locate last used location in RECORD.
    !
    i = len_trim(record)
    !
    !  Append comma.
    !
    if (0 < i) then
        i = i + 1
        record(i:i) = ','
    end if
    !
    !  Determine "width" of I4.
    !
    i4_len = i4_width(i4)
    !
    !  Create format for I4.
    !
    write ( fmat, '(a,i2,a)') '(i', i4_len, ')'
    !
    !  Write I4 to RECORD.
    !
    write ( record(i + 1:i + i4_len), fmat) i4

    return
end
subroutine csv_record_append_r4(r4, record)

    !*****************************************************************************80
    !
    !! CSV_RECORD_APPEND_R4 appends an R4 to a CSV record.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 November 2008
    !
    !  Author:
    !
    !    John Burkardt, modified by R.K. Preece 19 Aug 2019
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) R4, the value to be appended
    !
    !    Input/output, character ( len = * ) RECORD, the CSV record.
    !
    implicit none

    character ( len = 5) fmat
    character ( len = 7) fmat2
    integer ( kind = 4) i
    integer ( kind = 4) i4
    integer ( kind = 4) i4_len
    integer ( kind = 4) i4_width
    real ( kind = 4) r4
    character ( len = *) record
    !
    !  Locate last used location in RECORD.
    !
    i = len_trim(record)
    !
    !  Append comma.
    !
    if (0 < i) then
        i = i + 1
        record(i:i) = ','
    end if

    if (r4 == 0.0) then
        i = i + 1
        record(i:i) = '0'
    else if (r4 == real ( int ( r4), kind = 4)) then
        i4 = int ( r4)
        i4_len = i4_width(i4)
        write ( fmat, '(a,i2,a)') '(i', i4_len, ')'
        write ( record(i + 1:i + i4_len), fmat) i4
    else
! This branch for writing R4 variable modified by RK Preece, 19 August 2019
! Changed from g18.7 format to variable width f*.3 format
        i4 = int (r4*1000)
        i4_len = i4_width(i4) + 2
        if (i4_len .lt. 5) then
            if (r4 .lt. 0) then
                i4_len = 6
            else
                i4_len = 5
            end if
        end if
        write ( fmat2, '(a,i2,a)') '(f', i4_len, '.3)'
        write ( record(i + 1:i + i4_len), fmat2) r4
!        write ( record(i + 1:i + 8), '(f8.3)') r4
    end if

    return
end
subroutine csv_record_append_r8(r8, record)

    !*****************************************************************************80
    !
    !! CSV_RECORD_APPEND_R8 appends an R8 to a CSV record.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) R8, the value to be appended
    !
    !    Input/output, character ( len = * ) RECORD, the CSV record.
    !
    implicit none

    character ( len = 5) fmat
    integer ( kind = 4) i
    integer ( kind = 4) i4
    integer ( kind = 4) i4_len
    integer ( kind = 4) i4_width
    real ( kind = 8) r8
    character ( len = *) record
    !
    !  Locate last used location in RECORD.
    !
    i = len_trim(record)
    !
    !  Append comma.
    !
    if (0 < i) then
        i = i + 1
        record(i:i) = ','
    end if

    if (r8 == 0.0D+00) then
        i = i + 1
        record(i:i) = '0'
    else if (r8 == real ( int ( r8), kind = 8)) then
        i4 = int ( r8)
        i4_len = i4_width(i4)
        write ( fmat, '(a,i2,a)') '(i', i4_len, ')'
        write ( record(i + 1:i + i4_len), fmat) i4
    else
        write ( record(i + 1:i + 14), '(g14.6)') r8
    end if

    return
end
subroutine csv_record_append_s(s, record)

    !*****************************************************************************80
    !
    !! CSV_RECORD_APPEND_S appends a string to a CSV record.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 November 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) S, the string to be appended
    !
    !    Input/output, character ( len = * ) RECORD, the CSV record.
    !
    implicit none

    integer ( kind = 4) i
    character ( len = *) record
    character ( len = *) s
    integer ( kind = 4) s_len
    !
    !  Locate last used location in RECORD.
    !
    i = len_trim(record)
    !
    !  Append a comma.
    !
    if (0 < i) then
        i = i + 1
        record(i:i) = ','
    end if
    !
    !  Prepend a quote.
    !
    !i = i + 1
    !record(I:i) = '"'
    !
    !  Write S to RECORD.
    !
    s_len = len_trim(s)
    record(i + 1:i + s_len) = s(1:s_len)
    i = i + s_len
    !
    !  Postpend a quote
    !
    !i = i + 1
    !record(i:i) = '"'

    return
end
subroutine get_unit(iunit)

    !*****************************************************************************80
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 October 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) IUNIT, the free unit number.
    !
    implicit none

    integer ( kind = 4) i
    integer ( kind = 4) ios
    integer ( kind = 4) iunit
    logical lopen

    iunit = 0

    do i = 1, 99

        if (i /= 5 .and. i /= 6 .and. i /= 9) then

            inquire ( unit = i, opened = lopen, iostat = ios)

            if (ios == 0) then
                if (.not.lopen) then
                    iunit = i
                    return
                end if
            end if

        end if

    end do

    return
end
function i4_log_10(i)

    !*****************************************************************************80
    !
    !! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
    !
    !  Discussion:
    !
    !    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
    !
    !    An I4 is an integer ( kind = 4 ) value.
    !
    !  Example:
    !
    !        I  I4_LOG_10
    !    -----  --------
    !        0    0
    !        1    0
    !        2    0
    !        9    0
    !       10    1
    !       11    1
    !       99    1
    !      100    2
    !      101    2
    !      999    2
    !     1000    3
    !     1001    3
    !     9999    3
    !    10000    4
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 June 2003
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I, the number whose logarithm base 10
    !    is desired.
    !
    !    Output, integer ( kind = 4 ) I4_LOG_10, the integer part of the
    !    logarithm base 10 of the absolute value of X.
    !
    implicit none

    integer ( kind = 4) i
    integer ( kind = 4) i_abs
    integer ( kind = 4) i4_log_10
    integer ( kind = 4) ten_pow

    if (i == 0) then

        i4_log_10 = 0

    else

        i4_log_10 = 0
        ten_pow = 10

        i_abs = abs(i)

        do while (ten_pow <= i_abs)
            i4_log_10 = i4_log_10 + 1
            ten_pow = ten_pow * 10
        end do

    end if

    return
end
function i4_width(i)

    !*****************************************************************************80
    !
    !! I4_WIDTH returns the "width" of an I4.
    !
    !  Discussion:
    !
    !    The width of an integer is the number of characters necessary to print it.
    !
    !    The width of an integer can be useful when setting the appropriate output
    !    format for a vector or array of values.
    !
    !    An I4 is an integer ( kind = 4 ) value.
    !
    !  Example:
    !
    !        I  I4_WIDTH
    !    -----  -------
    !    -1234    5
    !     -123    4
    !      -12    3
    !       -1    2
    !        0    1
    !        1    1
    !       12    2
    !      123    3
    !     1234    4
    !    12345    5
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I, the number whose width is desired.
    !
    !    Output, integer ( kind = 4 ) I4_WIDTH, the number of characters
    !    necessary to represent the integer in base 10, including a negative
    !    sign if necessary.
    !
    implicit none

    integer ( kind = 4) i
    integer ( kind = 4) i4_log_10
    integer ( kind = 4) i4_width

    if (0 < i) then
        i4_width = i4_log_10(i) + 1
    else if (i == 0) then
        i4_width = 1
    else if (i < 0) then
        i4_width = i4_log_10(i) + 2
    end if

    return
end
subroutine timestamp()

    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 May 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
    implicit none

    character ( len = 8) ampm
    integer ( kind = 4) d
    integer ( kind = 4) h
    integer ( kind = 4) m
    integer ( kind = 4) mm
    character ( len = 9), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
    integer ( kind = 4) n
    integer ( kind = 4) s
    integer ( kind = 4) values(8)
    integer ( kind = 4) y

    call date_and_time(values = values)

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if (h < 12) then
        ampm = 'AM'
    else if (h == 12) then
        if (n == 0 .and. s == 0) then
        ampm = 'Noon'
    else
        ampm = 'PM'
        end if
    else
        h = h - 12
        if (h < 12) then
        ampm = 'PM'
    else if (h == 12) then
        if (n == 0 .and. s == 0) then
        ampm = 'Midnight'
    else
        ampm = 'AM'
        end if
        end if
    end if

    write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
    d, trim(month(m)), y, h, ':', n, ':', s, '.', mm, trim(ampm)

    return
end
!
!-------------------------------------------------------------------------------
! Begin functions for finding the determinant of an algebraic matrix
! These are calls from within an equation of the MAIN program
! In both functions, the array is passed to function A and the calculated determinant
! is returned to the calling equation
!-------------------------------------------------------------------------------
FUNCTION M33DET (A) RESULT (DET)
!***********************************************************************************************************************************
!  M33DET  -  Function to compute the determinant of a 3x3 matrix.
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!  Date:         July 22, 2005
!  Language:     Fortran-90
!  Version:      1.00a
!!***********************************************************************************************************************************

      IMPLICIT NONE
      real (kind = 8), DIMENSION(3,3), INTENT(IN)  :: A
      real (kind = 8) :: DET


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      RETURN

END FUNCTION M33DET
FUNCTION M44DET (A) RESULT (DET)
!***********************************************************************************************************************************
!  Function:     M44DET
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!  Date:         February 6, 2009
!  Language:     Fortran-90
!  Version:      1.00a
!  Description:  Computes the determinant a 4x4 matrix
!***********************************************************************************************************************************

      IMPLICIT NONE
      real (kind = 8), DIMENSION(4,4), INTENT(IN)  :: A
      real (kind = 8) :: DET

      DET =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)- &
             A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ &
             A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)- &
             A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ &
             A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

      RETURN

END FUNCTION M44DET


