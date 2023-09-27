MODULE initialization_mod

!--------------------------------------------------------------
! References:
! Launiainen, S., Guan, M., Salmivaara, A., Kieloaho, A.-J., 2019. 
! Modeling boreal forest evapotranspiration and water balance at stand and catchment scales: a spatial approach. 
! Hydrology and Earth System Sciences 23, 3457–3480. https://doi.org/10.5194/hess-23-3457-2019
!--------------------------------------------------------------

! Modules
  use netcdf             ! library for processing netcdf files
  use readvegpara_mod       ! module for reading parameter files in ASCII
  use readclim_mod       ! module for reading reading meteorological forcing data
  use readsoilpara_mod       ! module for reading soil properties (shared with yasso?)
  use yasso              ! module for reading soil properties (shared with yasso?)

  implicit none

  !Public member functions:
  public :: canopy_water_flux   ! require inputdata for p-hydro module
  public :: topmodel   !
  public :: soil_water
  public :: soil_water_retention_curve

contains
  
  SUBROUTINE initialization(phydro, soilwater, yasso)
  ! -----------------------------
  ! runs a timestep, updates saturation deficit and returns fluxes
  !      Args:
  !          R - recharge [m per unit catchment area] during timestep
  !      OUT:
  !          Qb - baseflow [m per unit area]
  !          Qr - returnflow [m per unit area]
  !          qr - distributed returnflow [m]
  !          fsat - saturated area fraction [-]
  !      Note: 
  !          R is the mean drainage [m] from bucketgrid.
  ! ------------------------------

  ! !ARGUMENTS:
    type(soilwater_state), intent(inout)    :: soilwater_state    
    real(r8)      , intent(in)    :: R       ! recharge [m per unit catchment area] during timestep
    real(r8)      , intent(out)   :: Qb      ! baseflow [m per unit area]
    real(r8)      , intent(out)   :: Qr      ! returnflow [m per unit area]
    real(r8)      , intent(out)   :: qr      ! distributed returnflow [m]
    real(r8)      , intent(out)   :: fsat    ! saturated area fraction [-]

  ! ! Local variables 
    real     ::    S0, s, Qb, S
     
    call readlai()
    
    call readveg_prop()

    call readctrl_namelist()

    call read_namelist()

    call initialization_phydro(phydro)

    ! Read soil properties for spafhy
    call readsoil_prop()

    call initialization_spafhy(soilwater)

    ! Read soil parameters for yasso
    call get_params(param_base, alpha_awen, beta12, decomp_pc, param_final)

    if () then
      ! A steady-state start 
      call initialization_yasso_steadystate(soilc)
    else
      ! A cold start
      cstate(:) = 0.0_r8
      nstate = 0.0_r8
    end if

  END SUBROUTINE initialization

  SUBROUTINE initialization_spafhy(pgen, pcpy, pbu, ptop, psoil, gisdata, ncf=True, ncf_file)
  ! ******************** sets up SpaFHy  **********************
  ! Normal SpaFHy run without parameter optimization
  !  1) gets parameters as input arguments
  !  2) reads GIS-data, here predefined format for Seurantaverkko -cathcments 
  !  3) creates CanopyGrid (cpy), BucketGrid (bu) and Topmodel (top) -objects 
  !      within Spathy-object (spa), and temporary outputs
  !  4) creates netCDF -file for outputs if 'ncf' = True 
  !  5) returns following outputs:
  !      spa - spathy object
  !      ncf - netcdf file handle (if ncf=True)
  !      ncf_file - ncf file name (if ncf=True)     
  !  IN:
  !      pgen, pcpy, pbu, psoil - parameter dictionaries
  !      gisdata - dict of 2d np.arrays containing gis-data with keys:
  !          cmask - catchment mask; integers within np.Nan outside
  !          LAI_conif [m2m-2]
  !          LAI_decid [m2m-2]
  !          hc, canopy closure [m]
  !          fc, canopy closure fraction [-]
  !          soil, soil type integer code 1-5
  !          flowacc - flow accumulation [units]
  !          slope - local surface slope [units]         
  !          cellsize - gridcell size
  !          lon0 - x-grid
  !          lat0 - y-grid
  !
  !  OUT:
  !      spa - spathy object
  ! -----------------------------------------

    ! !ARGUMENTS:
    type(soilwater_state), intent(inout)    :: soilwater_state    
    real(r8)      , intent(in)    :: R       ! recharge [m per unit catchment area] during timestep
    real(r8)      , intent(out)   :: Qb      ! baseflow [m per unit area]
    real(r8)      , intent(out)   :: Qr      ! returnflow [m per unit area]
    real(r8)      , intent(out)   :: qr      ! distributed returnflow [m]
    real(r8)      , intent(out)   :: fsat    ! saturated area fraction [-]
                                      ncf_file

  ! ! Local variables 
    real     ::    S0, s, Qb, S
      
  ! initial state of canopy storage [mm] and snow water equivalent [mm]
    w   = 0.0            ! canopy storage mm
    swe = 0.0            ! snow water equivalent m

  ! if running a regional simulation need to read spatial parameters.
  ! if (region) then
  !   call read_regionsoil()    
  !   call read_regionveg()
  ! end if

  bt=0
  et=(end_date - start_date)*dt

  !flowacc = gisdata['flowacc'].copy()
  !slope = gisdata['slope'].copy()        
        
        """
        flatten=True omits cells outside catchment
        """
        if flatten:
            ix = np.where(np.isfinite(cmask))
            cmask = cmask[ix].copy()
            flowacc = flowacc[ix].copy()
            slope = slope[ix].copy()
            
            for key in pcpy['state']:
                pcpy['state'][key] = pcpy['state'][key][ix].copy()
                        
            for key in soildata:
                soildata[key] = soildata[key][ix].copy()
                
            self.ix = ix  # indices to locate back to 2d grid       
   
        """--- initialize CanopyGrid ---"""
        self.cpy = CanopyGrid(pcpy, pcpy['state'], outputs=cpy_outputs)

        """--- initialize BucketGrid ---"""
        self.bu = BucketGrid(spara=soildata, outputs=bu_outputs)

        """ --- initialize Topmodel --- """
        self.top=Topmodel(ptop, self.GisData['cellsize']**2, cmask,
                          flowacc, slope, outputs=top_outputs)


   def __init__(self, cpara, state, outputs=False):
        """
        initializes CanopyGrid -object

        Args:
            cpara - parameter dict:
            state - dict of initial state
            outputs - True saves output grids to list at each timestep

        Returns:
            self - object
            
        NOTE:
            Currently the initialization assumes simulation start 1st Jan, 
            and sets self._LAI_decid and self.X equal to minimum values.
            Also leaf-growth & senescence parameters are intialized to zero.
        """
        epsi = 0.01 # small number
        
        self.Lat = cpara['loc']['lat']
        self.Lon = cpara['loc']['lon']

        # physiology: transpi + floor evap
        self.physpara = cpara['physpara']

        # phenology & LAI cycle
        self.phenopara = cpara['phenopara']
        
        # canopy parameters and state
        self.hc = state['hc'] + epsi
        self.cf = state['cf'] + epsi
        
        #self.cf = 0.1939 * ba / (0.1939 * ba + 1.69) + epsi
        # canopy closure [-] as function of basal area ba m2ha-1;
        # fitted to Korhonen et al. 2007 Silva Fennica Fig.2
        
        spec_para = cpara['spec_para']
        ptypes = {}
        LAI = 0.0

        for pt in list(spec_para.keys()):
            ptypes[pt] = spec_para[pt]
            ptypes[pt]['LAImax'] = state['LAI_' + pt]            

        self.ptypes = ptypes
        
        # compute gridcell average LAI and photosynthesis-stomatal conductance parameters:
        LAI = 0.0
        Amax = 0.0
        q50 = 0.0
        g1 = 0.0
        for pt in self.ptypes.keys():
            if self.ptypes[pt]['lai_cycle']:
                pt_lai = self.ptypes[pt]['LAImax'] * self.phenopara['lai_decid_min']
            else:
                pt_lai = self.ptypes[pt]['LAImax']
            
            LAI += pt_lai
            Amax += pt_lai * ptypes[pt]['amax']
            q50 += pt_lai * ptypes[pt]['q50']
            g1 += pt_lai * ptypes[pt]['g1']
        
        self.LAI = LAI + epsi
        self.physpara.update({'Amax': Amax / self.LAI, 'q50': q50 / self.LAI, 'g1': g1 / self.LAI})

        del Amax, q50, g1, pt, LAI, pt_lai
   
        # - compute start day of senescence: starts at first doy when daylength < self.phenopara['sdl']
        doy = np.arange(1, 366)
        dl = daylength(self.Lat, self.Lon, doy)

        ix = np.max(np.where(dl > self.phenopara['sdl']))
        self.phenopara['sso'] = doy[ix]  # this is onset date for senescence
        del ix
        
        # snow model
        self.wmax = cpara['interc']['wmax']
        self.wmaxsnow = cpara['interc']['wmaxsnow']
        self.Kmelt = cpara['snow']['kmelt']
        self.Kfreeze = cpara['snow']['kfreeze']
        self.R = cpara['snow']['r']  # max fraction of liquid water in snow

        # --- for computing aerodynamic resistances
        self.zmeas = cpara['flow']['zmeas']
        self.zground =cpara['flow']['zground'] # reference height above ground [m]
        self.zo_ground = cpara['flow']['zo_ground'] # ground roughness length [m]
        self.gsoil = self.physpara['gsoil']
        
        # --- state variables
        self.W = np.minimum(state['w'], self.wmax*self.LAI)
        self.SWE = state['swe']
        self.SWEi = self.SWE
        self.SWEl = np.zeros(np.shape(self.SWE))

        # deciduous leaf growth stage
        # NOTE: this assumes simulations start 1st Jan each year !!!
        self.DDsum = 0.0
        self.X = 0.0

        self._relative_lai = self.phenopara['lai_decid_min']
        self._growth_stage = 0.0
        self._senesc_stage = 0.0

        # phenological state
        self.fPheno = self.phenopara['fmin']
        
        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'PotInf': [], 'Trfall': [], 'Interc': [], 'Evap': [],
                            'ET': [], 'Transpi': [], 'Efloor': [], 'SWE': [],
                            'LAI': [], 'Mbe': [], 'LAIfract': [], 'Unload': []

       """
        Initializes BucketGrid:
        Args:
            REQUIRED:
            spara - dictionary of soil properties. keys - values np.arrays
                depth [m]
                poros [m3m-3]
                fc [m3m-3]
                wp [m3m-3]
                ksat [ms-1]
                beta [-]
                maxpond [-]
                org_depth [m]
                org_poros [m3m-3]
                org_fw [m3m-3]
                org_rw [m3m-3]
            
                pond_sto - initial pond storage [m]
                org_sat - initial saturation or organic layer [-]
                rootzone_sat - initial saturation of root zone [-]
            OPTIONAL:  
            outputs - True appends output grids to dict stored within object 

        CHANGES:
            05.05.2020 removed typo in watbal mbe computation and added outputs
        """
        
        """ set object properties. All will be 1d or 2d arrays of same shape """
        # above-ground pond storage [m]
        self.MaxPond = spara['maxpond']
        
        # top layer is interception storage, which capacity depends on its depth [m]
        # and field capacity
        self.D_top = spara['org_depth']     # depth, m3 m-3
        self.poros_top = spara['org_poros'] # porosity, m3 m-3
        self.Fc_top = spara['org_fc']       # field capacity m3 m-3
        self.rw_top = spara['org_rw']       # ree parameter m3 m-3
        self.MaxStoTop = self.Fc_top * self.D_top # maximum storage m 

        # root-zone layer properties
        self.D = spara['depth']             # depth, m
        self.poros = spara['poros']         # porosity, m3 m-3     
        self.Fc = spara['fc']               # field capacity, m3 m-3 
        self.Wp = spara['wp']               # wilting point, m3 m-3
        self.Ksat = spara['ksat']           # sat. hydr. cond., m s-1
        self.beta = spara['beta']           # hyd. cond. exponent, -
        # self.soilcode = spara['soilcode']   # soil type integer code
        
        self.MaxWatSto = self.D*self.poros  # maximum soil water storage, m

        """
        set buckets initial state: given as arrays
        """
        self.PondSto = np.minimum(spara['pond_sto'], self.MaxPond)
        
        # toplayer storage and relative conductance for evaporation
        self.WatStoTop = self.MaxStoTop * spara['org_sat']
        self.Wliq_top = self.poros_top *self.WatStoTop / self.MaxStoTop
        self.Ree = np.minimum(self.Wliq_top / self.rw_top, 1.0) # relative ecaporation rate (-)
        
        # root zone storage and relative extractable water
        self.WatSto = np.minimum(spara['rootzone_sat']*self.D*self.poros, self.D*self.poros)
        
        self.Wliq = self.poros*self.WatSto / self.MaxWatSto
        self.Wair = self.poros - self.Wliq
        self.Sat = self.Wliq/self.poros
        self.Rew = np.minimum((self.Wliq - self.Wp) / (self.Fc - self.Wp + eps), 1.0)
        
        # grid total drainage to ground water [m]
        self._drainage_to_gw = 0.0
        
        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'Infil': [], 'Retflow': [], 'Drain': [], 'Roff': [], 'ET': [],
            'Mbe': [], 'Wliq': [], 'PondSto': [], 'Wliq_top': [], 'Ree': []}

    def __init__(self, pp, cellarea, cmask, flowacc, slope, S_initial=None,
                 outputs=False):
        """
        sets up Topmodel for the catchment assuming homogenous
        effective soil depth 'm' and sat. hydr. conductivity 'ko'.
        This is the 'classic' version of Topmodel where hydrologic similarity\
        index is TWI = log(a / tan(b)).
        
        Args:
            pp - parameter dict with keys:
                dt - timestep [s]
                ko - soil transmissivity at saturation [m/s]
                m -  effective soil depth (m), i.e. decay factor of Ksat with depth
                twi_cutoff - max allowed twi -index
                so - initial catchment average saturation deficit (m)
            cmask - catchment mask, 1 = catchment_cell
            cellarea - gridcell area [m2]
            flowacc - flow accumulation per unit contour length (m)
            slope - local slope (deg)
            S_initial - initial storage deficit, overrides that in 'pp'
            outputs - True stores outputs after each timestep into dictionary
        """
        if not S_initial:
            S_initial = pp['so']

        self.dt = float(pp['dt'])
        self.cmask = cmask
        self.CellArea = cellarea
        dx = cellarea**0.5
        self.CatchmentArea = np.size(cmask[cmask == 1])*self.CellArea

        # topography
        self.a = flowacc*cmask  # flow accumulation grid
        self.slope = slope*cmask  # slope (deg) grid

        # effective soil depth [m]
        self.M = pp['m']
        # lat. hydr. conductivity at surface [m2/timestep]
        # self.To = pp['ko']*pp['m']*self.dt
        self.To = pp['ko']*self.dt
        
        """ 
        local and catchment average hydrologic similarity indices (xi, X).
        Set xi > twi_cutoff equal to cutoff value to remove tail of twi-distribution.
        This concerns mainly the stream network cells. 'Outliers' in twi-distribution are
        problem for streamflow prediction
        """
        slope_rad = np.radians(self.slope)  # deg to rad

        xi = np.log(self.a / dx / (np.tan(slope_rad) + eps))
        # apply cutoff
        clim = np.percentile(xi[xi > 0], pp['twi_cutoff'])
        xi[xi > clim] = clim
        self.xi = xi
  
        self.X = 1.0 / self.CatchmentArea*np.nansum(self.xi*self.CellArea)

        # baseflow rate when catchment Smean=0.0
        self.Qo = self.To*np.exp(-self.X)

        # catchment average saturation deficit S [m] is the only state variable
        s = self.local_s(S_initial)
        s[s < 0] = 0.0
        self.S = np.nanmean(s)

        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'S': [], 'Qb': [], 'Qr': [], 'Qt': [], 'qr': [],
                            'fsat': [], 'Mbe': [], 'R': []
                           }


    ! create netCDF output file
    call create_netcdf()
    
    dlat, dlon = np.shape(spa.GisData['cmask'])

    resultsfile = os.path.join(spa.pgen['results_folder'], spa.pgen['ncf_file'])

    ncf, ncf_file = initialize_netCDF(ID=spa.id, fname=resultsfile,
                                      lat0=spa.GisData['lat0'],
                                      lon0=spa.GisData['lon0'],
                                      dlat=dlat, dlon=dlon, dtime=None)
            
    print('********* created SpaFHy instance *********')
    #print('Loops total [s]: ', timeit.default_timer() - start_time)
    
    if ncf:
        return spa, ncf, ncf_file
    else:
        return spa


  
  END SUBROUTINE initialization_spafhy


  SUBROUTINE initialization_yasso_steadysate (param, flux_leafc_day, flux_rootc_day, flux_nitr_day, tempr_c, precip_day, tempr_ampl, totc_min, cstate, nstate)
  ! A simple algorithm to initialize the SOC pools into a steady state or a partial
  ! steady-state. First, the equilibrium SOC is evaluated. Then, if totc_min is > 0 and
  ! greater than the equilibrium, the deficit will be covered by increasing the H
  ! pool. The nitrogen pool is left unconstrained and is equal to the equilibrium N +
  ! the possible contribution from the extra H. 
    real, intent(in) :: param(:) ! parameter vector
    real, intent(in) :: flux_leafc_day ! carbon input with "leaf" composition per day
    real, intent(in) :: flux_rootc_day ! carbon input with "fineroot" composition per day
    real, intent(in) :: flux_nitr_day  ! organic nitrogen input per day
    real, intent(in) :: tempr_c
    real, intent(in) :: tempr_ampl
    real, intent(in) :: precip_day ! mm
    real, intent(in) :: totc_min ! see above
    real, intent(out) :: cstate(statesize_yasso)
    real, intent(out) :: nstate ! nitrogen
    
    real :: neg_c_input_yr(statesize_yasso)
    real :: matrix(statesize_yasso, statesize_yasso)
    real :: totc
    real :: decomp_h
    real :: cue
    real :: cupt_awen
    real :: nc_awen
    real :: growth_c
    real :: resp
    integer :: cue_iter
    real :: nc_som

    integer, parameter :: max_cue_iter = 10
    
    ! Carbon
    ! 
    neg_c_input_yr = -(flux_leafc_day * awenh_leaf + flux_rootc_day * awenh_fineroot) * 365.0
    call evaluate_matrix_mean_tempr(param, tempr_c, precip_day * days_yr,tempr_ampl,  matrix)
    ! Solve the equilibrium condition Ax + b = 0
    call solve(matrix, neg_c_input_yr, cstate)
    totc = sum(cstate)
    !print *, 'totc before adjust', totc, totc_min
    if (totc_min > 0.0 .and. totc < totc_min) then
       cstate(5) = cstate(5) + totc_min - totc
    end if
    !print *, 'totc after adjust', sum(cstate)

    ! Nitrogen
    !
    decomp_h = matrix(5,5) * cstate(5)
    cue = 0.43 ! initially
    resp = -sum(neg_c_input_yr) ! respiration equal to C input in equilibrium
    do cue_iter = 1, max_cue_iter
       cupt_awen = (resp - decomp_h) / (1.0 - cue)
       growth_c = cue * cupt_awen
       ! Solve nc_awen from the state equation (below) such that nstate becomes stationary:
       nc_awen = (1.0 / cupt_awen) * (nc_mb * cue * cupt_awen - nc_h_max*decomp_h + flux_nitr_day*days_yr)
       nstate = sum(cstate(1:4)) * nc_awen + nc_h_max * cstate(5)
       nc_som = nstate / sum(cstate)
       cue = max(min(0.43 * (nc_som / nc_mb) ** 0.6, 1.0), cue_min)
    end do
  
  END SUBROUTINE initialization_yasso_steadystate

end module initialization_mod
