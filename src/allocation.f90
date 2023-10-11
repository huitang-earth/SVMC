module allocation

implicit none

   type, public :: alloc_para_type
   !----------------------------
   ! allometric parameters
   !-----------------------------
     real(8) :: cratio_resp         ! fraction of respiration to gpp
     real(8) :: cratio_leaf         ! carbon ratio of leaf to npp
     real(8) :: cratio_root         ! carbon ratio of root to npp
     real(8) :: cratio_biomass      ! carbon ratio of biomass
     real(8) :: harvest_index       !
     real(8) :: turnover_cleaf      !
     real(8) :: turnover_croot      !
     real(8) :: sla                 ! specific leaf area, 

   end type alloc_para_type

   type, public :: management_para_type
   !----------------------------
   ! management parameters
   !-----------------------------
     real(8) :: cut_ratio          ! cutting ratio
     real(8) :: graze_ratio        ! grazing ratio
   end type management_para_type

   public alloc_hypothesis_1    ! estimate litter input from gpp directly, fixed ratio
   public alloc_hypothesis_2    ! estimate litter input from npp -> above- & below- ground biomass, only management 
!   public alloc_hypothesis_3    ! same to 2, but with continuous litter input?
!   public alloc_hypothesis_4    ! variable allometric co-efficient?
!   public alloc_hypothesis_5    ! variable litter qualitity co-efficient?
!   public cal_lai_from_leafc

contains

   subroutine alloc_hypothesis_1(gpp_day, npp_day, litter_cleaf, litter_croot, alloc_para)
      real(8), intent(in)    :: gpp_day      ! gpp (daily average),  kg C m-2 s-1
      real(8), intent(inout) :: npp_day      ! npp (daily average),  kg C m-2 s-1
      real(8), intent(inout) :: litter_cleaf ! carbon input with "leaf" composition per day
      real(8), intent(inout) :: litter_croot ! carbon input with "root" composition per day
      type(alloc_para_type), intent(inout) :: alloc_para   ! allometric parameters
      
      alloc_para%cratio_resp = 0.5
      alloc_para%cratio_leaf =0.5
      alloc_para%harvest_index = 0.5 

      litter_cleaf= gpp_day * (1-alloc_para%cratio_resp) * (1-alloc_para%harvest_index) * &
                               alloc_para%cratio_leaf * 3600 * 24
      litter_croot= gpp_day * (1-alloc_para%cratio_resp) * (1-alloc_para%harvest_index) * &
                              (1-alloc_para%cratio_leaf) * 3600 * 24

   end subroutine alloc_hypothesis_1

   subroutine alloc_hypothesis_2(gpp_day, npp_day, auto_resp, croot, cleaf, cstem, litter_cleaf, litter_croot, &
                                  abovebiomass, belowbiomass, yield, &
                                  lai, alloc_para, manage_para, pheno_stage, management_type)

      real(8), intent(in)    :: gpp_day      ! gpp (daily average),  kg C m-2 s-1
      real(8), intent(inout) :: npp_day      ! npp (daily average),  kg C m-2 s-1
      real(8), intent(inout) :: auto_resp      ! npp (daily average),  kg C m-2 s-1
      real(8), intent(inout) :: croot        ! root carbon
      real(8), intent(inout) :: cleaf        ! leaf carbon
      real(8), intent(inout) :: cstem        ! leaf carbon
      real, intent(inout) :: litter_cleaf ! carbon input with "leaf" composition per day
      real, intent(inout) :: litter_croot ! carbon input with "root" composition per day
      real(8), intent(inout) :: lai          ! leaf area index
      real(8), intent(inout) :: abovebiomass ! Abovegroud biomass (kg dry mass m-2)
      real(8), intent(inout) :: belowbiomass ! Abovegroud biomass (kg dry mass m-2)
      real(8), intent(inout) :: yield        ! yield
      type(alloc_para_type), intent(inout) :: alloc_para   ! allometric parameters
      type(management_para_type), intent(inout) :: manage_para   ! allometric parameters
      integer, intent(in)    :: pheno_stage, management_type

      ! local
      real(8) :: litter_cstem ! carbon input with "stem" composition per day

      alloc_para%cratio_resp    = 0.5
      alloc_para%cratio_leaf    = 0.4    ! later season grass, more in root and straw
      alloc_para%cratio_root    = 0.5
      alloc_para%harvest_index  = 0.5
      alloc_para%cratio_biomass = 0.45
      alloc_para%turnover_cleaf = 0.41/365 
      alloc_para%turnover_croot = 0.41/365  ! depend on phenological stage? scaled to NPP:NPPmax
      alloc_para%sla            = 30    ! m2 kg-1         
                                        ! leaf area in cm2 produced g−1 leaf dry weight plant−1 (500)

      if (pheno_stage .eq. 1) then    ! phenology more critic for cereal, sow, emergence, maturity, flower, grainfill
                                      ! temperature dependence can be good ...
                                      ! storage carbon pool?
                                      ! Use reversed LAI (remote sensed) to leafc carbon to estimate litter?

         if (management_type .eq. 0) then    ! no management, organic fertilizer, potential yields, Nitrogen (?) 
           npp_day = gpp_day * (1-alloc_para%cratio_resp) * 3600 * 24
           auto_resp = gpp_day * alloc_para%cratio_resp * 3600 * 24
           litter_cleaf= cleaf * alloc_para%turnover_cleaf
           litter_cstem= cstem * alloc_para%turnover_cleaf
           litter_croot= croot * alloc_para%turnover_croot
           cleaf   = cleaf + npp_day * alloc_para%cratio_leaf - litter_cleaf
           cstem   = cstem + npp_day * (1-alloc_para%cratio_leaf-alloc_para%cratio_root) - litter_cstem
           croot   = croot + npp_day * alloc_para%cratio_root - litter_croot
         end if

         if (management_type .eq. 1) then    ! cutting,  grass is special....
            manage_para%cut_ratio = 0.8            
            cleaf = (1-manage_para%cut_ratio) * cleaf
            cstem = (1-manage_para%cut_ratio) * cstem
            ! no litter input to soil
         end if

         if (management_type .eq. 2) then    ! grazing
            manage_para%graze_ratio = 0.1
            cleaf = (1-manage_para%graze_ratio) * cleaf
            cstem = (1-manage_para%graze_ratio) * cstem
            ! manure input
         end if

                                               

         litter_cleaf= litter_cleaf + litter_cstem ! combine leaf and stem together
      end if

      abovebiomass = (cleaf+cstem)/alloc_para%cratio_biomass   ! dry matter, 
                                                               ! To compare with observation need to convert to wet biomass
                                                               ! (divided by DM) and grain loss (LO) due to technical reasons.
      belowbiomass = croot/alloc_para%cratio_biomass

      ! yield        = abovebiomass * alloc_para%harvest_index                     
      lai          = cleaf/alloc_para%cratio_biomass * alloc_para%sla

      !call cal_lai_from_leafc(lai, leafc)

   end subroutine alloc_hypothesis_2


  !subroutine alloc_hypothesis_1(gpp_day, npp_day, croot, cleaf, litter_cleaf, litter_croot, lai, alloc_para)
  !    real(8), intent(in)    :: gpp_day      ! parameter vector
  !    real(8), intent(inout) :: npp_day      ! carbon input with "leaf" composition per day
  !    real(8), intent(inout) :: croot        ! carbon input with "leaf" composition per day
  !    real(8), intent(inout) :: cleaf        ! carbon input with "leaf" composition per day
  !    real(8), intent(inout) :: litter_cleaf ! carbon input with "leaf" composition per day
  !    real(8), intent(inout) :: litter_croot ! carbon input with "leaf" composition per day
  !    real(8), intent(inout) :: lai          ! carbon input with "leaf" composition per day
  !    real(8), intent(inout) :: alloc_para   ! carbon input with "leaf" composition per day    
      
  !    alloc_para%cratio_resp = 0.5
  !   alloc_para%cratio_leaf =0.5
  !    alloc_para%harvest_index = 0.5 

  !    litter_cleaf= gpp_day * (1-alloc_para%cratio_resp) * (1-alloc_para%harvest_index) * &
  !                             alloc_para%cratio_leaf * 3600 * 24
  !    litter_croot= gpp_day * (1-alloc_para%cratio_resp) * (1-alloc_para%harvest_index) * &
  !                            (1-alloc_para%cratio_leaf) * 3600 * 24

  ! end subroutine alloc_hypothesis_1

end module allocation

