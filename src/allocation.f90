module allocation

implicit none

   type, public :: alloc_para_type
   !----------------------------
   ! allometric parameters
   !-----------------------------
     real(8) :: cratio_resp         ! fraction of maintenance respiration (1/m2/s) at 20 degree
     real(8) :: cratio_leaf         ! carbon ratio of leaf to npp
     real(8) :: cratio_root         ! carbon ratio of root to npp
     real(8) :: cratio_biomass      ! carbon ratio of biomass
     real(8) :: harvest_index       !
     real(8) :: turnover_cleaf      ! turnover rate of leaf at 20 degree
     real(8) :: turnover_croot      ! turnover rate of root at 20 degree
     real(8) :: sla                 ! specific leaf area, 
     real(8) :: q10                 ! Q10 temperature coefficient (https://en.wikipedia.org/wiki/Q10_(temperature_coefficient))              

   end type alloc_para_type

   type, public :: management_data_type
   !----------------------------
   ! management parameters
   !-----------------------------
     integer :: management_type          ! management types
     real(8) :: management_c_input       ! 
     real(8) :: management_c_output      ! 
     real(8) :: management_n_input       ! 
     real(8) :: management_n_output      ! 
   end type management_data_type

   public alloc_hypothesis_1    ! estimate litter input from gpp directly, fixed ratio
   public alloc_hypothesis_2    ! estimate litter input from npp -> above- & below- ground biomass, only management 
!   public alloc_hypothesis_3    ! same to 2, but with continuous litter input?
!   public alloc_hypothesis_4    ! variable allometric co-efficient?
!   public alloc_hypothesis_5    ! variable litter qualitity co-efficient?
!   public cal_lai_from_leafc

contains

   subroutine readalloc_namelist(alloc_para)
  
      type(alloc_para_type), intent(inout) :: alloc_para
    
      ! Local variables
      real(8) :: cratio_resp         ! fraction of respiration to gpp
      real(8) :: cratio_leaf         ! carbon ratio of leaf to npp
      real(8) :: cratio_root         ! carbon ratio of root to npp
      real(8) :: cratio_biomass      ! carbon ratio of biomass
      real(8) :: harvest_index       !
      real(8) :: turnover_cleaf      !
      real(8) :: turnover_croot      !
      real(8) :: sla                 ! specific leaf area, 
      real(8) :: q10                 ! Q10 temperature coefficient 

      logical :: old
      integer :: readerror
      integer,parameter :: unitallocpara=5

      namelist /alloc_namelist/ &
       cratio_resp, &
       cratio_leaf, &
       cratio_root, &
       cratio_biomass, &
       harvest_index, &
       turnover_cleaf, &
       turnover_croot, &
       sla, &
       q10

      old=.false.

      ! Default setting of allometric parameters
      cratio_resp    = 0.4
      cratio_leaf    = 0.8        !0.4    ! later season grass, more in root and straw
      cratio_root    = 0.2
      cratio_biomass = 0.42
      harvest_index  = 0.5
      turnover_cleaf = 0.41/365 
      turnover_croot = 0.41/365  ! depend on phenological stage? scaled to NPP:NPPmax
      sla            = 10    ! m2 kg-1     
                             ! leaf area in cm2 produced g−1 leaf dry weight plant−1 (500)
                             ! derived from https://en.wikipedia.org/wiki/Specific_leaf_area
      q10            = 1     ! Commonly used in models

      ! Reading namelist
      open(unitallocpara, file='./alloc_namelist', status='old', form='formatted', err=999)
      read(unitallocpara,alloc_namelist,iostat=readerror)
      close(unitallocpara)

      print *, "harvest_index=", harvest_index

      alloc_para%cratio_resp    = cratio_resp
      alloc_para%cratio_leaf    = cratio_leaf
      alloc_para%cratio_root    = cratio_root 
      alloc_para%harvest_index  = harvest_index
      alloc_para%cratio_biomass = cratio_biomass
      alloc_para%turnover_cleaf = turnover_cleaf 
      alloc_para%turnover_croot = turnover_croot
      alloc_para%sla            = sla
      alloc_para%q10            = q10

      return

999   write(*,*) ' #### MODEL ERROR! FILE "alloc_namelist"    #### '
      write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
      stop

   end subroutine readalloc_namelist

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

   subroutine alloc_hypothesis_2(temp_day, gpp_day, npp_day, leaf_rdark_day, auto_resp, croot, cleaf, cstem, litter_cleaf, litter_croot, &
                                  compost, abovebiomass, belowbiomass, yield, &
                                  lai, alloc_para, manage_data, pheno_stage)

      real(8), intent(in)    :: temp_day     ! temperature (exponentially averaged),  celcius degree
      real(8), intent(in)    :: gpp_day      ! gpp (daily average),  kg C m-2 s-1
      real(8), intent(in)    :: leaf_rdark_day    ! leaf dark respiration (daily average), kg C m-2 s-1
      real(8), intent(inout) :: npp_day      ! npp (daily average),  kg C m-2 s-1
      real(8), intent(inout) :: auto_resp      ! npp (daily average),  kg C m-2 s-1
      real(8), intent(inout) :: croot        ! root carbon
      real(8), intent(inout) :: cleaf        ! leaf carbon
      real(8), intent(inout) :: cstem        ! leaf carbon
      real, intent(inout) :: litter_cleaf ! carbon input with "leaf" composition per day
      real, intent(inout) :: litter_croot ! carbon input with "root" composition per day
      real, intent(inout) :: compost      ! manure input (kg C m-2 s-1)
      real(8), intent(inout) :: lai          ! leaf area index, allometric
      real(8), intent(inout) :: abovebiomass ! Abovegroud biomass (kg dry mass m-2)
      real(8), intent(inout) :: belowbiomass ! Abovegroud biomass (kg dry mass m-2)
      real(8), intent(inout) :: yield        ! yield
      type(alloc_para_type), intent(inout) :: alloc_para   ! allometric parameters
      type(management_data_type), intent(inout) :: manage_data   ! allometric parameters
      integer, intent(in)    :: pheno_stage

      ! local
      real(8) :: litter_cstem ! carbon input with "stem" composition per day

      if (pheno_stage .eq. 1) then    ! phenology more critic for cereal, sow, emergence, maturity, flower, grainfill
                                      ! temperature dependence can be good ...
                                      ! storage carbon pool?
                                      ! Use reversed LAI (remote sensed) to leafc carbon to estimate litter?

         !if (manage_data%management_type .eq. 0) then    ! no management, organic fertilizer, potential yields, Nitrogen (?) 
         
         ! Allow maintenance respiration of root to be calculated from root carbon storage, not gpp.
           npp_day = (gpp_day - (croot + cstem) * (alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10)) &
                              - leaf_rdark_day )* 3600 * 24
           auto_resp = ((croot + cstem) * (alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10)) &
                              + leaf_rdark_day )* 3600 * 24

           litter_cleaf=cleaf * alloc_para%turnover_cleaf * alloc_para%q10 ** ((temp_day  - 20)/10)
           litter_cstem=cstem * alloc_para%turnover_cleaf * alloc_para%q10 ** ((temp_day  - 20)/10)
           litter_croot=croot * alloc_para%turnover_croot * alloc_para%q10 ** ((temp_day  - 20)/10)
           compost=0.0
           cleaf   = cleaf + npp_day * alloc_para%cratio_leaf - litter_cleaf
           cstem   = cstem + npp_day * (1-alloc_para%cratio_leaf-alloc_para%cratio_root) - litter_cstem
           croot   = croot + npp_day * alloc_para%cratio_root - litter_croot
         !end if

         if (manage_data%management_type .eq. 1) then    ! harvesting,  grass is special....           
            cleaf = cleaf - manage_data%management_c_output*3600*24*cleaf/(cleaf+cstem)
            cstem = cstem - manage_data%management_c_output*3600*24*cstem/(cleaf+cstem)
            ! no litter input to soil
         else if (manage_data%management_type .eq. 3) then    ! grazing
            cleaf = cleaf - manage_data%management_c_output*3600*24*cleaf/(cleaf+cstem)
            cstem = cstem - manage_data%management_c_output*3600*24*cstem/(cleaf+cstem)
            ! manure input
            compost= manage_data%management_c_input*3600*24
         else if (manage_data%management_type .eq. 4) then    ! organic materials
            litter_cleaf = litter_cleaf+manage_data%management_c_input*3600*24
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

