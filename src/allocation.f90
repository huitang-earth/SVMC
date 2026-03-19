module allocation

use readvegpara_mod

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
     real(8) :: invert_option       ! 0: no inversion from LAI
                                    ! 1: inversion for cratio_leaf
                                    ! 2: inversion for turnover_leaf           

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
      integer :: invert_option

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
       q10, &
       invert_option

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
      invert_option  = 0     ! no inversion

      ! Reading namelist
      open(unitallocpara, file='./alloc_namelist', status='old', form='formatted', err=999)
      read(unitallocpara,alloc_namelist,iostat=readerror)
      close(unitallocpara)

      print *, "invert_option=", invert_option

      alloc_para%cratio_resp    = cratio_resp
      alloc_para%cratio_leaf    = cratio_leaf
      alloc_para%cratio_root    = cratio_root 
      alloc_para%harvest_index  = harvest_index
      alloc_para%cratio_biomass = cratio_biomass
      alloc_para%turnover_cleaf = turnover_cleaf 
      alloc_para%turnover_croot = turnover_croot
      alloc_para%sla            = sla
      alloc_para%q10            = q10
      alloc_para%invert_option  = invert_option

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

   subroutine alloc_hypothesis_2(temp_day, gpp_day, npp_day, leaf_rdark_day, auto_resp, croot, cleaf, cstem, cgrain, & 
                                  litter_cleaf, litter_croot, compost, abovebiomass, belowbiomass, yield, &
                                  lai, alloc_para, grain_fill, manage_data, pheno_stage)

      real(8), intent(in)    :: temp_day     ! temperature (exponentially averaged),  celcius degree
      real(8), intent(in)    :: gpp_day      ! gpp (daily average),  kg C m-2 s-1
      real(8), intent(in)    :: leaf_rdark_day    ! leaf dark respiration (daily average), kg C m-2 s-1

      real(8), intent(inout) :: npp_day      ! npp (daily average),  kg C m-2 s-1
      real(8), intent(inout) :: auto_resp      ! npp (daily average),  kg C m-2 s-1
      real(8), intent(inout) :: croot        ! root carbon
      real(8), intent(inout) :: cleaf        ! leaf carbon
      real(8), intent(inout) :: cstem        ! leaf carbon
      real(8), intent(inout) :: cgrain       ! grain carbon
      real, intent(inout) :: litter_cleaf ! carbon input with "leaf" composition per day
      real, intent(inout) :: litter_croot ! carbon input with "root" composition per day
      real, intent(inout) :: compost      ! manure input (kg C m-2 s-1)
      real(8), intent(inout) :: lai          ! leaf area index, allometric
      real(8), intent(inout) :: abovebiomass ! Abovegroud biomass (kg dry mass m-2)
      real(8), intent(inout) :: belowbiomass ! Abovegroud biomass (kg dry mass m-2)
      real(8), intent(inout) :: yield        ! yield
      real(8), intent(inout) :: grain_fill   ! grain-filling C flux, kg C m-2 s-1
      type(alloc_para_type), intent(inout) :: alloc_para   ! allometric parameters
      type(management_data_type), intent(inout) :: manage_data   ! allometric parameters
      integer, intent(inout)    :: pheno_stage

      ! local
      real(8) :: litter_cstem                                  ! carbon input with "stem" composition per day
      real(8) :: gr_resp_leaf, gr_resp_stem, gr_resp_root, gr_resp_grain      ! growth respiration, kg C m-2 s-1 


      if (pheno_stage .eq. 1) then    ! phenology more critic for cereal, sow, emergence, maturity, flower, grainfill
                                      ! temperature dependence can be good ...
                                      ! storage carbon pool?
                                      ! Use reversed LAI (remote sensed) to leafc carbon to estimate litter?

         !if (manage_data%management_type .eq. 0) then    ! no management, organic fertilizer, potential yields, Nitrogen (?) 
         
         ! Allow maintenance respiration of root to be calculated from root carbon storage, not gpp.
         gr_resp_leaf=max(0.0, (gpp_day * alloc_para%cratio_leaf - leaf_rdark_day)*3600*24)*0.11
         gr_resp_stem=max((gpp_day * (1-alloc_para%cratio_leaf-alloc_para%cratio_root) &
                                  - cstem * (alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10)))*3600*24,0.0)*0.11
         gr_resp_root=max((gpp_day * alloc_para%cratio_root - grain_fill  &
                                  - croot * (alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10)))*3600*24,0.0)*0.11
         if (pft_type=="oat") then
            gr_resp_grain=max((grain_fill - cgrain*(0.1*alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10)))*3600*24, &
                              0.0)*0.11
         else
            gr_resp_grain=0.0 
         end if        

         npp_day = (gpp_day - (croot + cstem) * (alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10)) &
                            - cgrain * (0.1 * alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10)) &
                            - leaf_rdark_day )* 3600 * 24
         auto_resp = ((croot + cstem) * (alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10)) &
                            + cgrain * (0.1 * alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10)) &
                            + leaf_rdark_day )* 3600 * 24  &
                            + gr_resp_leaf &
                            + gr_resp_stem &
                            + gr_resp_root &
                            + gr_resp_grain

         if (alloc_para%invert_option .eq. 0) then
           litter_cleaf=cleaf * alloc_para%turnover_cleaf * alloc_para%q10 ** ((temp_day  - 20)/10)
         end if
         litter_cstem=cstem * alloc_para%turnover_croot * alloc_para%q10 ** ((temp_day  - 20)/10)
         litter_croot=croot * alloc_para%turnover_croot * alloc_para%q10 ** ((temp_day  - 20)/10)
         compost=0.0
         
         if (alloc_para%invert_option .eq. 0) then
            cleaf   = max(0.001, cleaf + gpp_day * alloc_para%cratio_leaf * 3600 * 24 &  
                               - litter_cleaf - leaf_rdark_day*3600*24 - gr_resp_leaf)
         end if
         
         ! Need to consider if cratio_leaf is too large for invert_option: 1 & 2
         cstem   = max(0.0, cstem + gpp_day * (1-alloc_para%cratio_leaf-alloc_para%cratio_root) *3600 *24 &
                     - cstem * (alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10))*3600*24 &
                     - litter_cstem &
                     - gr_resp_stem)
         croot   = max(0.0, croot + (gpp_day * alloc_para%cratio_root - grain_fill) * 3600 * 24 - litter_croot &
                     - croot * (alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10))*3600*24 &
                     - gr_resp_root)
         
         if (pft_type=="oat") then
             cgrain  = cgrain + grain_fill*3600*24  &
                     - cgrain * (0.1*alloc_para%cratio_resp * alloc_para%q10 ** ((temp_day - 20)/10))*3600*24 &
                     - gr_resp_grain
         else
            cgrain=0.0
         end if
         !end if

         if (manage_data%management_type .eq. 1) then    ! harvesting,  grass is special....        
            if (pft_type=="grass") then 
               if (alloc_para%invert_option .eq. 0) then    
                  cleaf = max(0.001, cleaf - manage_data%management_c_output*3600*24*cleaf/(cleaf+cstem))
               end if
               cstem = max(0.0, cstem - manage_data%management_c_output*3600*24*cstem/(cleaf+cstem))
               print *, "cstem=", cstem
               print *, "cleaf=", cleaf
      
               ! no litter input to soil for perennial forage grass, the plant parts remain alive.  
            else if (pft_type=="oat") then
               ! For cereal crop, the harvest yield is for grain carbon pool, which is separated from root (not leaf) in the model.  
               !if (alloc_para%invert_option .eq. 1) then 
                  ! Need to build grain C pool here
                  !croot = croot - manage_data%management_c_output*3600*24
                  ! After harvest, the remaining leaf, stem & root carbon go to litter carbon pools. 
                  litter_croot=croot
                  litter_cleaf=cleaf
                  litter_cstem=cstem
                  croot=0.0
                  cleaf=0.0
                  cstem=0.0
                  cgrain=0.0
                  npp_day=0.0
                  auto_resp=0.0
               !end if
            end if

         else if (manage_data%management_type .eq. 3) then    ! grazing
            if (alloc_para%invert_option .eq. 0) then
              cleaf = max(0.001, cleaf - manage_data%management_c_output*3600*24*cleaf/(cleaf+cstem))
            end if
            cstem = max(0.0, cstem - manage_data%management_c_output*3600*24*cstem/(cleaf+cstem))
            ! manure input from animal dung or urine.
            compost= manage_data%management_c_input*3600*24
         else if (manage_data%management_type .eq. 4) then    ! organic materials
            compost= manage_data%management_c_input*3600*24
            !litter_cleaf = litter_cleaf+manage_data%management_c_input*3600*24
         end if
                                               
         litter_cleaf= litter_cleaf + litter_cstem            ! combine leaf and stem together
      
      else if (pheno_stage .eq. 2) then
         ! needed only for dynamic LAI?
         ! Build a grain C pool here.                   
         litter_croot=croot
         litter_cleaf=cleaf
         litter_cstem=cstem
         croot=0.0
         cleaf=0.0
         cstem=0.0
         npp_day=0.0
         auto_resp=0.0
         pheno_stage=1  ! Get back to normal phenology after removing all the living biomass,
                        ! e.g. cover crop. If there is no harvest event for cover crop,
                        ! they will be moved to litter c pools in the spring of next year. 
      end if

      abovebiomass = (cleaf+cstem+cgrain)/alloc_para%cratio_biomass   ! dry matter, 
                                                               ! To compare with observation need to convert to wet biomass
                                                               ! (divided by DM) and grain loss (LO) due to technical reasons.
      belowbiomass = croot/alloc_para%cratio_biomass

      ! yield        = abovebiomass * alloc_para%harvest_index                     
      lai          = cleaf/alloc_para%cratio_biomass * alloc_para%sla
      
      !call cal_lai_from_leafc(lai, leafc)

   end subroutine alloc_hypothesis_2

   subroutine invert_alloc(delta_lai, alloc_para, leaf_rdark_day, temp_day, litter_cleaf, gpp_day, cleaf, &
                            cstem, manage_data, pheno_stage)

   ! This is for deriving allometric parameters: turnover_cleaf, cratio_leaf

      real(8), intent(in)    :: temp_day     ! temperature (exponentially averaged),  celcius degree
      real(8), intent(in)    :: gpp_day      ! gpp (daily average),  kg C m-2 s-1
      real(8), intent(in)    :: leaf_rdark_day   ! leaf dark respiration (daily average), kg C m-2 s-1
      real(8), intent(inout) :: cleaf        ! leaf carbon
      real(8), intent(inout) :: cstem        ! leaf carbon
      real(8), intent(inout) :: litter_cleaf ! carbon input with "leaf" composition per day
      real(8), intent(inout) :: delta_lai    ! leaf area index, allometric     
      type(alloc_para_type), intent(inout) :: alloc_para   ! allometric parameters
      type(management_data_type), intent(inout) :: manage_data   ! allometric parameters
      integer, intent(in)    :: pheno_stage

      ! local
      real(8) :: delta_cleaf      ! change of leaf carbon storage required for lai changes.
      real(8) :: gr_resp_leaf     ! growth respiration of leaf.

      if ( pheno_stage .eq. 1) then

         delta_cleaf = delta_lai/alloc_para%sla * alloc_para%cratio_biomass

         !HT: Here cratio_leaf is from the previous time step for simplicity.
         gr_resp_leaf=max(0.0, (gpp_day * alloc_para%cratio_leaf - leaf_rdark_day)*3600*24)*0.11

         ! Need to consider the managment here?
         if (alloc_para%invert_option .eq. 1) then
            litter_cleaf = cleaf * alloc_para%turnover_cleaf * alloc_para%q10 ** ((temp_day  - 20)/10)
            if (gpp_day .gt. 0.2e-8) then
               if (manage_data%management_type .eq. 1) then    ! harvesting,  grass is special....  
                  if (pft_type=="grass") then
                     alloc_para%cratio_leaf = (delta_cleaf + litter_cleaf + leaf_rdark_day * 3600*24 + gr_resp_leaf & 
                          + manage_data%management_c_output*3600*24*cleaf/(cleaf+cstem))/3600/24/gpp_day
                  else
                     alloc_para%cratio_leaf = (delta_cleaf + litter_cleaf + leaf_rdark_day * 3600*24 + gr_resp_leaf &
                           )/3600/24/gpp_day
                  end if

               else if (manage_data%management_type .eq. 3) then    ! grazing
                  if (pft_type=="grass") then
                     alloc_para%cratio_leaf = (delta_cleaf + litter_cleaf + leaf_rdark_day*3600*24 + gr_resp_leaf & 
                        + manage_data%management_c_output*3600*24*cleaf/(cleaf+cstem))/3600/24/gpp_day
                  else
                     alloc_para%cratio_leaf = (delta_cleaf + litter_cleaf + leaf_rdark_day*3600*24 + gr_resp_leaf &  
                          )/3600/24/gpp_day
                  end if 

               else
                  alloc_para%cratio_leaf = (delta_cleaf + litter_cleaf + leaf_rdark_day*3600*24 + gr_resp_leaf)/3600/24/gpp_day
               end if           
               
               alloc_para%cratio_leaf = min(0.9, max(0.1, alloc_para%cratio_leaf))
               alloc_para%cratio_root = 1-alloc_para%cratio_leaf
            end if
            cleaf=max(0.0, cleaf + delta_cleaf)

         else if (alloc_para%invert_option .eq. 2) then
            if (cleaf .gt. 0.00001) then
               if (manage_data%management_type .eq. 1) then   ! harvesting,  grass is special.... 
                  if (pft_type=="grass") then
                     alloc_para%turnover_cleaf = (gpp_day * alloc_para%cratio_leaf * 3600 * 24 - delta_cleaf - &
                           manage_data%management_c_output*3600*24*cleaf/(cleaf+cstem) - leaf_rdark_day*3600*24 - gr_resp_leaf)/ &
                              cleaf/(alloc_para%q10 ** ((temp_day  - 20)/10))
                  else
                     alloc_para%turnover_cleaf = (gpp_day * alloc_para%cratio_leaf * 3600 * 24 - delta_cleaf - &
                           leaf_rdark_day*3600*24 - gr_resp_leaf)/ &
                              cleaf/(alloc_para%q10 ** ((temp_day  - 20)/10))
                  end if

               else if (manage_data%management_type .eq. 3) then    ! grazing
                  if (pft_type=="grass") then
                      alloc_para%turnover_cleaf = (gpp_day * alloc_para%cratio_leaf * 3600 * 24 - delta_cleaf - &
                        manage_data%management_c_output*3600*24*cleaf/(cleaf+cstem) - leaf_rdark_day*3600*24 - gr_resp_leaf)/ &
                           cleaf/(alloc_para%q10 ** ((temp_day  - 20)/10))
                  else
                      alloc_para%turnover_cleaf = (gpp_day * alloc_para%cratio_leaf * 3600 * 24 - delta_cleaf - &
                           leaf_rdark_day*3600*24 - gr_resp_leaf)/ &
                              cleaf/(alloc_para%q10 ** ((temp_day  - 20)/10))
                  end if

               else
                  alloc_para%turnover_cleaf = (gpp_day * alloc_para%cratio_leaf * 3600 * 24 - delta_cleaf - &
                         leaf_rdark_day*3600*24 - gr_resp_leaf)/ &
                           cleaf/(alloc_para%q10 ** ((temp_day  - 20)/10))
               end if

            else
               cleaf=0.0
            end if         
            litter_cleaf= cleaf * alloc_para%turnover_cleaf * alloc_para%q10 ** ((temp_day  - 20)/10)
            cleaf=max(0.0,cleaf + delta_cleaf)
         end if
      
      end if

   end subroutine invert_alloc

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
  !    alloc_para%cratio_leaf =0.5
  !    alloc_para%harvest_index = 0.5 

  !    litter_cleaf= gpp_day * (1-alloc_para%cratio_resp) * (1-alloc_para%harvest_index) * &
  !                             alloc_para%cratio_leaf * 3600 * 24
  !    litter_croot= gpp_day * (1-alloc_para%cratio_resp) * (1-alloc_para%harvest_index) * &
  !                            (1-alloc_para%cratio_leaf) * 3600 * 24

  ! end subroutine alloc_hypothesis_1

end module allocation
