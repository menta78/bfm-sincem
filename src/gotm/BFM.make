#
# Makefile to build the BFM model inside GOTM
#
LIB     = $(LIBDIR)/libbio$(buildtype).a

# BFMDIR path is an environment variable
ifndef BFMDIR
  $(error The environment variable BFMDIR is not defined!)
endif
BFMSRC 	   = $(BFMDIR)/src/BFM
BFMGOTMSRC = $(BFMDIR)/src/gotm
BFMSHARE   = $(BFMDIR)/src/share

# Pelagic flags
INCLUDE_PELCO2=false
ifeq ($(INCLUDE_PELCO2),true)
  DEFINES += -DINCLUDE_PELCO2
endif
INCLUDE_SILT  = true
ifeq ($(INCLUDE_SILT),true)
  DEFINES += -DINCLUDE_SILT
endif

# Benthic ecosystem flags (activates compilation and macros, true by default)
INCLUDE_BEN = true
INCLUDE_BENCO2=false
INCLUDE_BENPROFILES=true
ifeq ($(INCLUDE_BEN),true)
  DEFINES += -DINCLUDE_BEN
  ifeq ($(INCLUDE_BENCO2),true)
    DEFINES += -DINCLUDE_BENCO2
  endif
  ifeq ($(INCLUDE_BENPROFILES),true)
    DEFINES += -DINCLUDE_BENPROFILES
  endif
endif

# this part is still experimental
# BFM subdirectories path (must be one single line)
#BFMSUBDIRS = ${BFMSRC}/Ben:${BFMSRC}/Bennut:${BFMSRC}/Forcing:${BFMSRC}/General:${BFMSRC}/include:${BFMSRC}/Light:${BFMSRC}/Oxygen:${BFMSRC}/PelB:${BFMSRC}/PelBen:${BFMGOTMSRC}:${BFMSHARE} 
# Automatic depend list (needs makedepf90)
#include ${BFMGOTMSRC}/.depend 

OBJ   = \
${LIB}(${BFMGOTMSRC}/bio_var.o)			\
${LIB}(bio_fluxes.o)				\
${LIB}(bio_template.o)				\
${LIB}(bio_npzd.o)				\
${LIB}(bio_iow.o)				\
${LIB}(bio_mab.o)				\
${LIB}(bio_fasham.o)				\
${LIB}(bio_sed.o)				\
${LIB}(${BFMSHARE}/string_functions.o)		\
${LIB}(${BFMSHARE}/init_cnps.o)			\
${LIB}(${BFMSHARE}/ResetFluxes.o)		\
${LIB}(${BFMGOTMSRC}/bio_bfm.o)			\
${LIB}(${BFMGOTMSRC}/trace_bdy.o)		\
${LIB}(${BFMSHARE}/init_var_bfm.o)	        \
${LIB}(${BFMSHARE}/init_benthic_bfm.o)	        \
${LIB}($(BFMSHARE)/make_flux_output.o)		\
${LIB}(${BFMGOTMSRC}/bfm_solver.o)		\
${LIB}(${BFMGOTMSRC}/bio_save.o)		\
${LIB}(${BFMGOTMSRC}/bio_save_bfm.o)		\
${LIB}(${BFMGOTMSRC}/prepare_bio_output.o)	\
${LIB}(${BFMGOTMSRC}/bio.o)

# gotm_error_msg needs to be here until dependencies are fixed
BFM_MOD = \
	${LIB}(${BFMGOTMSRC}/gotm_error_msg.o)				\
	${LIB}(${BFMSRC}/General/ModuleGlobalMem.o)			\
	${LIB}(${BFMSRC}/General/ModuleConstants.o)			\
	${LIB}(${BFMSRC}/General/ModuleGlobFun.o)			\
	${LIB}(${BFMSRC}/General/bfm_error_msg.o)			\
	${LIB}(${BFMSRC}/General/ModuleMem.o)				\
	${LIB}(${BFMSRC}/General/ModuleParam.o)				\
	${LIB}(${BFMSRC}/General/ModuleInterface.o)			\
	${LIB}(${BFMSRC}/Light/ModuleLightAdaptation.o)			\
	${LIB}(${BFMSRC}/Light/ModulePhotoAvailableRadiation.o)		\
	${LIB}(${BFMSRC}/Oxygen/ModuleWindOxReaeration_3.o)		\
	${LIB}(${BFMSRC}/PelB/ModuleMesoZoo.o)				\
	${LIB}(${BFMSRC}/PelB/ModuleMicroZoo.o)				\
	${LIB}(${BFMSRC}/PelB/ModulePelBac.o)				\
	${LIB}(${BFMSRC}/PelB/ModulePelChem.o)				\
	${LIB}(${BFMSRC}/PelB/ModulePelGlobal.o)			\
	${LIB}(${BFMSRC}/PelB/ModulePhyto.o)				\
	${LIB}(${BFMSRC}/PelBen/ModuleSettling.o)			\
	${LIB}(${BFMSRC}/PelBen/ModuleControlBennutBuffers.o)	\
	${LIB}(${BFMSRC}/Ben/ModuleBenOrganism.o)			\
	${LIB}(${BFMSRC}/Ben/ModuleBenBac.o)				\
	${LIB}(${BFMSRC}/Ben/ModuleFilterFeeder.o)			\
	${LIB}(${BFMSRC}/Ben/ModuleBioturbation.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicNutrient3.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenAmmonium.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenAnoxic.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenDenitriDepth.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNitrate.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutConstants.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutType.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutInterface.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutVariables.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenOxygen.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenPhosphate.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenQ1Transport.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenSilica.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicReturn1.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicReturn2.o)			\
	${LIB}(${BFMSRC}/Silt/ModuleSilt.o)        

#	${LIB}(${BFMSRC}/CO2/CO2System.o)			        \
#	${LIB}(${BFMSRC}/CO2/ModuleBenCO2Transport.o)		        \
#	${LIB}(${BFMSRC}/CO2/ModuleBenAlkalinity.o)		        \
#	${LIB}(${BFMSRC}/CO2/ModulePelCO2.o)                            \


BFM_OBJ = \
	${LIB}(${BFMGOTMSRC}/GetDelta.o)				\
	${LIB}(${BFMGOTMSRC}/D3toD1.o)					\
	${LIB}(${BFMGOTMSRC}/D2toD1.o)					\
	${LIB}(${BFMSRC}/General/AllocateMem.o)				\
	${LIB}(${BFMSRC}/General/set_var_info_bfm.o)			\
	${LIB}(${BFMSRC}/General/check_if_in_output.o)			\
	${LIB}(${BFMSRC}/General/Ecology.o)				\
	${LIB}(${BFMSRC}/General/InitBoxParams.o)			\
	${LIB}(${BFMSRC}/General/Initialize.o)				\
	${LIB}(${BFMSRC}/General/InitTransportStateTypes.o)		\
	${LIB}(${BFMSRC}/General/eTq.o)					\
	${LIB}(${BFMSRC}/General/CheckMassConservation.o)		\
	${LIB}(${BFMSRC}/General/CalcRiverConcentration.o)		\
	${LIB}(${BFMSRC}/Light/LightAdaptation.o)			\
	${LIB}(${BFMSRC}/Light/PhotoAvailableRadiation.o)		\
	${LIB}(${BFMSRC}/Oxygen/WindOxReaeration_3.o)			\
	${LIB}(${BFMSRC}/Oxygen/CalcOxygenSaturation_3.o)		\
	${LIB}(${BFMSRC}/Oxygen/CalcSchmidtNumberOx.o)			\
	${LIB}(${BFMSRC}/PelB/CalcChlorophylla.o)			\
	${LIB}(${BFMSRC}/PelB/CalcVerticalExtinction.o)			\
	${LIB}(${BFMSRC}/PelB/MicroZoo.o)				\
	${LIB}(${BFMSRC}/PelB/MesoZoo.o)				\
	${LIB}(${BFMSRC}/PelB/MicroZoo.o)				\
	${LIB}(${BFMSRC}/PelB/PelBac.o)					\
	${LIB}(${BFMSRC}/PelB/PelChem.o)				\
	${LIB}(${BFMSRC}/PelB/PelGlobal.o)				\
	${LIB}(${BFMSRC}/PelB/PelagicSystem.o)				\
	${LIB}(${BFMSRC}/PelB/Phyto.o)					\
	${LIB}(${BFMSRC}/PelBen/BentoPelCoup.o)				\
	${LIB}(${BFMSRC}/PelBen/PelForcingForBen.o)			\
	${LIB}(${BFMSRC}/PelBen/Sedimentation.o)			\
	${LIB}(${BFMSRC}/PelBen/Settling.o)				\
	${LIB}(${BFMSRC}/PelBen/ControlBennutBuffers.o)			\
	${LIB}(${BFMSRC}/PelBen/RecalcPenetrationDepth.o)		\
	${LIB}(${BFMSRC}/Ben/BenBac.o)					\
	${LIB}(${BFMSRC}/Ben/BenOrganism.o)				\
	${LIB}(${BFMSRC}/Ben/BenthicSystem.o)				\
	${LIB}(${BFMSRC}/Ben/BenGlobal.o)				\
	${LIB}(${BFMSRC}/Ben/Bioturbation.o)				\
	${LIB}(${BFMSRC}/Ben/FilterFeeder.o)				\
	${LIB}(${BFMSRC}/Ben/CorrectConcNearBed.o)			\
	${LIB}(${BFMSRC}/Bennut/BenAmmonium.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenAnoxic.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenDenitriDepth.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenNitrate.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenNitrogenShifting.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenOxygen.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenPhosphate.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenQ1Transport.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenSilica.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenthicNutrient2.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenthicNutrient3.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenthicReturn1.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenthicReturn2.o) 			\
	${LIB}(${BFMSRC}/Bennut/InitBenthicNutrient.o) 			\
	${LIB}(${BFMSRC}/Bennut/InitBenthicNutrient3.o) 		\
	${LIB}(${BFMSRC}/Bennut/BenProfiles.o) 				\
	${LIB}(${BFMSRC}/Bennut/CalculateFromSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/CalculateSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/CalculateShift.o) 			\
	${LIB}(${BFMSRC}/Bennut/LimitShift.o)	 			\
	${LIB}(${BFMSRC}/Bennut/LimitChange.o)	 			\
	${LIB}(${BFMSRC}/Bennut/FixProportionCoeff.o)	 		\
	${LIB}(${BFMSRC}/Bennut/CalculateTau.o) 			\
	${LIB}(${BFMSRC}/Bennut/CompleteSet.o) 				\
	${LIB}(${BFMSRC}/Bennut/InitializeSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/DefineSet.o) 				\
	${LIB}(${BFMSRC}/Bennut/GetInfoFromSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/PrintSet.o) 				\
	${LIB}(${BFMSRC}/Bennut/bess_exp.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessi0.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessi1.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessk0.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessk1.o) 				\
	${LIB}(${BFMSRC}/Bennut/calculate_equation.o) 			\
	${LIB}(${BFMSRC}/Bennut/calculatelayer.o) 			\
	${LIB}(${BFMSRC}/Bennut/funcalc.o) 				\
	${LIB}(${BFMSRC}/Bennut/input_para.o) 				\
	${LIB}(${BFMSRC}/Bennut/kfind.o) 				\
	${LIB}(${BFMSRC}/Bennut/lubksb.o) 				\
	${LIB}(${BFMSRC}/Bennut/ludcmp.o) 				\
	${LIB}(${BFMSRC}/Bennut/manage_coeff.o) 			\
	${LIB}(${BFMSRC}/Bennut/noutput.o) 				\
	${LIB}(${BFMSRC}/Bennut/qgaus_exp.o) 				\
	${LIB}(${BFMSRC}/Bennut/re_store.o) 				\
	${LIB}(${BFMSRC}/Bennut/set_max_sing.o) 			\
	${LIB}(${BFMSRC}/Bennut/svbksb.o) 				\
	${LIB}(${BFMSRC}/Bennut/svdcmp.o) 				\
	${LIB}(${BFMSRC}/Bennut/transfer.o) 				\
 	${LIB}(${BFMSRC}/Silt/Silt.o)        
#	${LIB}(${BFMSRC}/CO2/Alkalinity.o) 				\
#	${LIB}(${BFMSRC}/CO2/PelCO2.o) 					\
#	${LIB}(${BFMSRC}/CO2/CalcSchmidtNumberCO2.o) 			\
#	${LIB}(${BFMSRC}/CO2/SurfaceCO2Processes.o) 			\
#	${LIB}(${BFMSRC}/CO2/CalcCO2SatInField.o) 			\
#	${LIB}(${BFMSRC}/CO2/BenCO2Transport.o) 			\
#	${LIB}(${BFMSRC}/CO2/BenCO2Profiles.o)				\
#	${LIB}(${BFMSRC}/CO2/BenpH.o) 					\
#	${LIB}(${BFMSRC}/CO2/BenAlkalinity.o) 			        \


all: ${BFM_MOD} ${OBJ} ${BFM_OBJ}
	$(MOVE_MODULES_COMMAND)

$(BFM_MOD) : $(BFMSRC)/General/ModuleMem.F90
$(BFMSRC)/General/ModuleMem.o : $(BFMSRC)/General/ModuleMem.F90

${BFMSRC}/General/ModuleMem.F90: $(BFMSRC)/General/GlobalDefsBFM.model $(BFMSRC)/General/FluxFunctions.F90
	${BFMSRC}/scripts/GenerateGlobalBFMF90Code  $(DEFINES) -read ${BFMSRC}/General/GlobalDefsBFM.model \
		-from ${BFMSRC}/proto -to ${BFMSRC}/General -actions statemem allocmem netcdfmem \
		-to ${BFMSRC}/include -actions headermem  
${BFMSRC}/General/AllocateMem.F90: $(BFMSRC)/General/ModuleMem.F90
${BFMSRC}/General/set_var_info_bfm.F90: $(BFMSRC)/General/ModuleMem.F90

#${BFMGOTMSRC}/.depend:
#	makedepf90 -I ${BFMSUBDIRS} *.?90 > $@

realclean: clean
	$(RM) *.o *~
	$(RM) ${BFMSRC}/*/*.o ${BFMGOTMSRC}/*.o ${BFMSHARE}/*.o

distclean: realclean


#-----------------------------------------------------------------------
# Copyright (C) 2006 - the GOTM-team and the BFM-team
#-----------------------------------------------------------------------
