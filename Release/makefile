################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/stellar_state/subdir.mk
-include src/sod/subdir.mk
-include src/single_star/subdir.mk
-include src/rotating_disc/subdir.mk
-include src/radiation_test/subdir.mk
-include src/radiation/subdir.mk
-include src/poisson_test/subdir.mk
-include src/poisson/subdir.mk
-include src/oct_node/subdir.mk
-include src/multigrid/subdir.mk
-include src/hydro_rad_grid/subdir.mk
-include src/hydro_grid/subdir.mk
-include src/hydro_grav_grid/subdir.mk
-include src/hydro_FMM_grid/subdir.mk
-include src/euler_state/subdir.mk
-include src/blast_wave/subdir.mk
-include src/binary_star/subdir.mk
-include src/FMM/subdir.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: AstroHydro_AMR_MPI

# Tool invocations
AstroHydro_AMR_MPI: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Intel Intel(R) 64 C++ Linker'
	mpic++ -openmp -L"/home/dmarce1/AstroHydro_AMR_MPI/iso7_src" -L/home/dmarce1/lib -O3 -o "AstroHydro_AMR_MPI" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) AstroHydro_AMR_MPI
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
