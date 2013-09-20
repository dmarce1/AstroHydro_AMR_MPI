################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/hydro_rad_grid/hydro_rad_grid.cpp 

OBJS += \
./src/hydro_rad_grid/hydro_rad_grid.o 

CPP_DEPS += \
./src/hydro_rad_grid/hydro_rad_grid.d 


# Each subdirectory must supply rules for building sources it contributes
src/hydro_rad_grid/%.o: ../src/hydro_rad_grid/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpicc -O3 -xHOST -ipo -openmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


