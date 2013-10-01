################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/hydro_FMM_grid/hydro_FMM_grid.cpp \
../src/hydro_FMM_grid/moment.cpp \
../src/hydro_FMM_grid/taylor.cpp 

OBJS += \
./src/hydro_FMM_grid/hydro_FMM_grid.o \
./src/hydro_FMM_grid/moment.o \
./src/hydro_FMM_grid/taylor.o 

CPP_DEPS += \
./src/hydro_FMM_grid/hydro_FMM_grid.d \
./src/hydro_FMM_grid/moment.d \
./src/hydro_FMM_grid/taylor.d 


# Each subdirectory must supply rules for building sources it contributes
src/hydro_FMM_grid/%.o: ../src/hydro_FMM_grid/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpicc -DNDEBUG -fp-model fast=2 -O3 -xHOST -ipo -no-prec-div -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


