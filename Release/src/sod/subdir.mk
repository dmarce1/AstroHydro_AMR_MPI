################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/sod/exact_sod.cpp \
../src/sod/grid_sod.cpp 

OBJS += \
./src/sod/exact_sod.o \
./src/sod/grid_sod.o 

CPP_DEPS += \
./src/sod/exact_sod.d \
./src/sod/grid_sod.d 


# Each subdirectory must supply rules for building sources it contributes
src/sod/%.o: ../src/sod/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -opt-prefetch=3 -ipo -inline-level=2 -I/home/dmarce1/include -DNDEBUG -fp-speculation=off -fp-model strict -xHost -Wuninitialized -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


