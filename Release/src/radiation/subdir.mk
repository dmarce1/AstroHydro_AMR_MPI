################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/radiation/radiation.cpp 

OBJS += \
./src/radiation/radiation.o 

CPP_DEPS += \
./src/radiation/radiation.d 


# Each subdirectory must supply rules for building sources it contributes
src/radiation/%.o: ../src/radiation/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpicc -O3 -xHOST -ipo -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


