################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/euler_state/state.cpp 

OBJS += \
./src/euler_state/state.o 

CPP_DEPS += \
./src/euler_state/state.d 


# Each subdirectory must supply rules for building sources it contributes
src/euler_state/%.o: ../src/euler_state/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -g -I/home/dmarce1/include -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


