################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/poisson/poisson.cpp \
../src/poisson/poisson_static.cpp 

OBJS += \
./src/poisson/poisson.o \
./src/poisson/poisson_static.o 

CPP_DEPS += \
./src/poisson/poisson.d \
./src/poisson/poisson_static.d 


# Each subdirectory must supply rules for building sources it contributes
src/poisson/%.o: ../src/poisson/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -g -I/home/dmarce1/include -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


