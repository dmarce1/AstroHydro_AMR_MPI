################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/poisson_test/poisson_test.cpp \
../src/poisson_test/poisson_test_static.cpp 

OBJS += \
./src/poisson_test/poisson_test.o \
./src/poisson_test/poisson_test_static.o 

CPP_DEPS += \
./src/poisson_test/poisson_test.d \
./src/poisson_test/poisson_test_static.d 


# Each subdirectory must supply rules for building sources it contributes
src/poisson_test/%.o: ../src/poisson_test/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -opt-prefetch=3 -ipo -inline-level=2 -I/home/dmarce1/include -DNDEBUG -ansi-alias -fargument-noalias -fno-alias -fp-speculation=fast -fp-model fast=2 -xHost -Wuninitialized -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


