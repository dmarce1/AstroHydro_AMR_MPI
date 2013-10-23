################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/binary_star/binary_star.cpp \
../src/binary_star/binary_star_static.cpp 

OBJS += \
./src/binary_star/binary_star.o \
./src/binary_star/binary_star_static.o 

CPP_DEPS += \
./src/binary_star/binary_star.d \
./src/binary_star/binary_star_static.d 


# Each subdirectory must supply rules for building sources it contributes
src/binary_star/%.o: ../src/binary_star/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -opt-prefetch=3 -ipo -inline-level=2 -I/home/dmarce1/include -DNDEBUG -ansi-alias -fargument-noalias -fno-alias -fp-speculation=fast -fp-model fast=2 -xHost -Wuninitialized -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


