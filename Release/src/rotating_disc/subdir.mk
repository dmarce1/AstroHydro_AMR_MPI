################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/rotating_disc/rotating_disc.cpp \
../src/rotating_disc/state.cpp 

OBJS += \
./src/rotating_disc/rotating_disc.o \
./src/rotating_disc/state.o 

CPP_DEPS += \
./src/rotating_disc/rotating_disc.d \
./src/rotating_disc/state.d 


# Each subdirectory must supply rules for building sources it contributes
src/rotating_disc/%.o: ../src/rotating_disc/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -opt-prefetch=3 -ipo -inline-level=2 -I/home/dmarce1/include -DNDEBUG -ansi-alias -fargument-noalias -fno-alias -fp-model strict -xHost -Wuninitialized -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


