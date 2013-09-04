################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/single_star/lane_emden.cpp \
../src/single_star/single_star.cpp \
../src/single_star/single_star_static.cpp \
../src/single_star/wd.cpp 

OBJS += \
./src/single_star/lane_emden.o \
./src/single_star/single_star.o \
./src/single_star/single_star_static.o \
./src/single_star/wd.o 

CPP_DEPS += \
./src/single_star/lane_emden.d \
./src/single_star/single_star.d \
./src/single_star/single_star_static.d \
./src/single_star/wd.d 


# Each subdirectory must supply rules for building sources it contributes
src/single_star/%.o: ../src/single_star/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -ipo -inline-level=2 -I/home/dmarce1/include -I/usr/include/x86_64-linux-gnu/c++/4.7/ -DNDEBUG -xHost -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


