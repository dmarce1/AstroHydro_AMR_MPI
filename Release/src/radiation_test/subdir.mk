################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/radiation_test/radiation_test.cpp \
../src/radiation_test/radiation_test_static.cpp 

OBJS += \
./src/radiation_test/radiation_test.o \
./src/radiation_test/radiation_test_static.o 

CPP_DEPS += \
./src/radiation_test/radiation_test.d \
./src/radiation_test/radiation_test_static.d 


# Each subdirectory must supply rules for building sources it contributes
src/radiation_test/%.o: ../src/radiation_test/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -I/home/dmarce1/include -DNDEBUG -ipo -xHOST -O3 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


