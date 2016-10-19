################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/auxiliaries/Rand.cpp 

OBJS += \
./src/auxiliaries/Rand.o 

CPP_DEPS += \
./src/auxiliaries/Rand.d 


# Each subdirectory must supply rules for building sources it contributes
src/auxiliaries/%.o: ../src/auxiliaries/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++-4.8 -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


