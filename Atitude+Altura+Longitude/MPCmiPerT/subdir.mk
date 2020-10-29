################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
main.cpp

OBJS += \
main.o



# Each subdirectory must supply rules for building sources it contributes
%.o: %.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@g++ -std=c++11 -I ../include\
	    -I ./qpOASES/include\
	    -include ./include/qpOASES/include/qpOASES.hpp\
	    -include ./include/qpOASES/include/qpOASES/SubjectTo.hpp\
	    -include ./include/qpOASES/include/qpOASES/Bounds.hpp\
	    -include ./include/qpOASES/include/qpOASES/Constants.hpp\
	    -include ./include/qpOASES/include/qpOASES/ConstraintProduct.hpp\
	    -include ./include/qpOASES/include/qpOASES/Types.hpp\
	    -include ./include/qpOASES/include/qpOASES/UnitTesting.hpp\
	    -include ./include/qpOASES/include/qpOASES/Utils.hpp\
	    -include ./include/qpOASES/include/qpOASES/QProblem.hpp\
	    -include ./include/qpOASES/include/qpOASES/LapackBlasReplacement.hpp\
	    -include ./include/qpOASES/include/qpOASES/Constraints.hpp\
	    -include ./include/qpOASES/include/qpOASES/Flipper.hpp\
	    -include ./include/qpOASES/include/qpOASES/Indexlist.hpp\
	    -include ./include/qpOASES/include/qpOASES/Matrices.hpp\
	    -include ./include/qpOASES/include/qpOASES/MessageHandling.hpp\
	    -include ./include/qpOASES/include/qpOASES/Options.hpp\
	    -include ./include/qpOASES/include/qpOASES/QProblemB.hpp\
	    -include ./include/qpOASES/include/qpOASES/SparseSolver.hpp\
	    -include ./include/qpOASES/include/qpOASES/SQProblem.hpp\
	    -include ./include/qpOASES/include/qpOASES/SQProblemSchur.hpp\
	    -include ./include/qpOASES/src/SubjectTo.cpp\
	    -include ./include/qpOASES/src/Bounds.cpp\
	    -include ./include/qpOASES/src/QProblem.cpp\
	    -include ./include/qpOASES/src/BLASReplacement.cpp\
	    -include ./include/qpOASES/src/Constraints.cpp\
	    -include ./include/qpOASES/src/Flipper.cpp\
	    -include ./include/qpOASES/src/Indexlist.cpp\
	    -include ./include/qpOASES/src/LAPACKReplacement.cpp\
	    -include ./include/qpOASES/src/Matrices.cpp\
	    -include ./include/qpOASES/src/MessageHandling.cpp\
	    -include ./include/qpOASES/src/Options.cpp\
	    -include ./include/qpOASES/src/OQPinterface.cpp\
	    -include ./include/qpOASES/src/QProblemB.cpp\
	    -include ./include/qpOASES/src/SolutionAnalysis.cpp\
	    -include ./include/qpOASES/src/SparseSolver.cpp\
	    -include ./include/qpOASES/src/SQProblem.cpp\
	    -include ./include/qpOASES/src/SQProblemSchur.cpp\
	    -include ./include/qpOASES/src/Utils.cpp\
	    -include ./include/matrizesBANCADA.cpp\
	    -include ./include/Controlador/CONTROL.cpp\
	    -include ./include/Controlador/include/matrizesMPC.cpp\
	    -O2 -fwhole-program -c\
	    -o "src/$@" "$<"
	@echo 'Finished building: $<'
	@echo ''


