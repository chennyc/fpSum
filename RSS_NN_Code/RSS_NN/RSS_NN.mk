RSS_SRCS += \
../RSS_NN/benchmarks/Program.cpp \
../RSS_NN/benchmarks/Neural.cpp \
../RSS_NN/benchmarks/benchmark_nn.cpp \
../RSS_NN/utilities/init.cpp \
../RSS_NN/utilities/aes_ni.c \
../RSS_NN/utilities/util.cpp \
../RSS_NN/protocols/open.cpp \
../RSS_NN/protocols/matMult.cpp \
../RSS_NN/protocols/mult.cpp \
../RSS_NN/protocols/msb.cpp \
../RSS_NN/protocols/randBit.cpp \
../RSS_NN/protocols/lt.cpp \
../RSS_NN/protocols/neural_ops.cpp \
../RSS_NN/protocols/svm_ops.cpp \
../RSS_NN/protocols/edaBit.cpp \
../RSS_NN/protocols/conversion.cpp \
../RSS_NN/protocols/bitAdd.cpp  \
../RSS_NN/protocols/bitAddTrunc.cpp  \
../RSS_NN/protocols/bitLT.cpp  \
../RSS_NN/protocols/trunc.cpp  \
../RSS_NN/protocols/b2a.cpp  \
../RSS_NN/protocols/reveal.cpp \
../RSS_NN/protocols/prefixMult.cpp \
../RSS_NN/benchmarks/benchmark_fpSum.cpp \
../RSS_NN/protocols/int2mask.cpp \
../RSS_NN/protocols/fpSum.cpp \
../RSS_NN/protocols/allor.cpp \
../RSS_NN/protocols/b2u.cpp \
../RSS_NN/protocols/bitDec.cpp \
../RSS_NN/protocols/shift.cpp \
../RSS_NN/protocols/b2a3.cpp \
../RSS_NN/protocols/eqz.cpp \
../RSS_NN/protocols/superSum.cpp \
../RSS_NN/protocols/sa2fl.cpp \

RSS_OBJS += \
./RSS_NN/benchmarks/Program.o \
./RSS_NN/benchmarks/Neural.o \
./RSS_NN/benchmarks/benchmark_nn.o \
./RSS_NN/utilities/init.o \
./RSS_NN/utilities/aes_ni.o \
./RSS_NN/utilities/util.o \
./RSS_NN/protocols/open.o \
./RSS_NN/protocols/matMult.o \
./RSS_NN/protocols/mult.o \
./RSS_NN/protocols/msb.o \
./RSS_NN/protocols/randBit.o \
./RSS_NN/protocols/lt.o \
./RSS_NN/protocols/neural_ops.o \
./RSS_NN/protocols/svm_ops.o \
./RSS_NN/protocols/edaBit.o \
./RSS_NN/protocols/conversion.o \
./RSS_NN/protocols/bitAdd.o \
./RSS_NN/protocols/bitAddTrunc.o \
./RSS_NN/protocols/bitLT.o \
./RSS_NN/protocols/trunc.o \
./RSS_NN/protocols/b2a.o \
./RSS_NN/protocols/reveal.o \
./RSS_NN/protocols/prefixMult.o \
./RSS_NN/benchmarks/benchmark_fpSum.o \
./RSS_NN/protocols/int2mask.o \
./RSS_NN/protocols/fpSum.o \
./RSS_NN/protocols/allor.o \
./RSS_NN/protocols/b2u.o \
./RSS_NN/protocols/bitDec.o \
./RSS_NN/protocols/shift.o \
./RSS_NN/protocols/b2a3.o \
./RSS_NN/protocols/eqz.o \
./RSS_NN/protocols/superSum.o \
./RSS_NN/protocols/sa2fl.o \

CPP_DEPS += \
./RSS_NN/benchmarks/Program.d \
./RSS_NN/benchmarks/Neural.d \
./RSS_NN/benchmarks/benchmark_nn.d \
./RSS_NN/utilities/init.d \
./RSS_NN/utilities/aes_ni.d \
./RSS_NN/utilities/util.d \
./RSS_NN/protocols/open.d \
./RSS_NN/protocols/matMult.d \
./RSS_NN/protocols/mult.d \
./RSS_NN/protocols/msb.d \
./RSS_NN/protocols/randBit.d \
./RSS_NN/protocols/lt.d \
./RSS_NN/protocols/neural_ops.d \
./RSS_NN/protocols/svm_ops.d \
./RSS_NN/protocols/edaBit.d \
./RSS_NN/protocols/conversion.d \
./RSS_NN/protocols/bitAdd.d  \
./RSS_NN/protocols/bitAddTrunc.d  \
./RSS_NN/protocols/bitLT.d  \
./RSS_NN/protocols/trunc.d  \
./RSS_NN/protocols/b2a.d  \
./RSS_NN/protocols/prefixMult.d \
./RSS_NN/benchmarks/benchmark_fpSum.d \
./RSS_NN/protocols/int2mask.d \
./RSS_NN/protocols/fpSum.d \
./RSS_NN/protocols/allor.d \
./RSS_NN/protocols/b2u.d \
./RSS_NN/protocols/reveal.d \
./RSS_NN/protocols/bitDec.d \
./RSS_NN/protocols/shift.d \
./RSS_NN/protocols/b2a3.d \
./RSS_NN/protocols/eqz.d \
./RSS_NN/protocols/superSum.d \
./RSS_NN/protocols/sa2fl.d \

CFLAGS =  -maes -msse4.1 -mbmi2 -mavx -mavx2 -Ofast -g

#Supply the rules for building the source
RSS_NN/%.o: ../RSS_NN/%.c
	@echo 'Building file:'
	@echo 'Invoking: GCC C++ Compiler'
	gcc -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '
