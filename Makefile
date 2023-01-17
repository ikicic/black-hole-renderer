# CPP= clang++-3.5
CPP ?= g++
WARNINGFLAGS=-Wall -Wextra
CPPFLAGS=-D_GLIBCXX_USE_CXX11_ABI=0 -O3 -march=native -std=c++1z -fdiagnostics-color -g $(WARNINGFLAGS) #-ffast-math -fno-finite-math-only
BIN=otr
DEPDIR=.d
DEPFLAGS=-MT $@ -MMD -MP -MF
LIBS=-pthread
SOURCES= $(shell find src/bhr tests -type f -name "*.cpp")
OBJ= $(patsubst %,build/%,$(SOURCES:.cpp=.o))

RM= rm -f

$(shell mkdir -p $(DEPDIR) >/dev/null)

.PHONY: all clean run


all: $(BIN)

clean:
	${RM} $(OBJ) $(BIN) $(DEPDIR)/*.d

run:
	./$(BIN)

$(BIN): $(OBJ)
	$(CPP) $(OBJ) -o $(BIN) $(LIBS)

build/src/%.o: src/%.cpp $(DEPDIR)/src/%.d
	@mkdir -pv $(dir $(DEPDIR)/src/$*) $(dir build/src/$*)
	$(CPP) -I src $(DEPFLAGS) $(DEPDIR)/src/$*.Td -c $< -o $@ $(CPPFLAGS)
	@mv -f $(DEPDIR)/src/$*.Td $(DEPDIR)/src/$*.d

build/tests/%.o: tests/%.cpp $(DEPDIR)/tests/%.d
	@mkdir -pv $(dir $(DEPDIR)/tests/$*) $(dir build/tests/$*)
	$(CPP) -I src $(DEPFLAGS) $(DEPDIR)/tests/$*.Td -c $< -o $@ $(CPPFLAGS)
	@mv -f $(DEPDIR)/tests/$*.Td $(DEPDIR)/tests/$*.d

$(DEPDIR)/src/%.d: ;
.PRECIOUS: $(DEPDIR)/src/%.d

$(DEPDIR)/tests/%.d: ;
.PRECIOUS: $(DEPDIR)/tests/%.d

-include $(patsubst %,$(DEPDIR)/src/%.d,$(basename $(SOURCES:src/%=%)))
-include $(patsubst %,$(DEPDIR)/tests/%.d,$(basename $(SOURCES:tests/%=%)))
