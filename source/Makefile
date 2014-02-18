## This makefile is for gcc-4.1.2
ARNOLDINC := PATH_TO_ARNOLD_INCLUDE_DIR
ARNOLDLIB := PATH_TO_ARNOLD_INCLUDE_LIB

# if your g++ is located somewhere else
GXXLOC    := 
TARGETLOC := 

INCLUDES := -I$(ARNOLDINC)
LDFLAGS  := -L$(ARNOLDLIB) -lai
CXXFLAGS := -w -O2 -fPIC -DLINUX -D_BOOL -ffast-math -DGCC_VERSION=4.1.2
CXX	     := $(GXXLOC)g++

SOURCES  := jf_nested_dielectric.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)jf_nested_dielectric.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
