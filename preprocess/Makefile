
CXX=g++
CXXFLAGS += -Wall -Wno-unused-result -pedantic -O3 -mtune=native -std=c++11 -g -ggdb
CXXLIBS =
INCDIRS = -I./src
LIBDIRS =

OBJDIR = obj
SRCDIR = src

DEPS = $(shell find $(SRCDIR) -name '*.h')
SRCS = $(shell find $(SRCDIR) -name '*.cpp')
OBJS = $(patsubst $(SRCDIR)%.cpp, $(OBJDIR)%.o, $(SRCS))

all: $(OBJDIR) get_features

$(OBJDIR):
	mkdir -p $(OBJDIR)

get_features: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o get_features $(OBJS) $(CXXLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm $(OBJS) get_features
	rmdir $(OBJDIR)
