RM = rm -f

INC := -I../gclib
CC := g++
LDFLAGS =
BASEFLAGS := -std=c++11 $(INC) -Wall -Wextra -D_REENTRANT -fno-exceptions -fno-rtti
ifeq ($(findstring release,$(MAKECMDGOALS)),)
  CFLAGS := -g -DDEBUG $(BASEFLAGS)
  LDFLAGS += -g
else
  CFLAGS := -O2 -DNDEBUG $(BASEFLAGS)
endif

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER  := g++
LIBS := 
OBJS := $(INC)/GBase.o $(INC)/gdna.o $(INC)/codons.o \
 $(INC)/GFaSeqGet.o $(INC)/gff.o interval_tree.o $(INC)/GArgs.o

.PHONY : all
all:    mapped_nc

mapped_nc: mapped_nc.o $(OBJS)
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

# target for removing all object files

.PHONY : clean
clean:: 
	@${RM} $(OBJS) mapped_nc.o mapped_nc
	@${RM} core.*


