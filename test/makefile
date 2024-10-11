##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for core SMS++)                                        #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# module name
# NAME can be both main and test, the corresponding .cpp is compiled
# but the final executable is always called "test"
NAME = test

# basic directory
DIR = .

# here is where the specific MCFClass to be tested is decided
# debug switches
SW = -g -std=c++17 -ferror-limit=1 -Wno-enum-compare

# production switches
#SW = -O3 -DNDEBUG -std=c++17 -Wno-enum-compare

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

all default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# BundleSolver
BNDSLVSDR = ..
include $(BNDSLVSDR)/makefile-c

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(BNDSLVOBJ)

# libraries
MLIB =  $(BNDSLVLIB)

$(NAME): $(MOBJ) $(DIR)/$(NAME).o
	$(CC) -o $(NAME) $(DIR)/$(NAME).o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

# include directives
MINC =  $(BNDSLVINC)

# includes
MH =    $(BNDSLVH)

# compile command

$(DIR)/$(NAME).o: $(DIR)/$(NAME).cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) $(SW)

############################ End of makefile #################################
