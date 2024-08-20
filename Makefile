###############################################################################
#                                                                             #
# Makefile to build libaed-light                                              #
#                                                                             #
#  Developed by :                                                             #
#      AquaticEcoDynamics (AED) Group                                         #
#      School of Agriculture and Environment                                  #
#      The University of Western Australia                                    #
#                                                                             #
#      http://aquatic.science.uwa.edu.au/                                     #
#                                                                             #
#  Copyright 2023 - 2024  -  The University of Western Australia              #
#                                                                             #
#   AED is free software: you can redistribute it and/or modify               #
#   it under the terms of the GNU General Public License as published by      #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   AED is distributed in the hope that it will be useful,                    #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

LIBAEDLGT=aed-lighting
OUTLIB=lib$(LIBAEDLGT)

INCLUDES=-I../libaed-water/${incdir}  -I../libaed-water/${moddir}

include ../libaed-water/make_defs.inc

OBJS=${objdir}/aed_maths.o      \
     ${objdir}/aed_atmosphere.o \
     ${objdir}/aed_slingo.o     \
     ${objdir}/curtin_light.o   \
     ${objdir}/aed_oasim.o      \
     ${objdir}/aed_light.o      \
     ${objdir}/aed_lighting.o

include ../libaed-water/make_rules.inc
