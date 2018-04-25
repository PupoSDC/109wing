import os
import math
import copy
import sys

sys.path.append('/usr/lib/freecad/lib/');
import FreeCAD
import Part

# Reads a 2D array aerofoil point index, and converts it to a 3D point array
def readAerofoil( file, position, scale ):
    response = [];
    with open(file,'r') as file:
        for line in file:
            response.append( map( float, line.split() ) )
            response[-1][0] = response[-1][0] * scale + position[0];
            response[-1][1] = response[-1][1] * scale + position[1];
            response[-1].append( position[2] );
            response[-1] = App.Vector(
            	response[-1][0],
            	response[-1][1],
            	response[-1][2]
            );
    return response;

NACA_2411 = os.path.join(os.path.dirname(__file__), 'sources/NACA/2411')
NACA_2414 = os.path.join(os.path.dirname(__file__), 'sources/NACA/2414')

doc        = App.activeDocument();
root       = doc.addObject("Part::Polygon","Naca_2414");
root.Nodes = readAerofoil(NACA_2414, [ -0.958, 0.0,   0.0   ], 2.143 );
tip        = doc.addObject("Part::Polygon","Naca_2411");
tip.Nodes  = readAerofoil(NACA_2411, [ -0.473, 0.109, 4.215 ], 1.050 );

doc.addObject('Part::Loft','main-wing-section')
doc.ActiveObject.Sections= [ root, tip ]
doc.ActiveObject.Solid  = True;
doc.ActiveObject.Ruled  = False;
doc.ActiveObject.Closed = False;

doc.recompute();

