#!/usr/bin/env python
import os
import math
import copy

################################################################################
## BOCK MESH OBJECT DEFINITION #################################################
################################################################################

class BlockMesh:
    pointlist = [];
    hexblocks = [];
    boundlist = [];

    # Adds a point to a point list only if it does not exist previously.
    def addPoint(self, point ):
        point = [
            float( '%.6f' % point[0] ),
            float( '%.6f' % point[1] ),
            float( '%.6f' % point[2] )
        ];
        if point not in self.pointlist: self.pointlist.append(point)
        return self.pointlist.index(point);

    # Creates an HExablock between 8 points
    def addBlock(self, points, n_points, grade ):
        self.hexblocks.append([
            self.addPoint( points[0] ), self.addPoint( points[1] ),
            self.addPoint( points[2] ), self.addPoint( points[3] ),
            self.addPoint( points[4] ), self.addPoint( points[5] ),
            self.addPoint( points[6] ), self.addPoint( points[7] ),
            n_points, grade
        ]);

    # Adds a Face to the boundary list
    def addFace(self, name, points, type ):
        p_list = map( lambda x: self.addPoint(x), points );
        for face in self.boundlist:
            if face[0] == name:
                face[1].append( p_list );
                return;
        self.boundlist.append([name,[p_list],type]);

    # Prints out the blockMeshDict
    def write(self):
        blockMeshDict = open("system/blockMeshDict","w")

        # Write the heading
        blockMeshDict.write(
            "FoamFile\n{\n    version     2.0;\n" +
            "format      ascii;\n" +
            "class       dictionary;\n" +
            "object      blockMeshDict;\n}\n" +
            "convertToMeters 1;\n"
        );

        # Write down the vertices
        blockMeshDict.write( "\nvertices\n(\n")
        for point in self.pointlist:
            blockMeshDict.write(
                "    (" +
                repr( point[0] ) + " " +
                repr( point[1] ) + " " +
                repr( point[2] ) + ") // " +
                repr( self.pointlist.index(point) ) + "\n"
            );

        # Write down the blocks
        blockMeshDict.write( ");\n\nblocks\n(\n")
        for hex in self.hexblocks:
            blockMeshDict.write(
                "  hex ( " +
                repr(hex[0]) + " " +
                repr(hex[1]) + " " +
                repr(hex[2]) + " " +
                repr(hex[3]) + " " +
                repr(hex[4]) + " " +
                repr(hex[5]) + " " +
                repr(hex[6]) + " " +
                repr(hex[7]) + " ) (" +
                repr(hex[8][0]) + " " +
                repr(hex[8][1]) + " " +
                repr(hex[8][2]) + ") simpleGrading (" +
                repr(hex[9][0]) +" "+
                repr(hex[9][1]) +" "+
                repr(hex[9][2]) + ")\n"
            );

        # Write down the edges
        blockMeshDict.write( ");\n\nedges\n(\n")

        # Write down the boundaries
        blockMeshDict.write( ");\n\nboundary\n(\n")
        for boundary in self.boundlist:
            faces = "";
            for face in boundary[1]:
                faces += "( "
                for point in face:
                    faces += " " + repr(point) + " ";
                faces += ")";
            blockMeshDict.write(
                "    "
                + boundary[0] + " { type "
                + boundary[2] + "; faces ("
                + faces       + ");}\n"
            );

        #write down merged patches
        blockMeshDict.write( ");\n\nmergePatchPairs\n(\n);\n")
        blockMeshDict.close()
        return;

################################################################################
## FUNCTION DEFINITIONS ########################################################
################################################################################

# Reads a 2D array aerofoil point index, and converts it to a 3D point array
def readAerofoil( file, position, scale ):
    response = [];
    with open(file,'r') as file:
        for line in file:
            response.append( map( float, line.split() ) )
            response[-1][0] = response[-1][0] * scale + position[0];
            response[-1][1] = response[-1][1] * scale + position[1];
            response[-1].append(  position[2] );
    response.append([
        ( response[0][0] + response[-1][0] ) / 2.0,
        ( response[0][1] + response[-1][1] ) / 2.0,
        ( response[0][2] + response[-1][2] ) / 2.0
    ]);
    return response;

# Creates a wing section
def makeWingSection( block_mesh, section_1, section_2, radius, points, grade ):
    for i in range( 0, len(section_1) ):
        a_i = ( float(i+0) / (len(section_1) ) ) * 2.0 * math.pi;
        a_j = ( float(i+1) / (len(section_1) ) ) * 2.0 * math.pi;
        ps = [
            section_1[ i-1 ],
            section_1[ i+0 ],
            section_2[ i+0 ],
            section_2[ i-1 ],
            [ radius*math.cos(a_i), radius*math.sin(a_i), section_1[i-1][2] ],
            [ radius*math.cos(a_j), radius*math.sin(a_j), section_1[i+0][2] ],
            [ radius*math.cos(a_j), radius*math.sin(a_j), section_2[i+0][2] ],
            [ radius*math.cos(a_i), radius*math.sin(a_i), section_2[i-1][2] ]
        ];

        block_mesh.addBlock( ps, points, grade );
        block_mesh.addFace( "wing",       ps[0:4], "wall" );
        block_mesh.addFace( "freestream", ps[4:8], "patch" );
        if section_1[i][2] == 0.0:
            block_mesh.addFace("left", [ps[0],ps[1],ps[4],ps[5]], "symmetry");

#Creates a winglet section.
def makeWinglet( block_mesh, section, points, radius, n_points, grade ):
    reference = section[ len(section)/2+1 ];
    for i in range(0, len(section)/2 ):
        for j in range(0, points     ):
            # calculate  auxiliary values
            d_x     =    ( section[-1-i][0]   - section[i][0] );
            r_w     = abs( section[-1-i][1]   - section[i][1] ) / 2.0;
            c_w     =    ( section[-1-i][1]   + section[i][1] ) / 2.0;

            d_x_p   =    ( section[-2-i][0] - section[i+1][0] );
            r_w_p   = abs( section[-2-i][1] - section[i+1][1] ) / 2.0;
            c_w_p   =    ( section[-2-i][1] + section[i+1][1] ) / 2.0;

            phi     = math.pi * float(j)   / float( points );    # sideways angle ( 0 to 180 )
            theta   = math.pi * float(i)   / (len(section)-1);   # aerofoil angle ( 0 tp 90  )
            phi_p   = math.pi * float(j+1) / float( points );    # sideways angle ( 0 to 180 )
            theta_p = math.pi * float(i+1) / (len(section)-1);   # aerofoil angle ( 0 tp 90  )

            ps = [
                [ # first point on the aerofoil
                    section[i][0] + d_x * float(j) / float( points ),
                    c_w           + r_w * math.cos( phi ),
                    section[i][2] + r_w * math.sin( phi )
                ],
                [ # second point in the aerofoil
                    section[i+1][0] + d_x_p * float(j) / float( points ),
                    c_w_p           + r_w_p * math.cos( phi ),
                    section[i+1][2] + r_w_p * math.sin( phi )
                ],
                [ # extrude the point in the aerofoil
                    section[i+1][0] + d_x_p * float(j+1) / float( points ),
                    c_w_p           + r_w_p * math.cos( phi_p ),
                    section[i+1][2] + r_w_p * math.sin( phi_p )
                ],
                [ # next point in extruded aerofoil
                    section[i][0] + d_x * float(j+1) / float( points ),
                    c_w           + r_w * math.cos( phi_p ),
                    section[i][2] + r_w * math.sin( phi_p )
                ],
                [ # first point in away mesh
                    section[0][0] - (radius+section[0][0] ) * math.sin( theta ),
                                  + radius                  * math.cos( theta ) * math.cos( phi ),
                    section[i][2] + radius                  * math.cos( theta ) * math.sin( phi )
                ],
                [ # second point in away mesh
                    section[0][0]   - (radius+section[0][0] ) * math.sin( theta_p ),
                                    + radius                  * math.cos( theta_p ) * math.cos( phi ),
                    section[i+1][2] + radius                  * math.cos( theta_p ) * math.sin( phi )
                ],
                [ # third point in away mesh
                    section[0][0]   - (radius+section[0][0] ) * math.sin( theta_p ),
                                    + radius                  * math.cos( theta_p ) * math.cos( phi_p ),
                    section[i+1][2] + radius                  * math.cos( theta_p ) * math.sin( phi_p )
                ],
                [ # forth point in away mesh
                    section[0][0] - (radius+section[0][0] ) * math.sin( theta ),
                                  + radius                  * math.cos( theta ) * math.cos( phi_p ),
                    section[i][2] + radius                  * math.cos( theta ) * math.sin( phi_p )
                ]
            ];
            block_mesh.addBlock( ps, n_points, grade );

            if i == 0:
                a1 = [
                    ( section[0][0] + section[-1][0] ) / 2.0,
                    ( section[0][1] + section[-1][1] ) / 2.0,
                    ( section[0][2] + section[-1][2] ) / 2.0
                ];

                a2 = [
                    radius,
                    ( section[0][1] + section[-1][1] ) / 2.0,
                    ( section[0][2] + section[-1][2] ) / 2.0
                ];

                #must add C_w
                block_mesh.addBlock(
                    [
                        a1, # centre point
                        [ ps[0][0], ps[0][1], ps[0][2] ], # point 0
                        [ ps[3][0], ps[3][1], ps[3][2] ], # point 3
                        a1, # centre point
                        a2, # centre point away
                        [ radius,        ps[0][1], ps[0][2] ], # point 0 extrusion to wall
                        [ radius,        ps[3][1], ps[3][2] ], # point 3 extrusion to wall
                        a2, # centre point away
                    ],
                    [ n_points[1], n_points[1], n_points[2] ],
                    [ grade[1],    grade[1],    grade[2]    ]
                );
                block_mesh.addBlock(
                    [
                        [ ps[4][0], ps[4][1], ps[4][2] ], # point 4
                        [ radius,   ps[4][1], ps[4][2] ], # point 4 extrusion to wall
                        [ radius,   ps[7][1], ps[7][2] ], # point 7 extrusion to wall
                        [ ps[7][0], ps[7][1], ps[7][2] ], # point 7
                        [ ps[0][0], ps[0][1], ps[0][2] ], # copy point 0
                        [ radius,   ps[0][1], ps[0][2] ], # point 0 extrusion to wall
                        [ radius,   ps[3][1], ps[3][2] ], # point 3 extrusion to wall
                        [ ps[3][0], ps[3][1], ps[3][2] ], # copy point 3
                    ],
                    [ n_points[2], n_points[1], n_points[2] ],
                    [ grade[2],    grade[1],    1.0 / grade[2]    ]
                );
                block_mesh.addFace(
                    "freestream", [
                        [ ps[4][0], ps[4][1], ps[4][2] ],
                        [ radius,   ps[4][1], ps[4][2] ],
                        [ radius,   ps[7][1], ps[7][2] ],
                        [ ps[7][0], ps[7][1], ps[7][2] ]
                    ], "patch"
                );
                block_mesh.addFace(
                    "freestream", [
                        [ radius,   ps[4][1], ps[4][2] ],
                        [ radius,   ps[7][1], ps[7][2] ],
                        [ radius,   ps[0][1], ps[0][2] ],
                        [ radius,   ps[3][1], ps[3][2] ]
                    ], "patch"
                );
                block_mesh.addFace(
                    "wing", [
                        a1, # centre point
                        [ ps[0][0], ps[0][1], ps[0][2] ], # point 0
                        [ ps[3][0], ps[3][1], ps[3][2] ], # point 3
                        a1  # centre point
                    ], "wall"
                );
                block_mesh.addFace(
                    "freestream", [
                        a2, # centre point away
                        [ radius,        ps[0][1], ps[0][2] ], # point 0 extrusion to wall
                        [ radius,        ps[3][1], ps[3][2] ], # point 3 extrusion to wall
                        a2  # centre point away
                    ], "patch"
                );


            block_mesh.addFace(
                "wing", [ ps[0], ps[1], ps[2], ps[3] ], "wall"
            );
            block_mesh.addFace(
                "freestream", [ ps[4], ps[5], ps[6], ps[7] ], "patch"
            );

################################################################################
## Script ######################################################################
################################################################################

block_mesh = BlockMesh();

root = readAerofoil('constant/aerofoils/NACA_2414', [ 0.958, 0.0,   0.0   ], 2.143 );
tip  = readAerofoil('constant/aerofoils/NACA_2411', [ 0.470, 0.109, 4.215 ], 1.050 );

makeWingSection(
    block_mesh,
    root,
    tip,
    4.0,
    [1,1,20],    # 30,20
    [1,1,1000]  # 1000
);

#makeWinglet(
#    block_mesh,
#    tip,
#    15,
#    15.0,
#    [1,1,20],    # 1,20
#    [1,1,1000]  # 1000
#);

block_mesh.write();

