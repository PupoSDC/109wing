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
        if point not in self.pointlist: self.pointlist.append(point)
        return self.pointlist.index(point);

    # Creates an HExablock between 8 points
    def addBlock(self, points, n_points, grade ):
        self.hexblocks.append([
            self.addPoint( points[0] ), self.addPoint( points[1] ),
            self.addPoint( points[2] ), self.addPoint( points[3] ),
            self.addPoint( points[4] ), self.addPoint( points[5] ),
            self.addPoint( points[6] ), self.addPoint( points[7] ), n_points, grade
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
def readAerofoil(file, position, scale, rotation):
    with open(file,'r') as file:
        array = []
        for line in file:
            array.append( map( float, line.split() ) )
            array[-1].append(0.0);
    response = [];
    for point in array:
        response.append([
            point[0] * scale + position[0],
            point[1] * scale + position[1],
            point[2] + position[2]
        ]);
    return response;

# Creates a wing section
def makeWingSection( block_mesh, section_1, section_2, radius, n_points, grade ):
    for i in range( 0, len(section_1) - 1 ):
        alpha_i = ( 1.0 / 2.0 + float(i+0) / (len(section_1)-1) ) * math.pi;
        alpha_j = ( 1.0 / 2.0 + float(i+1) / (len(section_1)-1) ) * math.pi;
        points = [
            section_1[ i+0 ],
            section_1[ i+1 ],
            section_2[ i+1 ],
            section_2[ i+0 ],
            [
                section_1[0][0] + (radius+section_1[0][0]) * math.cos( alpha_i ),
                radius * math.sin( alpha_i ),
                section_1[i][2]
            ],
            [
                section_1[0][0] + (radius+section_1[0][0]) * math.cos( alpha_j ),
                radius * math.sin( alpha_j ),
                section_1[i][2]
            ],
            [
                section_2[0][0] + (radius+section_2[0][0]) * math.cos( alpha_j ),
                radius * math.sin( alpha_j ),
                section_2[i][2]
            ],
            [
                section_2[0][0] + (radius+section_2[0][0]) * math.cos( alpha_i ),
                radius * math.sin( alpha_i ),
                section_2[i][2]
            ]
        ];

        block_mesh.addBlock( points, n_points, grade );
        block_mesh.addFace(
            "wing", [ points[0], points[1], points[2], points[3] ], "patch"
        );
        block_mesh.addFace(
            "freestream", [ points[4], points[5], points[6], points[7] ], "patch"
        );

    p_0 = [ section_1[0][0], radius, section_1[0][2] ];
    p_1 = [ radius,          radius, section_1[0][2] ];
    p_2 = [ radius,          radius, section_2[0][2] ];
    p_3 = [ section_2[0][0], radius, section_2[0][2] ];
    p_4 = section_1[0];
    p_5 = [ radius, section_1[0][1], section_1[0][2] ];
    p_6 = [ radius, section_2[0][1], section_2[0][2] ];
    p_7 = section_2[0];
    p_8 = section_1[-1];
    p_9 = [ radius, section_1[-1][1], section_1[-1][2] ];
    p_a = [ radius, section_2[-1][1], section_2[-1][2] ];
    p_b = section_2[-1];
    p_c = [ section_1[0][0],  -1.0 * radius, section_1[0][2]  ];
    p_d = [ radius,           -1.0 * radius, section_1[-1][2] ];
    p_e = [ radius,           -1.0 * radius, section_2[-1][2] ];
    p_f = [ section_2[-1][0], -1.0 * radius, section_2[-1][2] ];

    block_mesh.addBlock(
        [ p_0, p_1, p_2, p_3, p_4, p_5, p_6, p_7   ],
        [ n_points[2], n_points[1], n_points[2]    ],
        [ grade[2],    grade[1],    1.0 / grade[2] ]
    );
    block_mesh.addBlock(
        [ p_4, p_5, p_6, p_7, p_8, p_9, p_a, p_b   ],
        [ n_points[2], n_points[1], 1        ],
        [ grade[2],    grade[1],    grade[2] ]
    );
    block_mesh.addBlock(
        [ p_8, p_9, p_a, p_b, p_c, p_d, p_e, p_f ],
        [ n_points[2], n_points[1], n_points[2] ],
        [ grade[2],    grade[1],    grade[2]    ]
    );
    block_mesh.addFace( "wing",       [ p_4, p_7, p_8, p_b ], "patch" );
    block_mesh.addFace( "freestream", [ p_0, p_1, p_2, p_3 ], "patch" );
    block_mesh.addFace( "freestream", [ p_f, p_e, p_d, p_c ], "patch" );
    block_mesh.addFace( "freestream", [ p_1, p_2, p_6, p_5 ], "patch" );
    block_mesh.addFace( "freestream", [ p_d, p_e, p_a, p_9 ], "patch" );
    block_mesh.addFace( "freestream", [ p_9, p_a, p_6, p_5 ], "patch" );

#Creates a winglet section.
def makeWinglet( block_mesh, section, points, radius, n_points, grade ):
    t_e_u = copy.deepcopy( section[ 0                  ]); t_e_u[2] += radius;
    t_e_d = copy.deepcopy( section[ -1                 ]); t_e_d[2] += radius;
    l_e_u = copy.deepcopy( section[ len(section)/2 - 1 ]); l_e_u[2] += radius;
    l_e_d = copy.deepcopy( section[ len(section)/2 + 1 ]); l_e_d[2] += radius;
    addPoint( t_e_u );
    addPoint( t_e_d );
    addPoint( l_e_u );
    addPoint( l_e_d );

################################################################################
## Script ######################################################################
################################################################################

block_mesh = BlockMesh();

# print readAerofoil('constant/aerofoils/NACA_2411');
root   = readAerofoil('constant/aerofoils/NACA_2414', [   0.0, 0.0,   0.0 ], 2.143, 0 );
tip    = readAerofoil('constant/aerofoils/NACA_2411', [ 0.293, 0.0, 4.215 ], 1.050, 0 );

makeWingSection(
    block_mesh,
    root,
    tip,
    15.0,
    [1,1,1], # 30,20
    [1,1,1]  # 1000
);

#makeWinglet(
#    block_mesh["pointlist"],
#    tip,
#    1,
#    5.0,
#    [1,1,1], # 30,20
#    [1,1,1]  # 1000
#);

block_mesh.write();
