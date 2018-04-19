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
            response[-1].append( position[2] );
    return response;

# Calculates the intermediate section between 2 sections.
def mixSections( section_1, section_2, percentage ):
    response = [];
    for i in range(0, len(section_1) ):
        response.append([
            section_1[i][0] + (section_2[i][0] - section_1[i][0]) * percentage,
            section_1[i][1] + (section_2[i][1] - section_1[i][1]) * percentage,
            section_1[i][2] + (section_2[i][2] - section_1[i][2]) * percentage
        ]);
    return response;

#distorts a wing section by a linear interpolation of 2 aerofoils.
def distortSection( section_1, section_2, p_start, p_end ):
    response = [];
    length  = section_1[ len(section_1) / 2 + 1 ][0] - section_1[0][0];
    d_start = section_1[ len(section_1) / 2 + 1 ][0] - p_start * length;
    d_end   = section_1[ len(section_1) / 2 + 1 ][0] - p_end   * length;
    for i in range(0, len(section_1) ):
        percentage = 0;
        if( section_1[i][0] > d_start ):
            percentage = (section_1[i][0] - d_start) / (d_end - d_start);
        response.append([
            section_1[i][0] + (section_2[i][0] - section_1[i][0]) * percentage,
            section_1[i][1] + (section_2[i][1] - section_1[i][1]) * percentage,
            section_1[i][2] + (section_2[i][2] - section_1[i][2]) * percentage
        ]);
    return response;

# Creates a wing section
def makeWingSection( block_mesh, sec_1, sec_2, p_cut, radius, points, grade ):
    length   = sec_1[ len(sec_1) / 2 + 1 ][0] - sec_1[0][0];
    d_cut    = sec_1[0][0] + p_cut * length;
    n_leads  = 0; k = 0;
    n_trails = 0;
    ps       = range(0,8);
    sp       = range(0,8);

    for i in range( 0, len(sec_1)/2 ):
        if sec_1[i][0] < 0.0:   n_leads  += 1;
        if sec_1[i][0] > d_cut: n_trails += 1;

    for i in range( 0, len(sec_1)/2 ):
        ps[0] = sec_1[i+0];
        ps[1] = sec_1[i+1];
        ps[2] = sec_2[i+1];
        ps[3] = sec_2[i+0];
        ps[4] = [ sec_1[i+0][0], radius,  sec_1[i+0][2]   ];
        ps[5] = [ sec_1[i+1][0], radius,  sec_1[i+1][2]   ];
        ps[6] = [ sec_2[i+1][0], radius,  sec_2[i+1][2]   ];
        ps[7] = [ sec_2[i+0][0], radius,  sec_2[i+0][2]   ];
        sp[0] = sec_1[ -i-2];
        sp[1] = sec_1[ -i-1];
        sp[2] = sec_2[ -i-1];
        sp[3] = sec_2[ -i-2];
        sp[4] = [ sec_1[-i-2][0], -radius, sec_1[-i-2][2] ];
        sp[5] = [ sec_1[-i-1][0], -radius, sec_1[-i-1][2] ];
        sp[6] = [ sec_2[-i-1][0], -radius, sec_2[-i-1][2] ];
        sp[7] = [ sec_2[-i-2][0], -radius, sec_2[-i-2][2] ];

        if sec_1[i+0][0] < 0.0:
            a_i = ( float(k) / n_leads ) * math.pi / 2.0; k += 1;
            ps[4][0] = -radius*math.sin(a_i); ps[4][1] =  radius*math.cos(a_i);
            ps[7][0] = -radius*math.sin(a_i); ps[7][1] =  radius*math.cos(a_i);
            sp[5][0] = -radius*math.sin(a_i); sp[5][1] = -radius*math.cos(a_i);
            sp[6][0] = -radius*math.sin(a_i); sp[6][1] = -radius*math.cos(a_i);

        if sec_1[i+1][0] < 0.0:
            a_i = ( float(k) / n_leads ) * math.pi / 2.0;
            ps[5][0] = -radius*math.sin(a_i); ps[5][1] =  radius*math.cos(a_i);
            ps[6][0] = -radius*math.sin(a_i); ps[6][1] =  radius*math.cos(a_i);
            sp[4][0] = -radius*math.sin(a_i); sp[4][1] = -radius*math.cos(a_i);
            sp[7][0] = -radius*math.sin(a_i); sp[7][1] = -radius*math.cos(a_i);

        # Create main outer wing blocks and free stream boundaries
        block_mesh.addBlock( ps[0:8], points, grade );
        block_mesh.addBlock( sp[0:8], points, grade );
        block_mesh.addFace( "freestream", ps[4:8], "patch" );
        block_mesh.addFace( "freestream", sp[4:8], "patch" );

        # Create inside of the wing blocks (if necessary)
        if( i < n_trails ):
            block_mesh.addBlock(
                [ sp[1], sp[0], sp[3], sp[2], ps[0], ps[1], ps[2], ps[3] ],
                [ points[0], points[1], points[2] ],
                [ 1 , 1 , 1 ]
            );
            block_mesh.addFace("wing", [ sp[1], sp[0], ps[1], ps[0] ], "wall" );
            block_mesh.addFace("wing", [ sp[3], sp[2], ps[3], ps[2] ], "wall" );
        if( i >= n_trails ):
            block_mesh.addFace("wing", ps[0:4], "wall" );
            block_mesh.addFace("wing", sp[0:4], "wall" );

        if( i == n_trails and i != 0 ):
            block_mesh.addFace("wing", [ sp[2], sp[1], ps[3], ps[0] ], "wall");

        if sec_1[i][2] == 0.0:
            block_mesh.addFace("left", [ps[0],ps[1],ps[4],ps[5]], "symmetry");
            block_mesh.addFace("left", [sp[0],sp[1],sp[4],sp[5]], "symmetry");

        # Create trailing edge Sections
        if( i == 0 ):
            tp_0 = [ radius, radius     , sec_1[0][2] ];
            tp_1 = [ radius, sec_1[0][1], sec_1[0][2] ];
            tp_2 = [ radius, -radius    , sec_1[0][2] ];
            te_0 = [ radius, radius     , sec_2[0][2] ];
            te_1 = [ radius, sec_2[0][1], sec_2[0][2] ];
            te_2 = [ radius, -radius    , sec_2[0][2] ];
            block_mesh.addBlock(
                [ ps[0], ps[3], te_1, tp_1, ps[4], ps[7], te_0, tp_0 ],
                [ points[1], points[2], points[2] ],
                [ grade[1] , grade[2] , grade[2]  ]
            );
            block_mesh.addBlock(
                [ sp[5], sp[6], te_2, tp_2, sp[1], sp[2], te_1, tp_1 ],
                [ points[1], points[2], points[2] ],
                [ grade[1] , grade[2] , 1.0/grade[2]  ]
            );

            if sec_1[0][2] == 0.0:
                block_mesh.addFace("left", [ps[0],ps[4],tp_0,tp_1], "symmetry");
                block_mesh.addFace("left", [tp_2,sp[5],sp[1],tp_1], "symmetry");

            block_mesh.addFace("freestream", [ps[7],te_0,tp_0,ps[4]], "patch");
            block_mesh.addFace("freestream", [te_1,tp_1,tp_0,te_0],   "patch");
            block_mesh.addFace("freestream", [te_2,tp_2,tp_1,te_1],   "patch");
            block_mesh.addFace("freestream", [sp[6],sp[5],tp_2,te_2], "patch");


##Creates a winglet section.
#def makeWinglet( block_mesh, section, radius, points, grade ):
#    for i in range(1, (len(section)-1)/2 ):
#        for j in range(0, points[0] ):
#
#            # calculate  auxiliary values for the first section
#            d_x_i =    ( section[-2-i][0] - section[i-1][0] );
#            r_w_i = abs( section[-2-i][1] - section[i-1][1] ) / 2.0;
#            c_w_i =    ( section[-2-i][1] + section[i-1][1] ) / 2.0;
#            alp_i = 2.0 * math.pi * float(i+0) / len(section);
#
#            # calculate auxiliary values for the second section
#            d_x_j =    ( section[-3-i][0] - section[i+0][0] );
#            r_w_j = abs( section[-3-i][1] - section[i+0][1] ) / 2.0;
#            c_w_j =    ( section[-3-i][1] + section[i+0][1] ) / 2.0;
#            alp_j = 2.0 * math.pi * float(i+1) / len(section);
#
#            # displacement indexes
#            ind_a = float(j) / float( points[0] );
#            ind_b = float(j+1) / float( points[0] );
#            phi_a = math.pi * ind_i;
#            phi_b = math.pi * ind_j;
#
#            ps = [
#                [ # firs point in the aerofoil
#                    section[i][0] + d_x_i * ind_i,
#                    c_w_i         + r_w_i * math.cos( phi_i ),
#                    section[i][2] + r_w_i * math.sin( phi_i )
#                ],
#                [ # second point in the aerofoil
#                    section[i+1][0] + d_x_j * ind_i,
#                    c_w_j           + r_w_j * math.cos( phi_i ),
#                    section[i+1][2] + r_w_j * math.sin( phi_i )
#                ],
#                [ # extrude the point in the aerofoil
#                    section[i+1][0] + d_x_j * ind_j,
#                    c_w_j           + r_w_j * math.cos( phi_j ),
#                    section[i+1][2] + r_w_j * math.sin( phi_j )
#                ],
#                [ # next point in extruded aerofoil
#                    section[i][0] + d_x_j * ind_j,
#                    c_w_j           + r_w_j * math.cos( phi_j ),
#                    section[i][2] + r_w_j * math.sin( phi_j )
#                ],
#                [ # first point in away mesh
#                    radius * math.sin( alp_i ),
#                    radius * math.cos( alp_i ) * math.cos( phi_i ),
#                    radius * math.cos( alp_i ) * math.sin( phi_i )
#                ],
#                [ # second point in away mesh
#                    radius * math.sin( alp_j ),
#                    radius * math.cos( alp_j ) * math.cos( phi_i ),
#                    radius * math.cos( alp_j ) * math.sin( phi_i )
#                ],
#                [ # third point in away mesh
#                    radius * math.sin( alp_j ),
#                    radius * math.cos( alp_j ) * math.cos( phi_j ),
#                    radius * math.cos( alp_j ) * math.sin( phi_j )
#                ],
#                [ # forth point in away mesh
#                    radius * math.sin( alp_i ),
#                    radius * math.cos( alp_i ) * math.cos( phi_j ),
#                    radius * math.cos( alp_i ) * math.sin( phi_j )
#                ]
#            ];
#            block_mesh.addPoint( ps[0] );
#            block_mesh.addPoint( ps[1] );
#            block_mesh.addPoint( ps[2] );
#            block_mesh.addPoint( ps[3] );
#            #block_mesh.addPoint( ps[4] );
#            #block_mesh.addPoint( ps[5] );
#            #block_mesh.addPoint( ps[6] );
#            #block_mesh.addPoint( ps[7] );
#            #block_mesh.addBlock( ps, [ 1, points[1], points[2] ] , grade );


################################################################################
## Script ######################################################################
################################################################################

block_mesh = BlockMesh();

root  = readAerofoil('constant/aerofoils/NACA/2414', [ -0.958, 0.0,   0.0   ], 2.143 );
tip_a = readAerofoil('constant/aerofoils/NACA/2411', [ -0.470, 0.109, 4.215 ], 1.050 );
tip_b = readAerofoil('constant/aerofoils/NACA/2411', [ -0.470, 0.109, 4.265 ], 1.050 );
tip_c = readAerofoil('constant/aerofoils/NACA/2411', [ -0.470, 0.109, 4.270 ], 1.050 );

mid_a = mixSections(root,tip_a,0.544);
mid_b = mixSections(root,tip_a,0.560);
mid   = distortSection( mid_a, mid_b, 0.75, 1.0 );
tip   = distortSection( tip_a, tip_b, 0.75, 1.0 );

makeWingSection(
    block_mesh,
    root,
    mid,
    0.0,
    3.0,
    [1,30,20],  # 30,20
    [1,1,1000]  # 1000
);

makeWingSection(
    block_mesh,
    mid,
    tip,
    0.25,
    3.0,
    [1,30,20],  # 30,20
    [1,1,1000]  # 1000
);

makeWingSection(
    block_mesh,
    tip,
    tip_c,
    0.0,
    3.0,
    [1,1,20],   # 30,20
    [1,1,1000]  # 1000
);

#akeWinglet(
#   block_mesh,
#   tip,
#   1.0,
#   [20,1,20],    # 1,20
#   [1,1,1000]  # 1000
#;

block_mesh.write();

os.system('blockMesh')
