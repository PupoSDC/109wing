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
    arclist   = [];

    # Adds a point to a point list only if it does not exist previously.
    def addPoint(self, point ):
        point = [
            float( '%.6f' % point[0] ),
            float( '%.6f' % point[1] ),
            float( '%.6f' % point[2] )
        ];
        if point not in self.pointlist: self.pointlist.append(point)
        return self.pointlist.index(point);

    # Creates one arc for the mesh between 2 points.
    def addArc( self, arc ):
       arc = [self.addPoint(arc[0]),self.addPoint(arc[1]),arc[2]];
       if arc not in self.arclist: self.arclist.append(arc)
       return self.arclist.index(arc);

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
                repr(hex[9][0]) + " " +
                repr(hex[9][1]) + " " +
                repr(hex[9][2]) + ") // " +
                repr( self.hexblocks.index(hex) ) + " \n"
            );

        # Write down the edges
        blockMeshDict.write( ");\n\nedges\n(\n")
        for arc in self.arclist:
            blockMeshDict.write(
                "    arc " +  repr(arc[0]) + " " +  repr(arc[1]) + " ("
                + repr(arc[2][0]) + " "
                + repr(arc[2][1]) + " "
                + repr(arc[2][2]) + ") // " +
                repr( self.arclist.index(arc) ) + " \n"
            );

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
    n_trails = 0;
    ps       = range(0,8);
    sp       = range(0,8);

    for i in range( 0, len(sec_1)/2 ):
        if sec_1[i][0] > d_cut: n_trails += 1;

    for i in range( 0, len(sec_1)/2 ):

        a_i = ( float(i+0) / float(len(sec_1)-1)*1.5 - 0.25 ) * math.pi / 1.0;
        a_j = ( float(i+1) / float(len(sec_1)-1)*1.5 - 0.25 ) * math.pi / 1.0;

        ps[0] = sec_1[i+0];
        ps[1] = sec_1[i+1];
        ps[2] = sec_2[i+1];
        ps[3] = sec_2[i+0];
        ps[4] = [ -radius*math.sin(a_i), radius*math.cos(a_i), sec_1[i+0][2] ];
        ps[5] = [ -radius*math.sin(a_j), radius*math.cos(a_j), sec_1[i+1][2] ];
        ps[6] = [ -radius*math.sin(a_j), radius*math.cos(a_j), sec_2[i+1][2] ];
        ps[7] = [ -radius*math.sin(a_i), radius*math.cos(a_i), sec_2[i+0][2] ];
        sp[0] = sec_1[ -i-2];
        sp[1] = sec_1[ -i-1];
        sp[2] = sec_2[ -i-1];
        sp[3] = sec_2[ -i-2];
        sp[4] = [ -radius*math.sin(a_j), -radius*math.cos(a_j), sec_1[-i-2][2] ];
        sp[5] = [ -radius*math.sin(a_i), -radius*math.cos(a_i), sec_1[-i-1][2] ];
        sp[6] = [ -radius*math.sin(a_i), -radius*math.cos(a_i), sec_2[-i-1][2] ];
        sp[7] = [ -radius*math.sin(a_j), -radius*math.cos(a_j), sec_2[-i-2][2] ];

        # Create main outer wing blocks and free stream boundaries
        block_mesh.addBlock( ps[0:8], points, grade );
        block_mesh.addBlock( sp[0:8], points, grade );
        block_mesh.addFace( "freestream", ps[4:8], "patch" );
        block_mesh.addFace( "freestream", sp[4:8], "patch" );

        # Create inside of the wing blocks (if necessary)
        if( i < n_trails ):
            block_mesh.addBlock(
                [ sp[1], sp[0], sp[3], sp[2], ps[0], ps[1], ps[2], ps[3] ],
                [ points[0], points[1], len(sec_1)/8 ],
                [ 1, 1, 1 ]
            );
            block_mesh.addFace("wing", [ sp[1], sp[0], ps[1], ps[0] ], "wall" );
            block_mesh.addFace("wing", [ sp[3], sp[2], ps[3], ps[2] ], "wall" );

        if( i >= n_trails ):
            block_mesh.addFace("wing", ps[0:4], "wall" );
            block_mesh.addFace("wing", sp[0:4], "wall" );

        if( i == n_trails and i != 0 ):
            block_mesh.addFace("wing", [ sp[2], sp[1], ps[3], ps[0] ], "wall");

        # Add left surfaces if first section
        if sec_1[i][2] == 0.0:
            block_mesh.addFace("left", [ps[0],ps[1],ps[4],ps[5]], "symmetry");
            block_mesh.addFace("left", [sp[0],sp[1],sp[4],sp[5]], "symmetry");

        # Create trailing edge Section
        if( i == 0 ):
            te_1 = [ radius, sp[1][1], sp[1][2] ];
            te_2 = [ radius, sp[2][1], sp[2][2] ];
            ac_1 = [ radius*math.cos(0.4), radius*math.sin(0.4), sp[1][2] ];
            ac_2 = [ radius*math.cos(0.4), radius*math.sin(0.4), sp[2][2] ];
            ac_3 = [ radius*math.cos(0.4),-radius*math.sin(0.4), sp[1][2] ];
            ac_4 = [ radius*math.cos(0.4),-radius*math.sin(0.4), sp[2][2] ];
            block_mesh.addBlock(
                [ ps[0], ps[0], ps[3], ps[3], te_1, ps[4] ,ps[7], te_2 ],
                [ len(sec_1)/8, points[1], points[2] ],
                grade
            );
            block_mesh.addBlock(
                [ sp[1], sp[1], sp[2], sp[2], sp[5], te_1, te_2, sp[6] ],
                [ len(sec_1)/8, points[1], points[2] ],
                grade
            );
            block_mesh.addArc( [ te_1, ps[4], ac_1 ]);
            block_mesh.addArc( [ te_2, ps[7], ac_2 ]);
            block_mesh.addArc( [ te_1, sp[5], ac_3 ]);
            block_mesh.addArc( [ te_2, sp[6], ac_4 ]);
            if sec_1[0][2] == 0.0:
                block_mesh.addFace("left",[ps[0],ps[0],ps[4],te_1],"symmetry");
                block_mesh.addFace("left",[sp[1],sp[1],te_1,sp[5]],"symmetry");
            block_mesh.addFace("freestream",[te_2,te_1,ps[4],ps[7]],"patch");
            block_mesh.addFace("freestream",[te_1,te_2,sp[6],sp[5]],"patch");

#Creates a winglet section.
def makeWinglet( block_mesh, section, radius, points, grade ):

    ps = range(0,8);

    for i in range(0, (len(section))/2 ):
        for j in range(0, points[0] ):

            # calculate  auxiliary values for the first section
            d_x_i =    ( section[-i-1][0] - section[i][0] );
            r_w_i = abs( section[-i-1][1] - section[i][1] ) / 2.0;
            c_w_i =    ( section[-i-1][1] + section[i][1] ) / 2.0;
            a_i = ( float(i+0) / float(len(section)-1)*1.5 - 0.25 ) * math.pi;

            # calculate auxiliary values for the second section
            d_x_j =    ( section[-i-2][0] - section[i+1][0] );
            r_w_j = abs( section[-i-2][1] - section[i+1][1] ) / 2.0;
            c_w_j =    ( section[-i-2][1] + section[i+1][1] ) / 2.0;
            a_j = ( float(i+1) / float(len(section)-1)*1.5 - 0.25 ) * math.pi;

            # displacement indexes
            ind_a = float(j)   / float( points[0] );
            ind_b = float(j+1) / float( points[0] );
            phi_a = math.pi * ind_a; # 0 to 180
            phi_b = math.pi * ind_b; # 0 to 180

            ps[0] = [ #first point in the aerofoil
                section[i][0] + d_x_i * ind_a,
                c_w_i         + r_w_i * math.cos( phi_a ),
                section[i][2] + r_w_i * math.sin( phi_a )
            ];
            ps[1] = [ # second point in the aerofoil
                section[i+1][0] + d_x_j * ind_a,
                c_w_j           + r_w_j * math.cos( phi_a ),
                section[i+1][2] + r_w_j * math.sin( phi_a )
            ];
            ps[2] = [ # extrude the point in the aerofoil
                section[i+1][0] + d_x_j * ind_b,
                c_w_j           + r_w_j * math.cos( phi_b ),
                section[i+1][2] + r_w_j * math.sin( phi_b )
            ];
            ps[3] = [ # next point in extruded aerofoil
                section[i][0] + d_x_i * ind_b,
                c_w_i         + r_w_i * math.cos( phi_b ),
                section[i][2] + r_w_i * math.sin( phi_b )
            ];
            ps[4] = [ # first point in away mesh
                - radius * math.sin( a_i ),
                radius * math.cos( a_i ) * math.cos( phi_a ),
                section[i][2] + radius * math.cos( a_i ) * math.sin( phi_a )
            ];
            ps[5] = [ # second point in away mesh
                - radius * math.sin(a_j),
                radius * math.cos(a_j) * math.cos( phi_a ),
                section[i+1][2] + radius * math.cos(a_j) * math.sin( phi_a )
            ];
            ps[6] = [ # third point in away mesh
                - radius * math.sin(a_j),
                radius * math.cos(a_j) * math.cos( phi_b ),
                section[i+1][2] + radius * math.cos(a_j) * math.sin( phi_b )
            ];
            ps[7] = [ # forth point in away mesh
                - radius * math.sin( a_i ),
                radius * math.cos( a_i ) * math.cos( phi_b ),
                section[i][2] + radius * math.cos( a_i ) * math.sin( phi_b )
            ];

            block_mesh.addBlock( ps, [ 1, points[1], points[2] ] , grade );
            block_mesh.addFace("freestream", [ps[4],ps[5],ps[6],ps[7]], "patch");
            block_mesh.addFace("wing",       [ps[3],ps[2],ps[1],ps[0]], "wall");

            # Create trailing edge Section
            if( i == 0 ):
                t_e = [ radius, ps[0][1], ps[0][2] ];
                a_c_1 = [
                    radius*math.sin(0.4),
                    radius*math.cos(0.4) * math.cos( phi_a ),
                    section[i][2] + radius * math.cos( 0.4 ) * math.sin( phi_a )
                ];
                a_c_2 = [
                    radius*math.sin(0.4),
                    radius*math.cos(0.4) * math.cos( phi_b ),
                    section[i][2] + radius * math.cos( 0.4 ) * math.sin( phi_b )
                ];
                block_mesh.addBlock(
                    [  ps[0], ps[0], ps[3], ps[3], t_e, ps[4], ps[7] ,t_e ],
                    [ len(section)/8, points[1], points[2] ],
                    grade
                );
                block_mesh.addFace("freestream",[t_e,t_e,ps[4],ps[7]],"patch");
                if( j > 0 ):
                    block_mesh.addArc( [ t_e, ps[4], a_c_1 ] );


################################################################################
## Script ######################################################################
################################################################################

block_mesh = BlockMesh();

root  = readAerofoil('constant/aerofoils/NACA/2414', [ -0.958, 0.0,   0.0   ], 2.143 );
tip_a = readAerofoil('constant/aerofoils/NACA/2411', [ -0.470, 0.109, 4.215 ], 1.050 );
tip_b = readAerofoil('constant/aerofoils/NACA/2411', [ -0.465, 0.109, 4.265 ], 1.045 );
tip_c = readAerofoil('constant/aerofoils/NACA/2411', [ -0.461, 0.109, 4.300 ], 1.029 );

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
    [1,15,20],  # 30,20
    [1,1,1000]
);

makeWingSection(
    block_mesh,
    mid,
    tip,
    0.25,
    3.0,
    [1,15,20],  # 30,20
    [1,1,1000]
);

makeWingSection(
    block_mesh,
    tip,
    tip_c,
    0,
    3.0,
    [1,1,20],  # 30,20
    [1,1,1000]
);

makeWinglet( #winglet radius should be 135
   block_mesh,
   tip_c,
   3.0,
   [20,1,20],    # 1,20
   [1,1,1000]  # 1000
);

block_mesh.write();

os.system('blockMesh')

# 1969 - 4215
#   63 - 135
