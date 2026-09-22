# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from abaqusConstants import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import csv

import numpy as np
import os

from odbAccess import openOdb

# pip install git+https://github.com/lcharleux/abapy.git

BEAM_SMALL_LENGTH = 35 #mm
BEAM_LARGE_LENGTH = 45 #mm 45
NODE_SPACING = 5

# Node Size
START_OF_A = -BEAM_SMALL_LENGTH - BEAM_LARGE_LENGTH/2
END_OF_C = BEAM_SMALL_LENGTH + BEAM_LARGE_LENGTH/2
BEAM_LENGTH = (END_OF_C-START_OF_A)
RADIUS = 4
DENSITY = 470

# size ratio percentage of one of the outer members (double shear) 
number_of_nodes = round(BEAM_LENGTH/NODE_SPACING + 1)


RATIO = BEAM_SMALL_LENGTH/(BEAM_SMALL_LENGTH + BEAM_LARGE_LENGTH + BEAM_SMALL_LENGTH)
BREAK_A = round(number_of_nodes*RATIO)
BREAK_B = round(number_of_nodes*(1-RATIO))

# Nodal parameters
TOTAL_FORCE = 600 * BREAK_A * (number_of_nodes - 2*BREAK_A)
# abaqus rounds forces to 3 decimal places so you need an int force amount

# force applied to centre 
SPRING_LENGTH = 10

model = mdb.models['Model-1']
assembly = model.rootAssembly


def create_line_part():

    sketch = model.ConstrainedSketch(
        name='__profile__',
        sheetSize=200.0)

    i = START_OF_A
    while i < END_OF_C:
        sketch.Line(
            point1=(i, 0.0),
            point2=(i + NODE_SPACING, 0.0))
        i = i + NODE_SPACING


    part = model.Part(
        dimensionality=TWO_D_PLANAR,
        name='Part-1',
        type=DEFORMABLE_BODY)

    part.BaseWire(sketch=sketch)



    del model.sketches['__profile__']

def assign_properties():
    model.Material(name='Steel')
    model.materials['Steel'].Elastic(table=((200000.0, 0.3), ))  # E in MPa, matches your existing units
    model.materials['Steel'].Plastic(table=(
        (200000, 0.000), 
        (28000, 0.005),
        (32000, 0.020),
        (36000, 0.050)
    ))
    model.CircularProfile(name='Profile-1', r=RADIUS)
    model.BeamSection(
        beamSectionOffset=(0.0, 0.0), 
        consistentMassMatrix=False, 
        integration=DURING_ANALYSIS, 
        material='Steel', 
        name='Section-1', 
        poissonRatio=0.0, 
        profile='Profile-1', 
        temperatureVar=LINEAR)

    part = model.parts['Part-1']

    part.Set(
        edges=part.edges[:],
        name='Set-1')

    model.parts['Part-1'].SectionAssignment(
        offset=0.0, 
        offsetField='', 
        offsetType=MIDDLE_SURFACE, 
        region=model.parts['Part-1'].sets['Set-1'], 
        sectionName='Section-1', 
        thicknessAssignment=FROM_SECTION)



def create_part_instances():
    assembly.DatumCsysByDefault(CARTESIAN)
    assembly.Instance(dependent=ON, name='centre_timber', part=model.parts['Part-1'])
    assembly.Instance(dependent=ON, name='right_timber', part=model.parts['Part-1'])
    assembly.Instance(dependent=ON, name='left_timber', part=model.parts['Part-1'])


    #assembly.LinearInstancePattern(direction1=(1.0, 0.0, 0.0), direction2=(0.0, 1.0, 0.0), instanceList=('centre_timber', ), number1=3, number2=1, spacing1=25.0, spacing2=100.0)
    assembly.rotate(angle=90.0, axisDirection=(0.0, 0.0, 1.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('left_timber', 'right_timber', 'centre_timber'))

    assembly.Instance(dependent=ON, name='steel_dowel', part=model.parts['Part-1'])

    
    distance = BEAM_LENGTH/2 + SPRING_LENGTH
    assembly.translate(instanceList=('right_timber',), vector=(BEAM_LENGTH/4, -distance, 0.0))
    assembly.translate(instanceList=('left_timber',), vector=(-BEAM_LENGTH/4, -distance, 0.0))
    assembly.translate(instanceList=('centre_timber', ), vector=(0.0, distance, 0.0))







def create_connector(
        name,
        instance_1,
        vertex_1_index,
        instance_2,
        vertex_2_index,

        # ------------------------------------------------------
        # Connector type
        # ------------------------------------------------------
        translational_type=CARTESIAN,
        rotational_type=ROTATION,

        # # ------------------------------------------------------
        # # Elastic behaviour direction 1 (HORIZONTAL)
        # # ------------------------------------------------------
        # components_1 =(1,),
        force_displacement_table_1 = (),
        # behavior_1 = NONLINEAR,
        # coupling_1 = UNCOUPLED,

        # # ------------------------------------------------------
        # # Elastic behaviour direction 2 (VERTICAL)
        # # ------------------------------------------------------
        # components_2 = (2,),
        force_displacement_table_2 = (),
        # behavior_2 = NONLINEAR,
        # coupling_2 = UNCOUPLED,

        # ------------------------------------------------------
        # Connector orientation
        # ------------------------------------------------------
        local_csys_1=None,
        local_csys_2=None,
        axis_1=AXIS_1,
        angle_1=0.0,
        axis_2=AXIS_1,
        angle_2=0.0,

        # ------------------------------------------------------
        # Section options
        # ------------------------------------------------------
        integration=UNSPECIFIED,
        extrapolation=CONSTANT,
        regularize=ON,
        default_tolerance=ON,
        regularization=0.03
    ):
    """
    Create a connector between two existing vertices.

    Parameters
    ----------
    name : str
        Base name used for the connector set and section.

    instance_1, instance_2 : str
        Names of the two instances.

    vertex_1_index, vertex_2_index : int
        Vertex indices on the two instances.

    translational_type : Abaqus SymbolicConstant
        e.g. CARTESIAN, AXIAL, SLOT, LINK, etc.

    rotational_type : Abaqus SymbolicConstant
        e.g. NONE, ROTATION, REVOLUTE, etc.

    components : tuple
        Connector components to which the elastic behaviour applies.

        Examples:
            (1,)       -> component 1
            (2,)       -> component 2
            (1, 2, 3) -> components 1, 2 and 3

    force_displacement_table : tuple
        Nonlinear force-displacement data.

        Example:
            (
                (-5000.0, -5.0),
                (-1000.0, -1.0),
                (0.0, 0.0),
                (1000.0, 1.0),
                (5000.0, 5.0),
            )

    behavior : Abaqus SymbolicConstant
        NONLINEAR, LINEAR or RIGID.

    coupling : Abaqus SymbolicConstant
        UNCOUPLED, COUPLED_POSITION or COUPLED_MOTION
        for nonlinear behaviour.

    local_csys_1, local_csys_2 :
        Optional datum coordinate systems defining connector
        orientations at the two connector points.

    Returns
    -------
    dict
        Information about the created connector.
    """

    # ==========================================================
    # Names
    # ==========================================================

    set_name = name + '-Set'
    section_name = name + '-Section'

    # ==========================================================
    # Get vertices
    # ==========================================================

    vertex_1 = assembly.instances[
        instance_1
    ].vertices[
        vertex_1_index
    ]

    vertex_2 = assembly.instances[
        instance_2
    ].vertices[
        vertex_2_index
    ]

    # ==========================================================
    # Get coordinates
    # ==========================================================

    p1 = vertex_1.pointOn[0]
    p2 = vertex_2.pointOn[0]

    midpoint = (
        0.5 * (p1[0] + p2[0]),
        0.5 * (p1[1] + p2[1]),
        0.5 * (p1[2] + p2[2])
    )

    # ==========================================================
    # Create connector wire
    # ==========================================================

    assembly.WirePolyLine(
        points=((vertex_1, vertex_2),),
        mergeType=SEPARATE,
        meshable=OFF
    )

    assembly.regenerate()

    # ==========================================================
    # Find newly-created connector edge
    # ==========================================================

    connector_edge = assembly.edges.findAt(
        (midpoint,)
    )

    # ==========================================================
    # Create connector set
    # ==========================================================

    connector_set = assembly.Set(
        name=set_name,
        edges=connector_edge
    )

    # ==========================================================
    # Create connector elasticity
    # ==========================================================

    # behaviorOptions=(
    #     ConnectorElasticity(
    #         behavior=NONLINEAR, coupling=COUPLED_POSITION, 
    #         table=((0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0)), 
    #         independentComponents=(1, 2, 5), components=(1, )), 
    #     ConnectorElasticity(
    #         behavior=NONLINEAR, 
    #         coupling=COUPLED_POSITION, 
    #         table=((0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0)), 
    #         independentComponents=(1, 2, 4), 
    #         components=(2, )
    #         ))



    # connector_elasticity_d1 = ConnectorElasticity(
    #     behavior=NONLINEAR,
    #     coupling=COUPLED_POSITION,
    #     table=((0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0)),
    #     components=(1, ),
    #     independentComponents=(1, 2, 5)
    # )

    connector_elasticity_d2 = ConnectorElasticity(
        behavior=NONLINEAR,
        coupling=COUPLED_POSITION,
        table=((0.0, 0.0, 0.0, 0.0), 
               (1.0, 1.0, 0, 0),
               (1.0, 0, 1.0, 0),
               (1.0, 0, 0, 1.0),
               ),
        components=(2, ),
        independentComponents=(1, 2, 6)

    )



    # ==========================================================
    # Create connector section
    # ==========================================================

    connector_section = model.ConnectorSection(
        name=section_name,
        translationalType=translational_type,
        rotationalType=rotational_type,
        behaviorOptions=(connector_elasticity_d2,),
        integration=integration,
        extrapolation=extrapolation,
        regularize=regularize,
        defaultTolerance=default_tolerance,
        regularization=regularization
    )

    # ==========================================================
    # Assign section
    # ==========================================================

    assembly.SectionAssignment(
        region=connector_set,
        sectionName=section_name
    )

    assembly.regenerate()

    # ==========================================================
    # Connector orientation
    # ==========================================================

    if local_csys_1 is not None or local_csys_2 is not None:

        assembly.ConnectorOrientation(
            region=connector_set,
            localCsys1=local_csys_1,
            axis1=axis_1,
            angle1=angle_1,
            localCsys2=local_csys_2,
            axis2=axis_2,
            angle2=angle_2
        )

    # ==========================================================
    # Return useful information
    # ==========================================================

    return {
        'name': name,
        'set_name': set_name,
        'section_name': section_name,
        'vertex_1': vertex_1,
        'vertex_2': vertex_2,
        'edge': connector_edge,
        'section': connector_section,
        'elasticity': connector_elasticity_d2,
    }





def spring_calc_0(l, rho, u, a = 3):
    # l = supporting length
    # d = diameter
    # rho = density
    # u = displacement
    # richard abbot curves
    #EC5 = 0.082 * (1 - 0.01 * RADIUS * 2) * rho
    f_h_inter = 0.1253 * rho - 20.32
    k_f_el_0 = 0.1374 * rho - 12.9
    k_f_pl_0 = 0.0047 * rho - 2
    f_h_0 = (k_f_el_0 - k_f_pl_0) * u / ((1 + ((k_f_el_0 - k_f_pl_0) * u / f_h_inter) ** a) ** (1 / a)) + k_f_pl_0 * u
    return f_h_0 * l * RADIUS * 2


def spring_calc_90(l, rho, u, a = 3):
    # l = supporting length
    # d = diameter
    # rho = density
    # u = displacement
    # richard abbot curves
    #EC5 = 0.082 * (1 - 0.01 * RADIUS * 2) * rho
    f_h_inter = 0.1108 * rho - 28.47
    k_f_el_90 = 0.0922 * rho - 18.2
    k_f_pl_90 = 0.0084 * rho - 2.21
    f_h_0 = (k_f_el_90 - k_f_pl_90) * u / ((1 + ((k_f_el_90 - k_f_pl_90) * u / f_h_inter) ** a) ** (1 / a)) + k_f_pl_90 * u
    return f_h_0 * l * RADIUS * 2


def apply_nodal_force(
        instance_name,
        set_name,
        load_name,
        node_index,
        force,
        step_name='Step-1'):

    instance = assembly.instances[instance_name]

    # ----------------------------------------------------------
    # Create set
    # ----------------------------------------------------------
    vertex = instance.vertices[node_index:node_index + 1]
    assembly.Set(
        name=set_name,
        vertices=(vertex,)
    )

    # ----------------------------------------------------------
    # Apply load
    # ----------------------------------------------------------

    model.ConcentratedForce(
        name=load_name,
        createStepName=step_name,
        region=assembly.sets[set_name],
        cf2=force,
        distributionType=UNIFORM,
        field='',
        localCsys=None
    )


def mesh_part():
    model.parts['Part-1'].seedPart(deviationFactor=0.2, minSizeFactor=0.1, size=1.0)
    model.parts['Part-1'].generateMesh()
    assembly.regenerate()



def set_boundary_conditions():
    assembly.Set(name='right_timber', vertices= assembly.instances['right_timber'].vertices[number_of_nodes - 1 : number_of_nodes])
    model.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, 
                        name= 'right_timber', region=assembly.sets['right_timber'], 
                        u1=0.0, 
                        u2=0.0, 
                        ur3=0.0)
    
    assembly.Set(name='left_timber', vertices= assembly.instances['left_timber'].vertices[number_of_nodes - 1 : number_of_nodes])
    model.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, 
                         name= 'left_timber', region=assembly.sets['left_timber'], 
                         u1=0.0, 
                         u2=0.0, 
                         ur3=0.0)
    
    assembly.Set(name='centre_timber', vertices= assembly.instances['centre_timber'].vertices[number_of_nodes - 1 : number_of_nodes])
    model.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, 
                         name='centre_timber', region=assembly.sets['centre_timber'], 
                         u1=0.0, 
                         u2=UNSET, 
                         ur3=0.0)

    centre = round(number_of_nodes/2)
    assembly.Set(name='steel_dowel', vertices=assembly.instances['steel_dowel'].vertices[centre - 1: centre])
    model.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                        'steel_dowel', region=assembly.sets['steel_dowel'], 
                        u1=UNSET, 
                        u2=UNSET, 
                        ur3=UNSET)
    assembly.regenerate()


def apply_load():
    for i in range(0, BREAK_A):
            apply_nodal_force(
                instance_name='steel_dowel',
                set_name='Load-Set-'+str(i),
                load_name='Load-'+str(i),
    
                node_index = i,
    
                force=-TOTAL_FORCE / BREAK_A / 2
                )
    
    for i in range(BREAK_A, BREAK_B):
            apply_nodal_force(
                instance_name='steel_dowel',
                set_name='Load-Set-'+str(i),
                load_name='Load-'+str(i),
    
                node_index = i,
    
                force=TOTAL_FORCE / (number_of_nodes - 2*BREAK_A)
                )
        
    for i in range(BREAK_B, number_of_nodes):
            apply_nodal_force(
                instance_name='steel_dowel',
                set_name='Load-Set-'+str(i),
                load_name='Load-'+str(i),
    
                node_index = i,
    
                force=-TOTAL_FORCE / BREAK_A / 2
                )


            
def set_step():
    model.StaticStep(initialInc=0.01, maxInc=0.5, maxNumInc=10000, 
    minInc=1e-11, name='Step-1', previous='Initial', nlgeom = OFF)
    #model['Model-1'].steps['Step-1'].setValues(nlgeom=ON)

def make_deflection_curve():
    DEFLECTION_MAX = 20
    DEFLECTION_MIN = -20
    INC = 0.1
    vals = np.arange(DEFLECTION_MIN, DEFLECTION_MAX, INC)
    displacement_table_0 = ()
    for i in vals:
        displacement_table_0 = displacement_table_0 + ((np.sign(i)*spring_calc_0(NODE_SPACING, DENSITY, abs(i)), i),)

    displacement_table_90 = ()
    for i in vals:
        displacement_table_90 = displacement_table_90 + ((np.sign(i)*spring_calc_0(0, DENSITY, abs(i)), i),)

    return (displacement_table_0, displacement_table_90)


def main_double_shear():
    create_line_part()

    assign_properties()

    create_part_instances()

    set_step()

    standard_displacement_table = make_deflection_curve()

    for i in range(0, BREAK_A):
        create_connector(
            name='Spring-' + str(i),
        
            instance_1='steel_dowel',
            vertex_1_index=i,

            instance_2='left_timber',
            vertex_2_index=number_of_nodes - 1,

            force_displacement_table_1=standard_displacement_table[1],
            force_displacement_table_2=standard_displacement_table[0]
            )
    
    for i in range(BREAK_A, BREAK_B):
        create_connector(
            name='Spring-' + str(i),

            instance_1='steel_dowel',
            vertex_1_index=i,

            instance_2='centre_timber',
            vertex_2_index=0,

            force_displacement_table_1=standard_displacement_table[1],
            force_displacement_table_2=standard_displacement_table[0]
            )

    for i in range(BREAK_B, number_of_nodes):
        create_connector(
            name='Spring-' + str(i),

            instance_1='steel_dowel',
            vertex_1_index=i,

            instance_2='right_timber',
            vertex_2_index=number_of_nodes - 1,

            force_displacement_table_1=standard_displacement_table[1],
            force_displacement_table_2=standard_displacement_table[0]
            )



    assembly.Set(name='Set-7', 
                 vertices=assembly.instances['steel_dowel'].vertices.getSequenceFromMask(('[#3fff ]', ), ))

    #apply_load()



    mesh_part()

    set_boundary_conditions()

    assembly.regenerate()

    model.parts['Part-1'].Set(edges=model.parts['Part-1'].edges[:], name='Set-2')
    model.parts['Part-1'].assignBeamSectionOrientation(method=N1_COSINES, n1=(0.0, 0.0, -1.0), region=model.parts['Part-1'].sets['Set-2'])

    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, name=
        'Job-1', nodalOutputPrecision=SINGLE, queue=None, resultsFormat=ODB, 
        scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

 

    model.HistoryOutputRequest(
    name='H-Output-Reactions',
    createStepName='Step-1',
    region=assembly.sets['left_timber'],   # or whichever support set
    variables=('RF1', 'RF2', 'RM3')
    )

    model.HistoryOutputRequest(
    name='H-Output-Displacement',
    createStepName='Step-1',
    region=assembly.sets['centre_timber'],   # or whichever support set
    variables=("U2",)
    )



    assembly.regenerate()
    model.rootAssembly.regenerate()


    mdb.models['Model-1'].rootAssembly.Set(name='Set-30', vertices=
        mdb.models['Model-1'].rootAssembly.instances['steel_dowel'].vertices.getSequenceFromMask(
        ('[#1800 ]', ), ))
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-5', region=mdb.models['Model-1'].rootAssembly.sets['Set-30'], u1=UNSET, 
        u2=15.0, ur3=UNSET)
    mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
    mdb.jobs['Job-1'].waitForCompletion()


    odb = openOdb('Job-1.odb')
    try:
        step = odb.steps['Step-1']

        # History Region key is usually something like 'Node ASSEMBLY.<node label>'
        # or for connectors, 'Element ASSEMBLY.<element label>'
        for region_key in step.historyRegions.keys():
            print(region_key)
            #print("frog")

        region = step.historyRegions['Node LEFT_TIMBER.' + str(str(number_of_nodes))]
        #print(region.historyOutputs)
        rf2_data = region.historyOutputs['RF2'].data   # list of (time, value) tuples

        region = step.historyRegions['Node CENTRE_TIMBER.' + str(number_of_nodes)]
        #print(region.historyOutputs)
        u2_data = region.historyOutputs['U2'].data   # list of (time, value) tuples

        print(rf2_data)
        print(u2_data)
        data = []
        for x, y in zip(rf2_data, u2_data):
            data.append({"rf2":x[1]*-2, "u2":y[1]})    
    

        with open("C:\\Users\\william\\aba_data\\test.csv", 'w', newline='') as csvfile:
            fieldnames = ['u2', 'rf2']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)
        # region = step.historyRegions['Node ASSEMBLY.116']   # example key, copy exact from above
        # rf2_data = region.historyOutputs['RF2'].data   # list of (time, value) tuples

        # for time, value in rf2_data:
        #     print(time, value)
    
    finally:
        odb.close()





def spring_test():

    create_line_part()

    assign_properties()

    create_part_instances()

    set_step()


    # create_connector(
    #     name='Spring-10',

    #     instance_1='steel_dowel',
    #     vertex_1_index=10,

    #     instance_2='centre_timber',
    #     vertex_2_index=0,

    #     force_displacement_table=standard_displacement_table)

    centre = round(number_of_nodes/2)
    create_connector(
        name='Spring-C',

        instance_2='steel_dowel',
        vertex_2_index=centre,

        instance_1='centre_timber',
        vertex_1_index=0
        )

    assembly.Set(name='Set-7', 
                 vertices=assembly.instances['steel_dowel'].vertices.getSequenceFromMask(('[#3fff ]', ), ))

        
    mesh_part()

    set_boundary_conditions()


    
    model.FieldOutputRequest(createStepName='Step-1', name=
        'F-Output-3', position=INTEGRATION_POINTS, rebar=EXCLUDE, region=
        mdb.models['Model-1'].rootAssembly.sets['Spring-C-Set'], sectionPoints=
        DEFAULT, variables=('CTF', 'CEF', 'CU', 'CUE', 'CUP'))

    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, name=
        'Job-1', nodalOutputPrecision=SINGLE, queue=None, resultsFormat=ODB, 
        scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

    model.parts['Part-1'].Set(edges=model.parts['Part-1'].edges[:], name='Set-2')
    model.parts['Part-1'].assignBeamSectionOrientation(method=N1_COSINES, n1=(0.0, 0.0, -1.0), region=model.parts['Part-1'].sets['Set-2'])


    assembly.regenerate()

    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name='Part-2', type=
    DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-2'].ReferencePoint(point=(0.0, 0.0, 0.0))
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-2-1', 
    part=mdb.models['Model-1'].parts['Part-2'])



    """
    mdb.models['Model-1'].ConnectorSection(name='ConnSect-267', rotationalType=
        EULER, translationalType=CARTESIAN)
    mdb.models['Model-1'].sections['ConnSect-267'].setValues(behaviorOptions=(
        ConnectorElasticity(behavior=NONLINEAR, coupling=COUPLED_POSITION, table=((
        0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0)), independentComponents=(1, 2, 
        5), components=(1, )), ConnectorElasticity(behavior=NONLINEAR, 
        coupling=COUPLED_POSITION, table=((0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 
        1.0)), independentComponents=(1, 2, 4), components=(2, ))))
    mdb.models['Model-1'].sections['ConnSect-267'].behaviorOptions[0].ConnectorOptions(
        )
    mdb.models['Model-1'].sections['ConnSect-267'].behaviorOptions[1].ConnectorOptions(
        )
    mdb.models['Model-1'].rootAssembly.ReferencePoint(point=
        mdb.models['Model-1'].rootAssembly.instances['steel_dowel'].InterestingPoint(
        mdb.models['Model-1'].rootAssembly.instances['steel_dowel'].edges[12], 
        MIDDLE))
    mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='RP-2', toName=
        'RP-2')
    mdb.models['Model-1'].rootAssembly.DatumCsysByThreePoints(coordSysType=
        CARTESIAN, isZ=True, origin=
        mdb.models['Model-1'].rootAssembly.instances['centre_timber'].vertices[0], 
        point1=
        mdb.models['Model-1'].rootAssembly.instances['steel_dowel'].InterestingPoint(
        mdb.models['Model-1'].rootAssembly.instances['steel_dowel'].edges[12], 
        MIDDLE))
    mdb.models['Model-1'].rootAssembly.WirePolyLine(mergeType=IMPRINT, meshable=
        False, points=((
        mdb.models['Model-1'].rootAssembly.instances['centre_timber'].vertices[0], 
        mdb.models['Model-1'].rootAssembly.referencePoints[92]), ))
    mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Wire-26', 
        toName='Wire-26')
    mdb.models['Model-1'].rootAssembly.Set(edges=
        mdb.models['Model-1'].rootAssembly.edges.getSequenceFromMask(('[#1 ]', ), )
        , name='Wire-26-Set-1')
    mdb.models['Model-1'].rootAssembly.SectionAssignment(region=
        mdb.models['Model-1'].rootAssembly.sets['Wire-26-Set-1'], sectionName=
        'Spring-0-Section')
    mdb.models['Model-1'].rootAssembly.sectionAssignments[25].getSet()
    mdb.models['Model-1'].rootAssembly.ConnectorOrientation(localCsys1=
        mdb.models['Model-1'].rootAssembly.datums[93], region=
        mdb.models['Model-1'].rootAssembly.allSets['Wire-26-Set-1'])"""
    




    """
    model.ConnectorSection(name='ConnSect-26', rotationalType=EULER, translationalType=CARTESIAN)

    model.sections['ConnSect-26'].setValues(behaviorOptions=(ConnectorElasticity(behavior=NONLINEAR, coupling=COUPLED_POSITION, table=((
        0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0)), independentComponents=(1, 2, 4), components=(1, )), ConnectorElasticity(behavior=NONLINEAR, 
        coupling=COUPLED_POSITION, table=((0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0)), independentComponents=(1, 2, 4), components=(2, ))))

    model.sections['ConnSect-26'].behaviorOptions[0].ConnectorOptions()

    model.sections['ConnSect-26'].behaviorOptions[1].ConnectorOptions()

    model.rootAssembly.ReferencePoint(point=model.rootAssembly.instances['steel_dowel'].InterestingPoint(model.rootAssembly.instances['steel_dowel'].edges[9], MIDDLE))

    model.rootAssembly.features.changeKey(fromName='RP-1', toName='RP-1')

    model.rootAssembly.DatumCsysByThreePoints(coordSysType=CARTESIAN, isZ=True, origin=model.rootAssembly.instances['steel_dowel'].InterestingPoint(model.rootAssembly.instances['steel_dowel'].edges[9], MIDDLE), point1=model.rootAssembly.instances['centre_timber'].vertices[0])

    model.rootAssembly.WirePolyLine(mergeType=IMPRINT, meshable=False, points=((model.rootAssembly.referencePoints[88], model.rootAssembly.instances['centre_timber'].vertices[0]), ))

    model.rootAssembly.features.changeKey(fromName='Wire-25', toName='Wire-25')

    model.rootAssembly.Set(edges=model.rootAssembly.edges.getSequenceFromMask(('[#1 ]', ), ), name='Wire-25-Set-1')

    model.rootAssembly.SectionAssignment(region=model.rootAssembly.sets['Wire-25-Set-1'], sectionName='Spring-0-Section')

    model.rootAssembly.sectionAssignments[24].getSet()

    model.rootAssembly.ConnectorOrientation(
        localCsys1=model.rootAssembly.datums[89], 
        region=model.rootAssembly.allSets['Wire-25-Set-1']
    )"""

    assembly.regenerate()






    # model.HistoryOutputRequest(
    #     name='H-Output-Springs',
    #     createStepName='Step-1',
    #     region=assembly.sets['Spring-C-Set'],  # or a set covering all springs
    #     variables=('CTF1', 'CTF2', 'CU1', 'CU2')
    # )

    # model.FieldOutputRequest(createStepName='Step-1', name=
    #     'F-Output-3', position=INTEGRATION_POINTS, rebar=EXCLUDE, region=
    #     mdb.models['Model-1'].rootAssembly.sets['Spring-C-Set'], sectionPoints=
    #     DEFAULT, variables=('CTF', 'CEF', 'CU', 'CUE', 'CUP'))

    # model.rootAssembly.Set(edges=
    #     mdb.models['Model-1'].rootAssembly.instances['left_timber'].edges.getSequenceFromMask(
    #     mask=('[#1 ]', ), )+\
    #     mdb.models['Model-1'].rootAssembly.instances['right_timber'].edges.getSequenceFromMask(
    #     mask=('[#1 ]', ), ), name='Set-4', vertices=
    #     mdb.models['Model-1'].rootAssembly.instances['left_timber'].vertices.getSequenceFromMask(
    #     mask=('[#1 #4000000 ]', ), )+\
    #     mdb.models['Model-1'].rootAssembly.instances['right_timber'].vertices.getSequenceFromMask(
    #     mask=('[#1 #4000000 ]', ), )+\
    #     mdb.models['Model-1'].rootAssembly.instances['centre_timber'].vertices.getSequenceFromMask(
    #     mask=('[#1 #4000000 ]', ), )
    # )
    # model.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    # distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    # 'BC-1', region=model.rootAssembly.sets['Set-4'], u1=0.0, 
    # u2=0.0, ur3=UNSET)

    # mdb.models['Model-1'].rootAssembly.Set(name='Set-5', vertices=
    #     mdb.models['Model-1'].rootAssembly.instances['steel_dowel'].vertices.getSequenceFromMask(
    #     ('[#0 #4000000 ]', ), ))
    # mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    #     distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    #     'BC-2', region=mdb.models['Model-1'].rootAssembly.sets['Set-5'], u1=UNSET, 
    #     u2=1.0, ur3=UNSET)





if __name__ == '__main__':
    #main_double_shear()
    spring_test()








# grave
"""def create_nonlinear_connector(
        name,
        instance_1,
        vertex_1_index,
        instance_2,
        vertex_2_index,
        force_displacement_table):
    
    # Create a nonlinear connector between two existing vertices.
    # force_displacement_table : tuple
    #     Example:
    #         (
    #             (-5000.0 N, -5.0 mm),
    #             (-1000.0 N, -1.0 mm),
    #             (    0.0,  0.0),
    #             ( 1000.0,  1.0),
    #             ( 5000.0,  5.0),
    #         )

    # ----------------------------------------------------------
    # Names
    # ----------------------------------------------------------

    set_name = name + '-Set'
    section_name = name + '-Section'

    # ----------------------------------------------------------
    # Get the two vertices
    # ----------------------------------------------------------

    vertex_1 = assembly.instances[
        instance_1
    ].vertices[
        vertex_1_index
    ]

    vertex_2 = assembly.instances[
        instance_2
    ].vertices[
        vertex_2_index
    ]

    # Get coordinates
    p1 = vertex_1.pointOn[0]
    p2 = vertex_2.pointOn[0]

    # ----------------------------------------------------------
    # Calculate midpoint of connector
    # ----------------------------------------------------------

    midpoint = (
        0.5 * (p1[0] + p2[0]),
        0.5 * (p1[1] + p2[1]),
        0.5 * (p1[2] + p2[2])
    )

    # ----------------------------------------------------------
    # Create connector wire
    # ----------------------------------------------------------

    assembly.WirePolyLine(
        points=((vertex_1, vertex_2),),
        mergeType=SEPARATE,
        meshable=OFF
    )

    assembly.regenerate()

    # ----------------------------------------------------------
    # Get the newly-created wire
    # ----------------------------------------------------------

    connector_edge = assembly.edges.findAt(
        (midpoint,)
    )

    # ----------------------------------------------------------
    # Create connector set
    # ----------------------------------------------------------

    assembly.Set(
        name=set_name,
        edges=connector_edge
    )

    # ----------------------------------------------------------
    # Create nonlinear elastic behavior
    #
    # Component 1 = global X
    # Component 2 = global Y
    # ----------------------------------------------------------

    connector_elasticity = ConnectorElasticity(
        behavior=RIGID,
        coupling=UNCOUPLED,
        components=(2,),
        table=force_displacement_table
    )

    # ----------------------------------------------------------
    # Create connector section
    # ----------------------------------------------------------

    model.ConnectorSection(
        name=section_name,
        translationalType=CARTESIAN,
        behaviorOptions=(connector_elasticity,),
        rotationalType=NONE
    )

    # ----------------------------------------------------------
    # Assign section
    # ----------------------------------------------------------

    assembly.SectionAssignment(
        region=assembly.sets[set_name],
        sectionName=section_name
    )

    assembly.regenerate()"""



"""
tests
"""



