import chimera, os, math, numpy, time
from chimera import runCommand as rc
from chimera import replyobj, Point, openModels, Molecule # for emitting status messages
from AddCharge import initiateAddCharges
import SurfaceColor as sc
from chimera import selection,MSMSModel
from numpy.linalg import norm

'''
This script can be opened in the Chimera GUI under File -> Open
'''

# gather the names of .pdb files in the folder
file_names = [fn for fn in os.listdir(".") if fn.endswith(".sdf")]
file_names = sorted(file_names)
final_array = []
hm_list = []
# loop through the files, opening, processing, and closing each in turn

def area(v1, v2, v3):
        return 0.5 * norm(numpy.cross(v2-v1, v3-v1))
        
        
for fn in file_names:
        
    #Open the file and extract the first line which contains the name of the molecule
    with open(fn) as myfile:
        head = [next(myfile) for x in xrange(1)]
        
        
#------------------------------------------------------------------
# Processing the molecule using Chimera
#------------------------------------------------------------------
    
    #Open the file in Chimera
    rc("open " + fn)
    
    #Add Gasteiger charges
    initiateAddCharges(method='gasteiger', nogui=True)
    rc('surface vertexDensity 2') # Add surface to molecule
    rc('coulombic  -10 red 0 white 10 blue') # Colr surface by coulombic charge
    rc('select') # Select the molecule from the model panel, also selects the surface for processing
    
    

    
#------------------------------------------------------------------
# Surface processing
#------------------------------------------------------------------
    #Surface processing section
    slist = [m for m in selection.currentGraphs() if isinstance(m, MSMSModel)] # initiates selected surface to object
    vertices = [] # empty list to populate with vertices
    s = slist[0] # takes the surface from the surface object
    
    p = s.surfacePieces[0] # Function of MSMSModel
    va, ta = p.geometry #Extract the list of vertices and triangles 
    rgba = p.vertexColors    
    '''
    coords = []
    for i, v in enumerate(va):
        x, y, z = v #Get cartesian coordinates of vertex
        point = (x, y, z) #Convert coords to tuple.....durh, probably could have just used the v vertex object here
        coords.append(point) #Add the distance to the coords list
    '''    
        
#------------------------------------------------------------------
# ELECTROSTATIC CALCULATION
#------------------------------------------------------------------

    '''Note: the color scheme in Chimera is set up so that R is 1 if there is negative charge and 
    B is 1 if there is any positive charge. Therefore, got to see which one is set to one to determine if we
    have a positive or negatively cahtged vertex. The second step is to scale the color to the maximum estst energy, 10 kcal/molecule
    

    rgba = p.vertexColors # Extract list of colors of each vertex into an array [Note: a is for transparency, usually set to 1]
    estats = []
    total_estat = 0
    for j in rgba:
        max_rgb = 1.0
        if j[2] == max_rgb:
            estat = abs((10*(1-j[0])))
            estats.append(estat)
            total_estat = total_estat + estat
        else:
            estat = abs((-10*(1-j[2])))
            estats.append(estat)
            total_estat = total_estat + estat

    ave_estat = total_estat/len(rgba)
    '''
    
#------------------------------------------------------------------
# Calculate area of all triangles (Equation 4 from Reisser paper)
#------------------------------------------------------------------

    '''
    From equation 4 of the paper, we need to take each triangle (ta) and calculate its area according to the function below.
    This is accomplished by first enumerating va and ta from p.geometry. Then for the three vertives in each ta, get the corresponding va and feed these 3 arrays into 
    the function.
    A cumulative total of areas is calculated
    '''
    
    cum_area = 0
    for n in ta:
        vert1 = n[0]
        vert2 = n[1]
        vert3 = n[2]
        triangle_area = area(va[vert1], va[vert2], va[vert3])
        cum_area = triangle_area + cum_area
    print ta[0]
    print "TA[0] - cum area"
    print "CUM_AREA"
    print cum_area
#------------------------------------------------------------------
# Calculate area-scaled vector of each vertex (Equation 3 from Reisser paper)
#------------------------------------------------------------------
    vecs = []
    max_rgb = 1.0
    total_estat = 0
    
    for n in ta:
        triangle_estat_sum = 0
        estats = []
        vert1 = n[0]
        vert2 = n[1]
        vert3 = n[2]
        triangle_area = area(va[vert1], va[vert2], va[vert3])
        scaled_area = triangle_area/(3*cum_area)
        r1 = numpy.multiply(va[vert1], scaled_area)
        r2 = numpy.multiply(va[vert2], scaled_area)
        r3 = numpy.multiply(va[vert3], scaled_area)
        vecs.append(r1)
        vecs.append(r2)
        vecs.append(r3)
        
        point1 = rgba[vert1]
        if point1[2] == max_rgb:
            estat1 = abs((10*(1-point1[0])))
            estat_sum1 = estat1*(triangle_area/3)
            estats.append(estat_sum1)
            
        else:
            estat1 = abs((-10*(1-point1[2])))
            estat_sum1 = estat1*(triangle_area/3)
            estats.append(estat_sum1)
           
            
        point2 = rgba[vert2]
        if point1[2] == max_rgb:
            estat2 = abs((10*(1-point2[0])))
            estat_sum2 = estat2*(triangle_area/3)
            estats.append(estat_sum2)
            
        else:
            estat2 = abs((-10*(1-point2[2])))
            estat_sum2 = estat2*(triangle_area/3)
            estats.append(estat_sum2)
            
            
        point3 = rgba[vert3]
        if point3[2] == max_rgb:
            estat3 = abs((10*(1-point3[0])))
            estat_sum3 = estat3*(triangle_area/3)
            estats.append(estat_sum3)
           
        else:
            estat3 = abs((-10*(1-point3[2])))
            estat_sum3 = estat3*(triangle_area/3)
            estats.append(estat_sum3)
                    
        triangle_estat_sum = sum(estats)
        total_estat = total_estat + triangle_estat_sum
        
    
    average_estat = total_estat/cum_area
    print "Ave Estat"
    print average_estat
#Equation 1
    hm = numpy.empty([1,3])
    for n in ta:
        triangle_estat_sum = 0
        estats = []
        vert1 = n[0]
        vert2 = n[1]
        vert3 = n[2]
        triangle_area = area(va[vert1], va[vert2], va[vert3])
        scaled_area = triangle_area/(3*cum_area)
        r1 = numpy.multiply(va[vert1], scaled_area)
        r2 = numpy.multiply(va[vert2], scaled_area)
        r3 = numpy.multiply(va[vert3], scaled_area)
        vecs.append(r1)
        vecs.append(r2)
        vecs.append(r3)
        
        point1 = rgba[vert1]
        if point1[2] == max_rgb:
            estat1 = abs((10*(1-point1[0])))
            estat_sum1 = estat1*(triangle_area/3)
            vect1 = numpy.multiply(average_estat-estat_sum1, r1)
            
        else:
            estat1 = abs((-10*(1-point1[2])))
            estat_sum1 = estat1*(triangle_area/3)
            vect1 = numpy.multiply(average_estat-estat_sum1, r1)
           
            
        point2 = rgba[vert2]
        if point1[2] == max_rgb:
            estat2 = abs((10*(1-point2[0])))
            estat_sum2 = estat2*(triangle_area/3)
            vect2 = numpy.multiply(average_estat-estat_sum2, r2)
            
        else:
            estat2 = abs((-10*(1-point2[2])))
            estat_sum2 = estat2*(triangle_area/3)
            vect2 = numpy.multiply(average_estat-estat_sum2, r2)
            
            
        point3 = rgba[vert3]
        if point3[2] == max_rgb:
            estat3 = abs((10*(1-point3[0])))
            estat_sum3 = estat3*(triangle_area/3)
            vect3 = numpy.multiply(average_estat-estat_sum3, r3)
           
        else:
            estat3 = abs((-10*(1-point3[2])))
            estat_sum3 = estat3*(triangle_area/3)
            vect3 = numpy.multiply(average_estat-estat_sum3, r3)

        triangle_moment = numpy.add(vect1, vect2, vect3)
        hm = hm+triangle_moment
        

    hm_final = str(numpy.linalg.norm(hm))
    hm_couple = (head, hm_final)
        
        
    print "HM"
    print hm
    print hm_final
    print 'len TA'
    print len(ta)
    hm_list.append(str(hm_couple))
    
    


#------------------------------------------------------------------
# Convert electrostatic and distance values to hydrophobic moments
#------------------------------------------------------------------
    '''
    all_values = zip(coords, estats)
    
    hm = numpy.empty([1,3])
    for k in all_values:
        o = (ave_estat-k[1])*numpy.array(k[0])
        hm = hm + o
        

    print hm

    
    '''
    rc("close all")
    '''
    hm_final = [numpy.linalg.norm(hm)]

    final_values = zip(head, hm_final)
    final_array.append(final_values)
    
    print final_array'''

hm_string = "".join(hm_list)
text_file = open('hm_lis-2.txt', 'w')
text_file.write(str(hm_list))
text_file.close()

