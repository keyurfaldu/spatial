import sys
import ConfigParser
from scipy import ndimage
from scipy import misc
import matplotlib.pyplot as plt
import os
import pyproj
import math
import rpy2.robjects as ro
import rpy2.robjects as robjects
from optparse import OptionParser


robjects.r('''source('spatial_plotter.R')''')

def load_image(filename):
    img = misc.imread(filename)
    return img

def is_adjacent_to_region(img, px, py, r, g, b):
    for (x, y) in [(px+1, py), (px+1, py+1), (px+1, py-1), (px, py+1), (px, py-1), (px-1, py), (px-1, py-1), (px-1, py+1)]:
        if img[x,y,0] == r and img[x,y,1] == g and img[x,y,2] == b:
            return True
    return False

def is_border(img, x, y):
    cutoff = 35 #Even if it is dark gray, consider it as black
    if img[x,y,0] <= cutoff and img[x, y, 1] <= cutoff and img[x, y, 2] <= cutoff:
        return True
    return False

def get_valid_adjacents(img, px, py, r, g, b):
    valid_adjacents = []
    for (x, y) in [(px+1, py), (px+1, py+1), (px+1, py-1), (px, py+1), (px, py-1), (px-1, py), (px-1, py-1), (px-1, py+1)]:
        if is_border(img, x, y) and is_adjacent_to_region(img, x, y, r, g, b):
            img[x,y][:3] = [255, 1, 1]    #Mark it red.
            valid_adjacents.append((x, y))
    return valid_adjacents

def floodFill(img,x,y,r,g,b):
    toFill = []
    toFill.append((x,y))
    nimg = img.copy()
    sp = (-1, -1)
    while toFill:
        (x,y) = toFill.pop()
        a,b,c = nimg[x,y][:3]
        if not (a,b,c) >= (250, 250, 250):
            if a == 0 and sp[0] < x:
                sp = (x,y)
            continue
        nimg[x,y][:3] = [r,g,b]
        toFill.append((x-1,y))
        toFill.append((x+1,y))
        toFill.append((x,y-1))
        toFill.append((x,y+1))
    return nimg, sp
    
def get_polygon(img, p, region_color):
    polygon = []
    r, g, b = region_color
    x, y = p
    nimg, sp = floodFill(img,x,y,r,g,b)
    x, y = sp
    adjacents = get_valid_adjacents(nimg, x, y, r, g, b)
    backtrack_adjacents = []
    while adjacents:
        np = adjacents.pop()
        polygon.append(np)
        x, y = np
 
        new_adjacents = get_valid_adjacents(nimg, x, y, r, g, b)
        while not new_adjacents and not adjacents and backtrack_adjacents and False:
            np = backtrack_adjacents.pop()
            x, y = np
            new_adjacents = get_valid_adjacents(nimg, x, y, r, g, b)            
        if new_adjacents:
            adjacents = new_adjacents
            backtrack_adjacents.extend(new_adjacents)
    return nimg, polygon

class BestEffortConversion:
    def __init__(self):
        #mercator 
        self.m_pixel = pyproj.Proj(init='EPSG:3857') 
        #lat long coordinates
        self.m_latlong = pyproj.Proj(init='EPSG:4326')
        #superimposed points
        self.points = []

    def register_superimposed_point(self, xy, lglt):
        pxpy = pyproj.transform(self.m_latlong, self.m_pixel, lglt[0], lglt[1])
        self.points.append((xy, lglt, pxpy))
        
    def map_point(self, x, y):
        #find the nearest point
        min_d = -1
        nearest_point = None
        i = -1
        for point in self.points:
            i += 1
            xy, lglt, pxpy = point
            px, py = pxpy
            d = (px - x)*(px - x) + (py - y)*(py - y)
            if min_d < 0:
                min_d = d
                nearest_point = i
            if d < min_d:
                min_d = d
                nearest_point = i
        
        extrapolated_points = []
        max_d = 0

        for p1 in self.points: 
            for p2 in self.points:
                if p1 == p2: continue
                x1, y1 =  p1[0]
                #lg1, lt1 = p1[1]
                lg1, lt1 = p1[2]
                x2, y2 =  p2[0]
                #lg2, lt2 = p2[1]
                lg2, lt2 = p2[2]
                        
                if not(abs(x1 - x2) > 200 and abs(y1 - y2) > 200):
                    continue
                
                d1 = self.ecludian(x, y, x1, y1)
                d2 = self.ecludian(x, y, x2, y2)
            	d = min(d1, d2)
            	
                lg = lg1 + (x - x1)*(lg2 - lg1)/(x2 - x1)
                lt = lt1 + (y - y1)*(lt2 - lt1)/(y2 - y1)
                if d > max_d:
                    max_d = d
                extrapolated_points.append((lg, lt, d))

        if max_d == 0:
            return  self.points[nearest_point][1]   
        coeff = 100.0/max_d
        wlg = 0
        wlt = 0
        denominator = 0
        for (lg, lt, d) in extrapolated_points:
            d = d+1
            scale = 1/((coeff*d)*(coeff*d))
            denominator += scale
            wlg += lg*scale
            wlt += lt*scale
        wlg = wlg/denominator
        wlt = wlt/denominator
        pt = pyproj.transform(self.m_pixel, self.m_latlong, wlg, wlt)

        return pt
        
    def ecludian(self, x1, y1, x2, y2):
        return math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))


class SpatialDataGenerator:
    def __init__(self, config_file):
        #load config
        self.config = ConfigParser.ConfigParser()
        self.config.read(config_file)
        
        #Load main section
        self.map = eval(self.config.get("MAIN","Map"))
        self.image_file = eval(self.config.get("MAIN","Mask_png"))
        self.num_regions = eval(self.config.get("MAIN","Regions"))
        self.out_file = eval(self.config.get("MAIN","Output_spatial_file"))
        self.num_superimposed_points = eval(self.config.get("MAIN","Num_superimposed_points"))
      
        self.xlab = "Longitude"
        self.ylab = "Lattitude"
        self.title = eval(self.config.get("MAIN","Title"))

        #Load Image
        self.img = load_image(self.image_file)
        
        self.BEC = BestEffortConversion()
        for i in xrange(1, self.num_superimposed_points +1):
            xy_lglt = eval(self.config.get("SUPERIMPOSE", "P%s"%i))
            self.BEC.register_superimposed_point(xy_lglt[0], xy_lglt[1])
            
        
    def generate_spatial_data(self):
        out_data = []
        try:
            for i in xrange(1, self.num_regions+1):
                region = "R%s"%(i)
                region_code = "%s-%s"%(self.map, region)
                region_name = self.config.get(region, "Name")
                            
                XY = eval(self.config.get(region, "XY"))
              
                if XY == (0, 0):
                    continue
                nimg, polygon = get_polygon(self.img, (XY[1], XY[0]), (100, 100, 255))
                
                count = 0
                
                for (x, y) in polygon:
                    count += 1
                    lt, lg = self.BEC.map_point(y,x)
                    out_data.append(",".join(map(str,[lt, lg, count,"FALSE",'"1"',"%s.1"%region_code, region_code, region_name])))  
        except Exception, e:
            print "Exception: %s, Region: %s"%(e, region)
            
        o = open(self.out_file, 'w')
        for l in out_data:
            o.write(l)
            o.write("\n")
        o.close()
        
    def plot_map(self):
        image_file = "%s.png"%self.out_file.rsplit('.')[0]
        r_plot_map = robjects.globalenv["plot.map"]
        r_plot_map(self.out_file, image_file, self.xlab, self.ylab, self.title)
 
def main(options, args):
    if not options.config_file:
        print "Error: Please specify config file"
        sys.exit(0)
    SDG = SpatialDataGenerator(options.config_file)
    SDG.generate_spatial_data()
    SDG.plot_map()
     
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-c", "--config-file", help="config file, example specified in data/delhi.cfg")
    (options, args) = parser.parse_args()
    main(options,args)     
            

