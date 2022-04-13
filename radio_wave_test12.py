import argparse
import rasterio 
import numpy as np

def compute_receiver_h(resolution, max_height):
    receiver_h = np.arange(resolution, max_height + resolution, resolution)
    return receiver_h

def compte_FSPL(input_data, resolution, receiver_h, transmitter_x, transmitter_y, transmitter_h, transmitter_freq):
    tx_terrain_h = input_data[transmitter_y, transmitter_x]
    tx_actual_h = tx_terrain_h + transmitter_h
    rows, columns = len(input_data), input_data.shape[1]
    FSPL = np.empty((rows, columns, len(receiver_h)))
    for h in range(len(receiver_h)):
        rx_actual_h = input_data + receiver_h[h]
        for y in range(rows):
            for x in range(columns):
                d = np.sqrt(resolution**2*((x-transmitter_x)**2+(y-transmitter_y)**2)+(rx_actual_h[y][x]-tx_actual_h)**2)
                FSPL[y][x][h] = 20*(np.log10(d/1000) + np.log10(transmitter_freq/1E+6)) + 32.45
    return FSPL

def get_terrain_height_profile(input_data, resolution, receiver_x, receiver_y, transmitter_x, transmitter_y):
    link = line_to_coords(transmitter_y, transmitter_x, receiver_y, receiver_x)
    terrain_height_profile = [input_data[l] for l in link]
    return link, terrain_height_profile

def identify_knife_edge(link, terrain_height_profile, transmitter_h, receiver_h_actual):
    link_, link_h=[], []
    dif_h = terrain_height_profile[0] - receiver_h_actual
    for l in range(1, len(link)-1):
        ref_h = dif_h*(len(link)-l)/len(link)
        if terrain_height_profile[l] - terrain_height_profile[l-1] < 0 and terrain_height_profile[l] - terrain_height_profile[l+1] >= 0 and terrain_height_profile[l] > ref_h:
            link_.append(link[l])  
            link_h.append(terrain_height_profile[l])
    return link_, link_h
    
def compute_KED(input_data, resolution, max_height, transmitter_h, transmitter_freq, transmitter_x, transmitter_y, receiver_x, receiver_y, link_, link_h):
    transmitter_h_actual = input_data[transmitter_y, transmitter_x] + transmitter_h
    receiver_h_actual = input_data[receiver_y, receiver_x] + max_height
    d = np.sqrt(resolution**2*((receiver_x-transmitter_x)**2+(receiver_y-transmitter_y)**2)+(receiver_h_actual-transmitter_h_actual)**2)
    FSPL = 20*(np.log10(d/1000) + np.log10(transmitter_freq/1E+6)) + 32.45
    L, ind= 0, list(np.arange(len(link_)))
    while len(ind) > 1:
        vlist=[]      
        for l in range(len(ind)):
            d1=np.sqrt(resolution**2*((link_[l][0]-transmitter_x)**2+(link_[l][1]-transmitter_y)**2)+(link_h[l]-transmitter_h_actual)**2)
            d2=np.sqrt(resolution**2*((link_[l][0]-receiver_x)**2+(link_[l][1]-receiver_y)**2)+link_h[l]**2+(link_h[l]-receiver_h_actual)**2)
            wavelength=3E+8/transmitter_freq
            vlist.append(link_h[l]*np.sqrt(2*(d1+d2)/(wavelength*d1*d2)))
        ind=np.argsort(vlist)[::-1]+1
        ind=list(np.arange(ind[0]))
        if vlist[ind[0]] >= -0.7:
            J=6.9+20*np.log10(np.sqrt((vlist[ind[0]]-0.1)**2+1)+vlist[ind[0]]-0.1)
        else:
            J=0
        L+=J
    KED=FSPL+L
    return KED
    
def line_to_coords(x,y,x2,y2):
    steep = 0
    coords = []
    dx = abs(x2 - x)
    if (x2 - x) > 0: sx = 1
    else: sx = -1
    dy = abs(y2 - y)
    if (y2 - y) > 0: sy = 1
    else: sy = -1
    if dy > dx:
        steep = 1
        x,y = y,x
        dx,dy = dy,dx
        sx,sy = sy,sx
    d = (2 * dy) - dx
    for i in range(0,dx):
        if steep: coords.append((y,x))
        else: coords.append((x,y))
        while d >= 0:
            y = y + sy
            d = d - (2 * dx)
        x = x + sx
        d = d + (2 * dy)
    coords.append((x2,y2))
    return coords   
    
if __name__ == '__main__':    
    # create a parser object
    parser = argparse.ArgumentParser(description = "")
    # "free space path loss model" "random terrain data generator"
    parser.add_argument("Model_Name", type = str)
    parser.add_argument("input_file", type = str)
    parser.add_argument("resolution", type = int)
    parser.add_argument("rows", type = int)
    parser.add_argument("columns", type = int)
    parser.add_argument("max_height", type = int)
    parser.add_argument("transmitter_x", type = int)
    parser.add_argument("transmitter_y", type = int)
    parser.add_argument("transmitter_h", type = int)
    parser.add_argument("transmitter_freq", type = int)
    parser.add_argument("output_file", type = str)
    
    args = parser.parse_args()
    Model_Name=args.Model_Name
    input_file=args.input_file
    resolution = args.resolution
    rows = args.rows
    columns = args.columns
    max_height = args.max_height
    transmitter_x = args.transmitter_x
    transmitter_y = args.transmitter_y    
    transmitter_h = args.transmitter_h
    transmitter_freq = args.transmitter_freq
    output_file = args.output_file
   
    input_data = rasterio.open(input_file)
    input_data = input_data.read(1)  
    
    receiver_h=compute_receiver_h(resolution, max_height)
    FSPL = compte_FSPL(input_data, resolution, receiver_h, transmitter_x, transmitter_y, transmitter_h, transmitter_freq)
        
    transmitter_h = 0
    receiver_x, receiver_y = columns-1, rows-1
    receiver_h_actual = input_data[receiver_y, receiver_x]+max_height
    link, terrain_height_profile = get_terrain_height_profile(input_data, resolution, receiver_x, receiver_y, transmitter_x, transmitter_y)
    link_, link_h = identify_knife_edge(link, terrain_height_profile, transmitter_h, receiver_h_actual)
    KED = compute_KED(input_data, resolution, max_height, transmitter_h, transmitter_freq, transmitter_x, transmitter_y, receiver_x, receiver_y, link_, link_h)

    print("receiver height above terrain:","\n", receiver_h)
    print("FSPL model free space loss table map for any receiver location and height:","\n", FSPL)
    print("KED model diffraction loss:","\n", KED)
    