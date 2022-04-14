#%% import necessary packages; note rasterio is required to read image from the terrain.bil file
import argparse
import rasterio 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#%% a function to generate a random terrain image of fixed numbers of rows and columns, and minmum and maximum height. 

def generator(rows, columns, min_terrain_height, max_terrain_height):
    r,c=np.arange(rows), np.arange(columns)
    r,c=np.meshgrid(r,c)
    r,c=np.transpose(r), np.transpose(c)
    c_=np.array([list(reversed(k)) for k in c])
    r_=np.array(list(reversed(r)))
    terrain=np.random.rand(1)*r+np.random.rand(1)*c+np.random.rand(1)*r_+np.random.rand(1)*c_
    terrain=(terrain-terrain.min())/(terrain.max()-terrain.min())*(max_terrain_height-min_terrain_height)+min_terrain_height
    return terrain

# %% the main function, parsed using the Python argparse module, to execute the main function using 
# the IDE's CLI, please follow "Run" -> "Configuration per file" -> tick box "command line options" and 
# input your syntax line for all positional arguments in order    
if __name__ == '__main__':    
    # create a parser object
    parser = argparse.ArgumentParser(description = "terrain image generation mode")
    # "free space path loss model" "random terrain data generator"
    parser.add_argument("Model_Name", type = str, help="path loss model name")
    parser.add_argument("input_file", type = str, help="the image file name")
    parser.add_argument("resolution", type = int, help="row and column resolution of the image (m)")
    parser.add_argument("rows", type = int, help="number of rows")
    parser.add_argument("columns", type = int, help="number of columns")
    parser.add_argument("min_terrain_height", type = int, help="minimum terrain height (m)")
    parser.add_argument("max_terrain_height", type = int, help="maximum terrain height (m)")
    parser.add_argument("save_file", type = str, help="save generated filename")
    
    args = parser.parse_args()
    Model_Name=args.Model_Name
    input_file=args.input_file
    resolution = args.resolution
    rows = args.rows
    columns = args.columns
    min_terrain_height = args.min_terrain_height
    max_terrain_height = args.max_terrain_height
    save_file = args.save_file
   
    input_data = rasterio.open(input_file)
    input_data = input_data.read(1)
    
    terrain = generator(rows, columns, min_terrain_height, max_terrain_height)
    plt.imshow(terrain, cmap='gray')
    # save the generated file so it can be tested with path loss models
    matplotlib.image.imsave('{}.tiff'.format(save_file), terrain)