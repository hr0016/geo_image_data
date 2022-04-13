import argparse
import rasterio 
import numpy as np
import matplotlib.pyplot as plt

def generator(rows, columns, min_terrain_height, max_terrain_height):
    r,c=np.arange(rows), np.arange(columns)
    r,c=np.meshgrid(r,c)
    r,c=np.transpose(r), np.transpose(c)
    c_=np.array([list(reversed(k)) for k in c])
    r_=np.array(list(reversed(r)))
    terrain=np.random.rand(1)*r+np.random.rand(1)*c+np.random.rand(1)*r_+np.random.rand(1)*c_
    terrain=(terrain-terrain.min())/(terrain.max()-terrain.min())*(max_terrain_height-min_terrain_height)+min_terrain_height
    return terrain
    
if __name__ == '__main__':    
    # create a parser object
    parser = argparse.ArgumentParser(description = "")
    # "free space path loss model" "random terrain data generator"
    parser.add_argument("Model_Name", type = str)
    parser.add_argument("input_file", type = str)
    parser.add_argument("resolution", type = int)
    parser.add_argument("rows", type = int)
    parser.add_argument("columns", type = int)
    parser.add_argument("min_terrain_height", type = int)
    parser.add_argument("max_terrain_height", type = int)
    
    args = parser.parse_args()
    Model_Name=args.Model_Name
    input_file=args.input_file
    resolution = args.resolution
    rows = args.rows
    columns = args.columns
    min_terrain_height = args.min_terrain_height
    max_terrain_height = args.max_terrain_height
   
    input_data = rasterio.open(input_file)
    input_data = input_data.read(1)
    
    terrain = generator(rows, columns, min_terrain_height, max_terrain_height)
    plt.imshow(terrain)