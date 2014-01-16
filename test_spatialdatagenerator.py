from SpatialDataGenerator import SpatialDataGenerator

def run():
    SDG = SpatialDataGenerator("data/delhi.cfg")
    SDG.generate_spatial_data()
    SDG.plot_map()
    
if __name__ == '__main__':
    run()