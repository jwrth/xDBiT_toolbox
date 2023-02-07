from matplotlib.colors import ListedColormap, rgb2hex

class color_palettes:
    '''
    Class containing a collection of custom color palettes.
    '''
    def __init__(self):
        # palette for colorblind people. From: https://gist.github.com/thriveth/8560036
        self.colorblind = ListedColormap(
            ['#377eb8', '#ff7f00', '#4daf4a',
             '#f781bf', '#dede00', '#a65628', 
             '#984ea3', '#999999', '#e41a1c'], name="colorblind")
        
        # palette from Caro. Optimized for colorblind people.
        self.caro = ListedColormap(['#3288BD','#440055', '#D35F5F', '#A02C2C','#225500', '#66C2A5', '#447C69'], name="caro")
        
        # from https://thenode.biologists.com/data-visualization-with-flying-colors/research/
        self.okabe_ito = ListedColormap(["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"], name="okabe_ito")
        self.tol_bright = ListedColormap(["#EE6677", "#228833", "#4477AA", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"], name="tol_bright")
        self.tol_muted = ListedColormap(["#88CCEE", "#44AA99", "#117733", "#332288", "#DDCC77", "#999933", "#CC6677", "#882255", "#AA4499", "#DDDDDD"], name="tol_muted")
        self.tol_light = ListedColormap(["#BBCC33", "#AAAA00", "#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#DDDDDD"], name="tol_light")
        
        
def cmap2hex(cmap):
    '''
    Generate list of hex-coded colors from cmap.
    '''
    hexlist = [rgb2hex(cmap(i)) for i in range(cmap.N)]
    return hexlist