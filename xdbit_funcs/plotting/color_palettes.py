class color_palettes:
    '''
    Class containing a collection of custom color palettes.
    '''
    def __init__(self):
        # palette for colorblind people. From: https://gist.github.com/thriveth/8560036
        self.colorblind = ['#377eb8', '#ff7f00', '#4daf4a',
                           '#f781bf', '#a65628', '#984ea3',
                           '#999999', '#e41a1c', '#dede00']
        
        # palette from Caro. Optimized for colorblind people.
        self.caro = ['#3288BD','#440055', '#D35F5F', '#A02C2C','#225500', '#66C2A5', '#447C69']