import os
from pathlib import Path
import pandas as pd
import cv2


def save_regprops(regprops, savedir, save_images=False, image_dir=None):
    '''
    Function to save the region properties of an image.
    '''
    save_images = False

    # create feature directory if necessary
    Path(savedir).mkdir(parents=True, exist_ok=True)

    if save_images:
        if image_dir is None:
            image_dir = os.path.join(savedir, "images")
        Path(image_dir).mkdir(parents=True, exist_ok=True)

    for idx, props in regprops.items():
        print("Saving {}...".format(idx))
        # save binary regionprops
        savefile = os.path.join(savedir, "featureprops-{}-binary.csv".format(idx))
        props['binary'].to_csv(savefile, index=False)
        
        # save intensity based regionprops
        labels = props['intensity'].keys()
        for label in labels:
            print("\t{}".format(label))
            # extract dataframe
            data = props['intensity'][label]
            nfill = len(str(len(data))) # needed for leading zeros later
            
            if save_images:
                # save images
                img_dir = os.path.join(image_dir, "{}-{}".format(idx, label))
                Path(img_dir).mkdir(parents=True, exist_ok=True)
                for i, img in data['image_intensity'].iteritems():
                    ii = str(i).zfill(nfill)
                    image_file = os.path.join(img_dir, "nucleus-{}-{}-{}.tif".format(idx, label, ii))
                    cv2.imwrite(image_file, img)
                
            # save dataframe without image column
            savefile = os.path.join(savedir, "featureprops-{}-{}.csv".format(idx, label))
            data.drop('image_intensity', axis=1).to_csv(savefile, index=False)

def load_regprops(feature_files):
    regprops = {}

    for file in feature_files:
        # get informations from file name
        idx = os.path.basename(file).split("-")[1]
        label = os.path.basename(file).split("-")[2].rstrip(".csv")
        
        # read dataframe
        df = pd.read_csv(file)
        
        # check if this index already exists in dictionary
        if idx not in regprops:
            regprops[idx] = {}
            regprops[idx]['intensity'] = {}
        
        if label == 'binary':
            regprops[idx][label] = df
        else:
            regprops[idx]['intensity'][label] = df

    return regprops
