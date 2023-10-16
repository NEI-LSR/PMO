''' Setup '''
import pandas as pd
import numpy as np
from psychopy import visual, monitors, event
import pickle

### PARAMETERS -- change this!
isphotometer = 1
ismeg = 1
isfullscreen = 1
isinitmeasure = 1

refreshrate = 60
view_dist = 75
screen_width = 42

def run_calibration(sdict,repeats):

    ## HELPER FUNCTIONS
    def setup_window(background_color = 0.5, color_space = 'rgb1', isfullscreen = 0, ismeg = ismeg):
        if ismeg:
            #mon = monitors.Monitor('MEG_20221028')
            win = visual.Window(color=background_color,colorSpace=color_space,units='pix',checkTiming=True,fullscr=isfullscreen)
        else: 
            win = visual.Window(color=background_color,colorSpace=color_space,units='pix',checkTiming=True,fullscr=isfullscreen)

        return win 

    def deg_to_pix(dva,win,screen_width,view_dist):
        size_cm = view_dist*np.tan(np.deg2rad(dva/2))*2
        pix_per_cm = win.size[0]/screen_width
        size_pix = size_cm*pix_per_cm
        return int(np.round(size_pix))

    def make_circle(win, color = [1,0,0], stim_size = 200, colorSpace = 'rgb1'):
        stim = visual.Circle(win, radius=0.5, units='pix', fillColor=color, size=stim_size, colorSpace=colorSpace)
        return stim

    def get_photometer():
        from psychopy.hardware.pr import PR655
        photom = PR655('COM7')
        return photom

    def draw_stim_getspectra(win,stim,photometer):
        stim.draw()
        win.flip()    
        lum = photometer.getLum()
        nm, power = photometer.getLastSpectrum(parse=True)
        
        return nm, power, lum


    # Setup
    if isphotometer:
        photom = get_photometer()

    win = setup_window(background_color = [0.5,0.5,0.5], color_space = 'rgb1', isfullscreen = isfullscreen, ismeg = ismeg)
    stim_size_deg = 15
    stim_size = deg_to_pix(stim_size_deg,win,screen_width,view_dist)

    # make stimuli in the right colour
    for i in range(len(all_colours)):
        col = (all_colours.to_numpy()/255)[i]
        circle = make_circle(win, color = col, stim_size = stim_size, colorSpace = 'rgb1')

        # loop over the colours and if a spectrophotometer is connected, save the spectra
        if isphotometer:
            nm, power, lum= draw_stim_getspectra(win,circle,photometer=photom)
            print(lum)
            tmp = {}
            tmp['nm'] = nm
            tmp['power'] = power
            tmp['lum'] = lum
            sdict[i] = tmp   

        else:
            for _ in range(int(1/(1/refreshrate))):
                circle.draw()
                win.flip()

    # save stuff
    if isinitmeasure:
        output = open(f'./measurements/data_init_{repeats}.pkl', 'wb')
        pickle.dump(sdict, output)
        output.close()
    else:
        output = open(f'./measurements/data_verification.pkl', 'wb')
        pickle.dump(sdict, output)
        output.close()

    win.close()


### RUN FUNCTION
if isinitmeasure:
    all_colours = pd.read_csv('colours_rgb_init.csv',header=None)
else:
    all_colours = pd.read_csv('RGB_stims.csv')


if isinitmeasure:
    for repeats in range(2):
        sdict = dict.fromkeys(range(len(all_colours)))
        run_calibration(sdict,repeats)
        print(sdict)
else:
    sdict = dict.fromkeys(range(len(all_colours)))
    run_calibration(sdict,repeats=0)






