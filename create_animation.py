
import analyze_r3d_functions as a3d
import matplotlib.pyplot as plt
import os

def create_movieinp():
    # Create scripts to create a 360 degree rotation
    # Create angle lists, the movie.inp file

    AUcm = 1.49598e13 # cm
    halfimagesize = 15*AUcm

    with open('../r3dsims/190_rotation/movie.inp', 'w') as fmovie:
        fmovie.write('1\n')
        fmovie.write('360\n')

        for angle in range(0,360):
            fmovie.write(f'0. 0. 0. {halfimagesize} {halfimagesize} 0. {angle}. 0.\n')


def create_pngfiles():
    # Load and create png-files of all images in 360-degree rotation
    # Apply new gamma-scaled imageplotter

    pathstardust =  '../r3dresults/st28gm06n052_staranddust_1/190movie_01um/'

    for nn in range(1,361):

        fig,ax,fluxtotal = a3d.plot_images(
            path=pathstardust,
            images=[f'image_{nn:04d}.out']
        )
        fig.savefig(f'{pathstardust}image_{nn:03d}.png', dpi=300, facecolor="white")


def create_giffile():
    os.system('cd ../r3dresults/st28gm06n052_staranddust_1/190movie_01um/')
    os.system("ffmpeg -framerate 60 -pattern_type glob -i '*.png' st28gm06n052_190_01um_rotate.gif")

