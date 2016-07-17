#!/usr/bin/python

from math import pi, sin, cos
from subprocess import call
import sys

def safe_call(*args, **kwargs):
    if call(*args, **kwargs) != 0:
        print "ABORTING THETA ANIMATION"
        sys.exit()

def deg_to_rad(x):
    return x / 180. * pi

def animation_curve(x):
    """Given x from 0 to 1, return a smooth variation f(x) from 0 to 1,
    where f'(0) = f'(1) = 0."""
    return (sin(pi * (x - .5)) + 1) / 2

theta_zoom = 10
theta0 = 1 * theta_zoom
theta1 = 90 * theta_zoom
phi0 = 0
phi1 = 0

# WIDTH = 512
# HEIGHT = WIDTH
FOVY = 1.5
WIDTH = 1920
HEIGHT = WIDTH * 9 / 16
# FOVY = 0.75

FPS = 25
SECONDS = 8


DISTANCE = 1000 * .1

FRAME_COUNT = FPS * SECONDS

safe_call(["make"], cwd="../")
safe_call(["cp", "../otr", "otr_copy"])

# frames = [k for k in xrange(FRAME_COUNT) if k % 2 == 1]
frames = range(FRAME_COUNT)
approxbytes = WIDTH * HEIGHT * (12 + 4) * len(frames)
print "It will take approx. {}MB.".format((approxbytes + 1048575) / 1048576)

for k in frames:
    # theta = theta0 + (theta1 - theta0) * animation_curve(x)
    # theta = theta0 + (theta1 - theta0) * (sin(pi * (x - .5)) + 1) / 2
    # phi = phi0 + (phi1 - phi0) * (sin(pi * (x - .5)) + 1) / 2
    theta = theta0 + (theta1 - theta0) * float(k) / (FRAME_COUNT - 1)
    if theta < 1:
        theta = 1.0
    elif 90 * theta_zoom - 1 <= theta <= 90 * theta_zoom + 1:
        theta = 90.0 * theta_zoom - 1.0
    theta = deg_to_rad(theta / theta_zoom)
    phi = deg_to_rad(phi0 + (phi1 - phi0) * float(k) / (FRAME_COUNT - 1))

    camx = DISTANCE * sin(theta) * cos(phi)
    camy = DISTANCE * sin(theta) * sin(phi)
    camz = DISTANCE * cos(theta)

    # phi += pi / 2
    upx = -cos(theta) * cos(phi)
    upy = -cos(theta) * sin(phi)
    upz = sin(theta)

    filename_tga = "./output/animation/theta%05d.tga" % k
    filename_png = "./output/animation/theta%05d.png" % k
    filename_tgf = "./output/animation/float/theta%05d.tgf" % k
    args = [
        "./scripts/otr_copy",
        "--width", str(WIDTH),
        "--height", str(HEIGHT),
        "--fovy", str(FOVY),
        "--threads", "8",
        # "--ortho",
        # "--ver_range_km", "80",
        "--recursive",
        "--max_extra_recursive_depth", "1",
        "--camera_position_cart", "{},{},{}".format(camx, camy, camz),
        "--camera_up_cart", "{},{},{}".format(upx, upy, upz),
        "--no-cache",
        # "--output_float", filename_tgf,
        "--output_image", filename_tga,
        "--debug", "0",
        # "--horizontal_flip",
    ]
    print "Frame #{}/{}:   ".format(k + 1, FRAME_COUNT), " ".join(args)
    print "   <Total> theta deg: {}".format(theta * 180 / pi)
    safe_call(args, cwd="../")
    # safe_call(["convert", "." + filename_tga, "." + filename_png])
    # safe_call(["rm", "." + filename_tga])

if len(frames) >= FRAME_COUNT / 3:
    safe_call([
            "avconv",
            "-y",
            "-r", str(FPS),
            "-start_number", "1",
            "-i", "theta%05d.tga",
            "-b:v", "2000k",
            "test.mp4"
        ],
        cwd="../output/animation")
