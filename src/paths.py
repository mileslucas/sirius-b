import os


right_here = os.path.abspath(os.path.dirname(__file__))
base = os.path.join(right_here, "..", "data")
raw = os.path.join(base, "raw", "2020feb04")
flat = os.path.join(base, "flat")
sky = os.path.join(base, "sky")
sci = os.path.join(base, "sci")
output = os.path.join(base, "processed")
