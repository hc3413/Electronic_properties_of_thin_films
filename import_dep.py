# Standard Libraries
import os
import datetime
import pathlib
from pathlib import Path
import math
from dataclasses import dataclass, field
import readline

# Data Handling
import numpy as np
import pandas as pd

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.gridspec as gridspec
import seaborn as sns
from cycler import cycler
import scienceplots

# Scientific Computing
import scipy
from scipy.optimize import fsolve
from scipy.stats import linregress, zscore
import scipy.signal

# PowerPoint Automation
from pptx import Presentation
from pptx.util import Inches, Pt, Cm
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.oxml.ns import qn
from pptx.oxml import parse_xml







